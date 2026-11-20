# -*- coding: utf-8 -*-
"""D5: posterior CH4 from D1 NetCDF + D4 regional TDI corrections.

Inputs:
    fluxes/gcp2025_flux_prior/gcp25_inv1.nc (348 months, 1998-2026)
    results2025/<CASE>/s057_p30/s057_ctl_inv1/predicted.sources
    inv_dir/mask/mask80.bin and ch4_regn_pulse.grd
Outputs:
    results2025/<CASE>/flux2d/p30/s057_ctl_inv1_noLC.nc
    results2025/<CASE>/flux2d/p30/s057_ctl_inv1_noLC.gt3
    results2025/<CASE>/diagnostics/d5/p30_ctl_inv1_noLC/

No prior/posterior GRD is created. A Fortran process receives monthly posterior
fields on stdin and performs the conservative 1-degree -> T42 GT3 conversion.
No LossCorr is calculated, read or applied.
"""
import argparse
import calendar
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.io import netcdf_file

from _set_case import SetCase


# ============================================================
# === D5: posterior flux assembly
class PostFlux(SetCase):

    def __init__(self, case=None):
        super().__init__()
        if case is not None and self.hcase != str(case).upper():
            raise ValueError('Case mismatch: ' + self.hcase + ' != ' + str(case))
        self.case = self.hcase
        self.root = Path(self.inv_wrk_dir).resolve()
        self.diag_root = self.root / 'diagnostics' / 'd5'
        self.source_code = Path(__file__).resolve().parent / 'fortran' / 'df5_post_flux.f90'
        override = os.environ.get('CH4_EXE_POST', '').strip()
        default = Path(self.fort_exe_dir) / 'df5_post_flux'
        self.executable = Path(override or default).expanduser().resolve()

    # --------------------------------------------------------
    # --- All existing project paths, no old stage_impl or prior GRD
    def get_paths(self, prior='p30', unc='ctl', model='inv1'):
        if prior not in self.unpcases or unc not in self.unxcases or model not in self.invcases:
            raise ValueError('Invalid D5 experiment: ' + str((prior, unc, model)))
        run = self.root / (self.icase + '_' + prior) / (self.icase + '_' + unc + '_' + model)
        stem = self.icase + '_' + unc + '_' + model + '_noLC'
        out = self.root / 'flux2d' / prior
        diag = self.diag_root / (prior + '_' + unc + '_' + model + '_noLC')
        return {
            'source': run / 'predicted.sources',
            'prior': Path(self.flx_apr_dir).resolve() / self.flx_apr_nc,
            'mask': Path(self.inv_dir).resolve() / 'mask' / 'mask80.bin',
            'basis': Path(self.inv_dir).resolve() / 'mask' / 'ch4_regn_pulse.grd',
            'nc': out / (stem + '.nc'),
            'gt3': out / (stem + '.gt3'),
            'diag': diag,
            'budget': diag / 'monthly_budget.csv',
            'regional': diag / 'regional_budget.csv',
        }

    @staticmethod
    def check_file(filename):
        path = Path(filename)
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError('D5 input missing/empty: ' + str(path))
        if path.is_symlink():
            raise ValueError('D5 input cannot be a symlink: ' + str(path))

    # --------------------------------------------------------
    # --- Read D4 TDI native coefficients (Tg CH4 per 365-day year)
    def read_sources(self, filename):
        tdi = np.empty((self.ntime, 54), dtype=np.float64)
        with Path(filename).open('r', encoding='ascii', errors='replace') as f:
            for r in range(54):
                label = f.readline().strip()
                if label != 'land_' + format(r + 1, '02d'):
                    raise ValueError('Invalid D4 region label: ' + repr(label))
                header = f.readline()
                if not header:
                    raise ValueError('Truncated D4 header: region ' + str(r + 1))
                for m in range(self.ntime):
                    row = f.readline().replace(',', ' ').split()
                    if len(row) < 2:
                        raise ValueError('Truncated D4 month: ' + str((r + 1, m + 1)))
                    stamp = float(row[0].replace('D', 'E').replace('d', 'e'))
                    value = float(row[1].replace('D', 'E').replace('d', 'e'))
                    year, month = self.yr_s + m // 12, m % 12 + 1
                    expect = year + (month - 0.5) / 12
                    if abs(stamp - expect) > 0.001:
                        raise ValueError('D4 calendar mismatch: region ' + str(r + 1) +
                                         ', expected ' + str((year, month)))
                    if not np.isfinite(value):
                        raise ValueError('Nonfinite TDI regional coefficient')
                    tdi[m, r] = value
            # --- Legacy D5 and D4 diagnostics consume exactly 54 land regions.
            #     TDI may append further records; they are not regional corrections.
            trailing = [line.strip() for line in f if line.strip()]
            if trailing:
                print('D5 TDI: ignored', len(trailing),
                      'trailing nonblank records after land_54; first:', repr(trailing[0][:120]))
        if np.max(np.abs(tdi)) < 1.0e-12:
            raise ValueError('ZERO predicted.sources: refusing posterior output')
        return tdi

    # --------------------------------------------------------
    # --- Read record 1 of legacy direct-access Fortran files.
    #     The original df5_ch4_grd.f90 uses READ(unit, REC=1), regardless
    #     of how many complete grid records the file contains.
    @staticmethod
    def read_direct_grid(filename):
        path = Path(filename)
        cells = 360 * 180
        record_bytes = cells * np.dtype('<f4').itemsize
        size = path.stat().st_size
        if size < record_bytes or size % record_bytes:
            raise ValueError('Invalid Fortran direct-access grid: ' + str(path) +
                             '; size=' + str(size) + ' bytes, expected a whole number of ' +
                             str(record_bytes) + '-byte records (360x180 float32).')
        with path.open('rb') as file:
            grid = np.fromfile(file, dtype='<f4', count=cells)
        if grid.size != cells:
            raise ValueError('Short first grid record: ' + str(path))
        print('D5 input grid:', path.name, '| bytes:', size,
              '| 1-degree records:', size // record_bytes, '| selected: 1')
        return grid.reshape(180, 360).astype(np.float64)

    # --------------------------------------------------------
    # --- Decode original direct-access mask and land basis, not GRD prior.
    def read_regions(self, paths, tdi):
        mask = self.read_direct_grid(paths['mask'])
        basis = self.read_direct_grid(paths['basis'])
        if not np.isfinite(mask).all() or not np.isfinite(basis).all():
            raise ValueError('Invalid mask or basis: NaN/Inf')
        if (basis < 0).any():
            raise ValueError('Negative regional pulse basis')
        index = np.rint(mask).astype(np.int16)
        valid = (index >= 1) & (index <= 54) & (np.abs(index - mask) < 0.01)
        if not np.any(valid):
            raise ValueError('Mask contains no recognized regions')
        # --- Match original df5_ch4_grd.f90: only 1..54 are inverted land regions.
        #     mask80.bin may contain other (e.g., 55..80) region IDs: leave their
        #     prior flux unchanged, just as the legacy Fortran does.
        ignored = (~valid) & (np.abs(mask) > 0.01)
        if np.any(ignored):
            values, counts_ignored = np.unique(mask[ignored], return_counts=True)
            order = np.argsort(counts_ignored)[::-1][:8]
            summary = ', '.join(str(float(values[i])) + ':' + str(int(counts_ignored[i])) for i in order)
            print('D5 regional mask: ignored', int(np.count_nonzero(ignored)),
                  'cells outside regions 1..54; most common values:', summary)
        region = np.where(valid, index, 0)
        lat_edges = np.deg2rad(np.arange(-90, 91, dtype=np.float64))
        area_lat = 6371000.0**2 * (2 * np.pi / 360) * np.diff(np.sin(lat_edges))
        area = area_lat[:, None]
        mass_year = basis * area * 365 * 86400 * 1.0e-9
        total = np.bincount(region.ravel(), weights=mass_year.ravel(), minlength=55)
        counts = np.bincount(region.ravel(), minlength=55)
        for r in range(1, 55):
            if total[r] <= 0 and np.max(np.abs(tdi[:, r - 1])) > 1.e-12:
                raise ValueError('Missing basis for nonzero TDI region ' + str(r))
        print('D5 regional mask cells:', int(np.sum(valid)), 'unassigned:', int(np.sum(~valid)))
        return region, basis, area_lat, total, counts

    # --------------------------------------------------------
    # --- Verify D1 coordinates + time before streaming any data to Fortran
    def check_prior_nc(self, ds):
        var = ds.get('flux_ch4_prior')
        if var is None or var.dims != ('time', 'lat', 'lon'):
            raise ValueError('D1 NetCDF must contain flux_ch4_prior(time,lat,lon)')
        if var.attrs.get('units') not in ('kg-CH4/m2/s', 'kg CH4 m-2 s-1', 'kg CH4/m2/s'):
            raise ValueError('Unexpected D1 prior units: ' + str(var.attrs.get('units')))
        if var.shape != (348, 180, 360):
            raise ValueError('Expected D1 prior shape (348,180,360), got ' + str(var.shape))
        expect_time = pd.date_range('1998-01-01', '2026-12-01', freq='MS')
        dates = pd.to_datetime(ds.time.values)
        if not np.array_equal(dates.values, expect_time.values):
            raise ValueError('D1 NetCDF has wrong monthly calendar; expected 1998-2026')
        lat = np.asarray(ds.lat.values, dtype=np.float64)
        lon = np.asarray(ds.lon.values, dtype=np.float64)
        ilat = np.argsort(lat)
        ilon = np.argsort((lon + 180) % 360 - 180)
        normlon = (lon[ilon] + 180) % 360 - 180
        if not np.allclose(lat[ilat], np.arange(-89.5, 90, 1), atol=1.e-4, rtol=0):
            raise ValueError('D1 prior latitude centers are not standard 1-degree grid')
        if not np.allclose(normlon, np.arange(-179.5, 180, 1), atol=1.e-4, rtol=0):
            raise ValueError('D1 prior longitude centers are not standard 1-degree grid')
        if (self.yr_s, self.yr_e, self.ntime) != (1999, 2025, 324):
            # Tests may override the period, but the input calendar must still be complete.
            if self.ntime < 1 or self.yr_s < 1998 or self.yr_s > 2026:
                raise ValueError('Invalid D5 interval')
        return ilat, ilon

    # --------------------------------------------------------
    # --- Compile only the Fortran 1-degree to T42 stream converter
    def build_fortran(self):
        if not self.source_code.is_file():
            raise FileNotFoundError('Missing D5 Fortran: ' + str(self.source_code))
        if os.environ.get('CH4_EXE_POST', '').strip():
            self.check_file(self.executable)
            return
        if self.executable.is_file() and self.executable.stat().st_mtime >= self.source_code.stat().st_mtime:
            return
        compiler = os.environ.get('FC', '').strip()
        if not compiler:
            compiler = 'gfortran' if shutil.which('gfortran') else 'ifort'
        if not shutil.which(compiler):
            raise RuntimeError('Fortran compiler unavailable: ' + compiler)
        flags = ['-O2', '-ffree-line-length-none'] if 'gfortran' in Path(compiler).name else ['-O2', '-free']
        self.executable.parent.mkdir(parents=True, exist_ok=True)
        cmd = [compiler] + flags + ['-o', str(self.executable), str(self.source_code)]
        print('D5 build:', ' '.join(cmd))
        subprocess.run(cmd, check=True, cwd=str(self.source_code.parent))

    # --------------------------------------------------------
    # --- Output NetCDF3 compatible with D1; no prior GRD anywhere
    def write_nc_header(self, outfile):
        nc = netcdf_file(str(outfile), mode='w', version=2)
        nc.createDimension('time', self.ntime)
        nc.createDimension('lat', 180)
        nc.createDimension('lon', 360)
        tim = nc.createVariable('time', 'i', ('time',))
        tim.units = ('days since ' + str(self.yr_s) + '-01-01 00:00:00').encode('ascii')
        tim.calendar = b'proleptic_gregorian'
        dates = pd.date_range(str(self.yr_s) + '-01-01', periods=self.ntime, freq='MS')
        tim[:] = np.array([(day - dates[0]).days for day in dates], dtype=np.int32)
        lat = nc.createVariable('lat', 'f', ('lat',))
        lon = nc.createVariable('lon', 'f', ('lon',))
        lat.units, lon.units = b'degrees_north', b'degrees_east'
        lat[:] = np.arange(-89.5, 90, 1, dtype=np.float32)
        lon[:] = np.arange(-179.5, 180, 1, dtype=np.float32)
        flux = nc.createVariable('flux_ch4_post', 'f', ('time', 'lat', 'lon'))
        flux.units = b'kg-CH4/m2/s'
        flux.long_name = b'CH4 posterior emissions without loss correction'
        nc.source_prior = b'gcp25_inv1.nc'
        nc.source_inversion = b'D4 predicted.sources (54 regional series)'
        nc.loss_correction = b'NOT APPLIED'
        nc.inversion_case = self.case.encode('ascii')
        return nc, flux

    @staticmethod
    def quote(value):
        value = str(value)
        if "'" in value or '\n' in value or '\r' in value or len(value) > 1024:
            raise ValueError('Invalid Fortran namelist path')
        return "'" + value + "'"

    def write_config(self, paths, staged):
        contents = [
            '&post_config',
            ' out_gt3 = ' + self.quote(staged['gt3']),
            ' budget_file = ' + self.quote(staged['gt3_budget']),
            ' initial_year = ' + str(self.yr_s),
            ' nmonths = ' + str(self.ntime),
            '/',
        ]
        staged['config'].write_text('\n'.join(contents) + '\n', encoding='ascii')

    # --------------------------------------------------------
    # --- Annual APR + CORR = POST budget, following the legacy check_gtotal.
    #     Monthly values are actual calendar-month masses (Tg CH4/month).
    #     Annual totals are sums over all months (Tg CH4/calendar year).
    def diagnostics(self, paths):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        budget = pd.read_csv(paths['budget'])
        regions = pd.read_csv(paths['regional'])
        if len(budget) != self.ntime or len(regions) != 54 * self.ntime:
            raise ValueError('D5 incomplete diagnostic budgets')
        if not np.isfinite(budget.select_dtypes(include='number').values).all():
            raise ValueError('Nonfinite D5 budget')
        if not np.allclose(budget.post_Tg_month, budget.prior_Tg_month + budget.correction_Tg_month,
                           rtol=3.e-7, atol=1.e-7):
            raise ValueError('Posterior monthly mass mismatch')
        if not np.allclose(budget.gt3_Tg_month, budget.post_Tg_month, rtol=3.e-7, atol=1.e-7):
            raise ValueError('GT3 monthly mass mismatch')
        if not np.allclose(regions.tdi_Tg_per_365d_year, regions.gridded_Tg_per_365d_year,
                           rtol=3.e-7, atol=2.e-8):
            raise ValueError('D5 regional budget mismatch')

        # --- Keep the original monthly budget. Annual values are sums, not December snapshots.
        month_counts = budget.groupby('year')['month'].nunique()
        annual = budget.groupby('year', as_index=False).agg(
            APR_Tg_year=('prior_Tg_month', 'sum'),
            CORR_Tg_year=('correction_Tg_month', 'sum'),
            POST_Tg_year=('post_Tg_month', 'sum'),
            GT3_Tg_year=('gt3_Tg_month', 'sum'),
            negative_cells_max=('negative_cells', 'max'),
        )
        annual['balance_Tg_year'] = annual.POST_Tg_year - annual.APR_Tg_year - annual.CORR_Tg_year
        annual.to_csv(paths['diag'] / 'annual_budget.csv', index=False)
        if not np.allclose(annual.APR_Tg_year + annual.CORR_Tg_year, annual.POST_Tg_year,
                           rtol=3.e-7, atol=1.e-7):
            raise ValueError('Posterior annual mass mismatch')

        core = regions[(regions.year >= 2000) & (regions.year <= 2024)]
        stats = core.groupby('region')['tdi_Tg_per_365d_year'].agg(['mean', 'std']).reset_index()
        stats.to_csv(paths['diag'] / 'regional_2000_2024.csv', index=False)
        top = stats.assign(absm=lambda x: x['mean'].abs()).nlargest(12, 'absm')

        fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
        axes[0, 0].plot(annual.year, annual.APR_Tg_year, label='APR')
        axes[0, 0].plot(annual.year, annual.POST_Tg_year, label='POST')
        axes[0, 0].set_ylabel('Tg CH4 / calendar year')
        axes[0, 0].set_title('Annual APR and POST')
        axes[0, 0].legend()
        axes[0, 1].plot(annual.year, annual.CORR_Tg_year, label='CORR')
        axes[0, 1].axhline(0, color='gray', lw=0.8)
        axes[0, 1].set_ylabel('Tg CH4 / calendar year')
        axes[0, 1].set_title('Annual D4 correction: POST - APR')
        dates = pd.to_datetime(dict(year=budget.year, month=budget.month, day=1))
        axes[1, 0].plot(dates, budget.prior_Tg_month, label='APR', alpha=0.7)
        axes[1, 0].plot(dates, budget.post_Tg_month, label='POST')
        axes[1, 0].set_ylabel('Tg CH4 / calendar month')
        axes[1, 0].legend()
        axes[1, 1].barh(top.region.astype(str), top['mean'])
        axes[1, 1].axvline(0, color='gray', lw=0.8)
        axes[1, 1].set_xlabel('Tg CH4 / 365-day year')
        axes[1, 1].set_title('Top regional signed corrections')
        for ax in axes.flat:
            ax.grid(alpha=0.2)
        fig.suptitle('D5 POST: ' + self.case + ' / no LossCorr / ' + self.icase)
        fig.savefig(paths['diag'] / 'post_overview.png', dpi=140)
        plt.close(fig)

        # --- Print the same three quantities as the original df5_ch4_grd check_gtotal.
        table = [
            'D5 annual global CH4 fluxes [' + self.case + ', Tg CH4/calendar year]',
            '  Year         APR         CORR         POST          GT3',
        ]
        for row in annual.itertuples(index=False):
            table.append(f'  {int(row.year):4d} {row.APR_Tg_year:11.3f} {row.CORR_Tg_year:+12.3f} {row.POST_Tg_year:12.3f} {row.GT3_Tg_year:12.3f}')
        complete_years = month_counts.index[month_counts == 12]
        core_annual = annual[(annual.year >= 2000) & (annual.year <= 2024) & annual.year.isin(complete_years)]
        if len(core_annual):
            mean_apr = float(core_annual.APR_Tg_year.mean())
            mean_corr = float(core_annual.CORR_Tg_year.mean())
            mean_post = float(core_annual.POST_Tg_year.mean())
            table.append('')
            table.append('Core 2000-2024 (complete years only), mean Tg CH4/year:')
            table.append(f'APR = {mean_apr:.4f}   CORR = {mean_corr:+.4f}   POST = {mean_post:.4f}')
            table.append(f'Check APR + CORR - POST = {mean_apr + mean_corr - mean_post:+.3e} Tg/year')
        else:
            table.append('Core 2000-2024: no complete calendar years in this test window')
        table.append('APR + CORR = POST; LossCorr NOT APPLIED')
        print('\n'.join(table))

        lines = [
            'D5 POST / ' + self.case + ' / ' + self.icase,
            'Loss correction: NOT APPLIED',
            'Prior: ' + str(paths['prior']),
            'D4 sources: ' + str(paths['source']),
            'Posterior: ' + str(paths['nc']),
            'Period: ' + str(self.yr_s) + '-' + str(self.yr_s + (self.ntime - 1) // 12),
            '',
            *table,
            '',
            'Max relative NC->GT3 budget error: ' + str(float(np.max(np.abs(budget.gt3_rel_error)))),
            'Max monthly negative cells: ' + str(int(budget.negative_cells.max())),
            'Global/regional monthly budgets: PASS',
            '1999 is inversion boundary, excluded from core means',
            'D4 native TDI coefficients: Tg CH4 per 365-day year',
            'Monthly budget columns: Tg CH4/calendar month',
            'Annual totals: Tg CH4/calendar year',
            'NC and GT3 fields: kg CH4 m-2 s-1',
        ]
        (paths['diag'] / 'diagnostic_report.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')
        print('D5 max NC->GT3 relative budget error:', float(np.max(np.abs(budget.gt3_rel_error))))
        print('D5 max monthly negative cells:', int(budget.negative_cells.max()))
        print('D5 global/regional monthly budgets: PASS')

    # --------------------------------------------------------
    # --- Stream 324 posterior months: NetCDF -> Python -> Fortran stdin -> GT3
    def run(self, prior='p30', unc='ctl', model='inv1', dry_run=False):
        paths = self.get_paths(prior, unc, model)
        for key in ('source', 'prior', 'mask', 'basis'):
            self.check_file(paths[key])
        print('D5 POST:', self.case, self.icase, prior, unc, model)
        print('D1 prior NetCDF:', paths['prior'])
        print('D4 sources:', paths['source'])
        print('Posterior:', paths['nc'])
        print('GT3:', paths['gt3'])
        print('LossCorr: NOT APPLIED')
        if dry_run:
            return
        if sys.byteorder != 'little':
            raise RuntimeError('D5 Fortran streaming expects little-endian host')
        tdi = self.read_sources(paths['source'])
        region, basis, area_lat, totals, counts = self.read_regions(paths, tdi)
        self.build_fortran()
        paths['nc'].parent.mkdir(parents=True, exist_ok=True)
        paths['diag'].mkdir(parents=True, exist_ok=True)
        flux_scale = 365 * 86400 * 1.e-9
        with xr.open_dataset(paths['prior'], decode_times=True) as prior_ds:
            ilat, ilon = self.check_prior_nc(prior_ds)
            # --- Do not modify source prior or D4 files
            with tempfile.TemporaryDirectory(prefix='_stage_', dir=paths['diag']) as td:
                td = Path(td)
                staged = {
                    'nc': td / 'posterior.nc',
                    'gt3': td / 'posterior.gt3',
                    'gt3_budget': td / 'gt3_budget.csv',
                    'config': td / 'post_config.nml',
                }
                self.write_config(paths, staged)
                logfile = paths['diag'] / 'fortran.log'
                rows = []
                reg_rows = []
                with logfile.open('w', encoding='utf-8') as log:
                    nc, nc_flux = self.write_nc_header(staged['nc'])
                    try:
                        process = subprocess.Popen([str(self.executable), str(staged['config'])],
                                                   stdin=subprocess.PIPE, stdout=log, stderr=subprocess.STDOUT,
                                                   cwd=str(td))
                        try:
                            for m in range(self.ntime):
                                year, month = self.yr_s + m // 12, m % 12 + 1
                                prior_index = (year - 1998) * 12 + month - 1
                                field = np.asarray(prior_ds['flux_ch4_prior'].isel(time=prior_index).values,
                                                   dtype=np.float64)
                                field = field[np.ix_(ilat, ilon)]
                                if not np.isfinite(field).all():
                                    raise ValueError('Invalid D1 prior flux: ' + str((year, month)))
                                factor = np.zeros(55, dtype=np.float64)
                                good = totals[1:] > 0
                                factor[1:][good] = tdi[m, good] / totals[1:][good]
                                correction = basis * factor[region]
                                post64 = field + correction
                                post = np.asarray(post64, dtype=np.float32, order='C')
                                if not np.isfinite(post).all():
                                    raise ValueError('Nonfinite posterior flux')
                                seconds = calendar.monthrange(year, month)[1] * 86400
                                scale = area_lat[:, None] * seconds * 1.e-9
                                prior_mass = float(np.sum(field * scale))
                                post_mass = float(np.sum(post.astype(np.float64) * scale))
                                correction_mass = float(np.sum(tdi[m, :]) * seconds / (365 * 86400))
                                actual_corr = float(np.sum(correction * scale))
                                abs_total = float(np.sum(np.abs(post) * scale))
                                tol = 3.e-7 * max(abs_total, 1.e-6)
                                if abs(prior_mass + correction_mass - post_mass) > tol:
                                    raise RuntimeError('D5 mass budget mismatch: ' + str((year, month)))
                                if abs(actual_corr - correction_mass) > 3.e-8 * max(1., abs(correction_mass)):
                                    raise RuntimeError('D5 gridded correction mismatch')
                                region_mass = np.bincount(region.ravel(),
                                      weights=(correction * area_lat[:, None] * flux_scale).ravel(), minlength=55)
                                if not np.allclose(region_mass[1:], tdi[m, :], rtol=3.e-7, atol=2.e-8):
                                    raise RuntimeError('D5 region balance mismatch: ' + str((year, month)))
                                for r in range(1, 55):
                                    reg_rows.append((year, month, r, tdi[m, r - 1], region_mass[r], counts[r]))
                                rows.append((year, month, prior_mass, correction_mass, post_mass,
                                             int((post < 0).sum()), float(post.min()), float(post.max())))
                                nc_flux[m, :, :] = post
                                process.stdin.write(post.tobytes(order='C'))
                            process.stdin.close()
                            if process.wait() != 0:
                                raise RuntimeError('D5 Fortran converter failed; see ' + str(logfile))
                        except Exception:
                            if process.poll() is None:
                                process.kill()
                            process.wait()
                            raise
                    finally:
                        nc.close()
                if staged['gt3'].stat().st_size != self.ntime * 33808:
                    raise RuntimeError('D5 GT3 unexpected size')
                nc_verify = xr.open_dataset(staged['nc'])
                try:
                    if nc_verify['flux_ch4_post'].shape != (self.ntime, 180, 360):
                        raise RuntimeError('D5 NetCDF output wrong shape')
                finally:
                    nc_verify.close()
                budget = pd.DataFrame(rows, columns=['year', 'month', 'prior_Tg_month',
                            'correction_Tg_month', 'post_Tg_month', 'negative_cells', 'min_post', 'max_post'])
                fbudget = pd.read_csv(staged['gt3_budget'])
                if not np.array_equal(budget[['year', 'month']].values, fbudget[['year', 'month']].values):
                    raise RuntimeError('D5 GT3 calendar mismatch')
                if not np.allclose(fbudget.input_Tg_month, budget.post_Tg_month, rtol=3.e-7, atol=1.e-7):
                    raise RuntimeError('D5 NC/GT3 input budget mismatch')
                budget['gt3_Tg_month'] = fbudget['gt3_Tg_month']
                budget['gt3_rel_error'] = fbudget['gt3_rel_error']
                region_df = pd.DataFrame(reg_rows, columns=['year', 'month', 'region',
                         'tdi_Tg_per_365d_year', 'gridded_Tg_per_365d_year', 'land_grid_cells'])
                budget.to_csv(paths['budget'], index=False)
                region_df.to_csv(paths['regional'], index=False)
                self.diagnostics(paths)
                staged['nc'].replace(paths['nc'])
                staged['gt3'].replace(paths['gt3'])
        print('D5 POST: PASS')


# ============================================================
# === _main.py entry point

def run_stage(case='CYC', dry_run=False):
    old = os.environ.get('GCP_CASE')
    os.environ['GCP_CASE'] = str(case).upper()
    try:
        PostFlux(case=case).run(dry_run=dry_run)
    finally:
        if old is None:
            os.environ.pop('GCP_CASE', None)
        else:
            os.environ['GCP_CASE'] = old


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True, choices=['CYC', 'INCA'])
    parser.add_argument('--prior', default='p30', choices=['p30', 'p50', 'p99'])
    parser.add_argument('--unc', default='ctl', choices=['ctl', 'ux2', 'ux4'])
    parser.add_argument('--model', default='inv1')
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--diagnostics-only', action='store_true', help='Reprint APR/CORR/POST from existing D5 budget CSV, without rebuilding outputs')
    args = parser.parse_args()
    if args.diagnostics_only and args.dry_run:
        parser.error('--diagnostics-only cannot be combined with --dry-run')
    os.environ['GCP_CASE'] = args.case
    obj = PostFlux(case=args.case)
    if args.diagnostics_only:
        paths = obj.get_paths(prior=args.prior, unc=args.unc, model=args.model)
        obj.check_file(paths['budget'])
        obj.check_file(paths['regional'])
        obj.diagnostics(paths)
    else:
        obj.run(prior=args.prior, unc=args.unc, model=args.model, dry_run=args.dry_run)
