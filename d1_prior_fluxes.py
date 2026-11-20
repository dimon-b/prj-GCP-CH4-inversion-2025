# -*- coding: utf-8 -*-
"""
D1: GCP/GMB2025 prior CH4 emissions.

All Python operations of this stage live here; no stage_impl or
_fortran_stage_helpers imports.

Workflow:
  1. GCP_Prior_CH4_fluxes.nc -> gcp25_inv1.nc / gcp25_inv1_full.nc
  2. Regrid (df1_regrid_ch4_t42) -> prior GT3 (T42)
  3. Regions (df1_aggregate_ch4_regions) -> regional prior CH4 fluxes
  4. Uncertainty (df1_write_prior_uncertainty) -> prior uncertainty NetCDF (p30/p50/p99)

Run from c_gcpv3:
  python _main.py --case INCA --steps d1_prior
  python _main.py --case INCA --steps d1_prior --dry-run
  python d1_prior_fluxes.py --case INCA [--dry-run]

The prior inputs and D1 Fortran results are shared between CYC and INCA.
Existing outputs are checked and reused, never silently overwritten.
"""

import argparse
import importlib.util
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

import _set_case


# ============================================================
# === D1 prior fluxes
class PriorFluxes(_set_case.SetCase):

    def __init__(self, case=None):
        super().__init__()
        if case is not None and self.hcase.upper() != case.upper():
            raise ValueError('SetCase case mismatch: ' + self.hcase + ' != ' + case +
                             '. Check GCP_CASE / CH4_CASE in _set_case.py.')
        self.case = self.hcase.upper()
        # --- D1 production time span; source GCP data remain 2000-2024.
        self.yr_s = self.d1_year_s
        self.yr_e = self.d1_year_e
        if list(self.years_nc) != [2000, 2024]:
            raise ValueError('D1 expects original GCP years_nc = [2000, 2024]')
        self.config_dir = Path(self.inv_wrk_dir).expanduser().resolve() / 'config'

    # --------------------------------------------------------
    # --- Prepare prior NetCDF files
    def get_inp_flx(self):
        inp = Path(self.flx_inp_dir) / self.flx_inp_nc
        prior = Path(self.flx_apr_dir) / self.flx_apr_nc
        full = Path(self.flx_apr_dir) / self.flx_apr_nc_full

        print('\nD1: Prepare prior NetCDF')
        print('Input:', inp)
        print('Prior:', prior)
        print('Full :', full)

        existing = [prior.exists(), full.exists()]
        if all(existing):
            self.check_prior_nc(prior)
            self.check_prior_nc(full)
            if self.is_fortran_nc(prior) and self.is_fortran_nc(full):
                print('Prior NetCDF: reuse verified NetCDF3 64-bit files')
                return
            # --- Previous h5netcdf output uses NC_STRING attributes.
            # --- NetCDF-Fortran 4.5.2 cannot read those as CHARACTER.
            if not inp.is_file():
                raise FileNotFoundError('Cannot rebuild incompatible prior NetCDF: ' + str(inp))
            backups = [p.with_name(p.name + '.netcdf4_backup') for p in (prior, full)]
            if any(p.exists() for p in backups):
                raise FileExistsError('NetCDF4 backup already exists; inspect prior outputs: ' + str(backups))
            print('Prior NetCDF: old NetCDF4 encoding is incompatible with D1 regridding')
            print('Archiving both old files and rebuilding as NetCDF3 64-bit offset')
            prior.rename(backups[0])
            try:
                full.rename(backups[1])
            except Exception:
                backups[0].rename(prior)
                raise
            existing = [False, False]
        if any(existing):
            raise RuntimeError('Incomplete prior NetCDF pair: ' + str(prior) + ' / ' + str(full))
        if not inp.is_file():
            raise FileNotFoundError(inp)

        # --- Preserve the existing GCP input reading method.
        input_engine = 'h5netcdf' if importlib.util.find_spec('h5netcdf') else None
        with xr.open_dataset(inp, engine=input_engine, decode_times=False) as ds:
            flux = self.manage_ds(ds)
            self.check_total(flux)
            self.write_nc(flux, prior, full)

        self.check_prior_nc(prior)
        self.check_prior_nc(full)

    # --------------------------------------------------------
    # --- Expand seasonal sectors and compute total prior flux
    def manage_ds(self, ds):
        if 'time' not in ds.sizes:
            raise ValueError('GCP input has no time dimension')

        ntime = ds.sizes['time']
        nseason = ds.sizes.get('time_climato', 0)
        expanded = {}

        # --- Decode original CF time, preserving its calendar.
        decoded_time = xr.decode_cf(ds[['time']])['time']
        months = pd.to_datetime(decoded_time.values)
        target = pd.date_range(str(self.years_nc[0]) + '-01-01',
                               str(self.years_nc[1]) + '-12-01', freq='MS')
        if (len(months) != len(target) or
                not np.array_equal(months.to_period('M').values, target.to_period('M').values)):
            raise ValueError('Unexpected GCP source monthly coverage; expected ' +
                             str(target[0]) + ' .. ' + str(target[-1]))
        # --- Mid-month timestamps in the input are valid monthly observations.
        # --- Normalize to first-of-month for NetCDF/Fortran consistency.
        decoded_time = xr.DataArray(target.values, dims=('time',), name='time')

        if nseason:
            if nseason != 12 or ntime % nseason:
                raise ValueError('Expected 12 time_climato months and full years')
            month_index = xr.DataArray(np.arange(ntime) % nseason, dims='time')
            for v in ds.data_vars:
                if 'time_climato' in ds[v].dims:
                    if ds[v].dims != ('time_climato', 'lat', 'lon'):
                        raise ValueError('Unsupported climatology dimensions for ' + v + ': ' + str(ds[v].dims))
                    # --- Same Jan..Dec repetition as the original np.tile logic.
                    expanded[v] = ds[v].isel(time_climato=month_index).drop_vars('time_climato', errors='ignore')
                    expanded[v] = expanded[v].assign_coords(time=decoded_time)

        time_vars = {}
        for v in ds.data_vars:
            if 'time' in ds[v].dims and 'time_climato' not in ds[v].dims:
                time_vars[v] = ds[v].assign_coords(time=decoded_time)

        if not expanded and not time_vars:
            raise ValueError('No GCP fluxes on time/time_climato dimensions')

        ds_out = xr.Dataset({**time_vars, **expanded}).assign_coords(time=decoded_time)
        flux_vars = [v for v in ds_out.data_vars
                     if v.startswith('flux_ch4_') and
                     v not in ('flux_ch4_prior', 'flux_ch4_prior_soil0')]
        if 'flux_ch4_soils' not in flux_vars or not flux_vars:
            raise ValueError('Missing required flux_ch4_soils or sector flux variables')

        for var in flux_vars:
            if ds_out[var].dims != ('time', 'lat', 'lon'):
                raise ValueError('Unexpected dimensions: ' + var + ' ' + str(ds_out[var].dims))
            print(var, 'NaN cells:', int(ds_out[var].isnull().sum().item()))

        # --- Missing sector fluxes other than explicitly masked soil uptake are fatal.
        for var in flux_vars:
            if var != 'flux_ch4_soils':
                if not bool(np.isfinite(ds_out[var]).all()):
                    raise ValueError('Missing/non-finite sector flux: ' + var)
        # --- Preserve soil NaN -> zero and original sum of sector fluxes.
        ds_out['flux_ch4_soils'] = ds_out['flux_ch4_soils'].fillna(0.0)
        ds_out['flux_ch4_prior'] = sum(ds_out[v] for v in flux_vars)
        ds_out['flux_ch4_prior_soil0'] = ds_out['flux_ch4_prior'] - ds_out['flux_ch4_soils']
        self.check_sector_balance(ds_out, flux_vars)

        return ds_out

    # --------------------------------------------------------
    # --- Sector-by-sector pointwise conservation, including negative soil sink
    def check_sector_balance(self, ds, flux_vars):
        total = sum(ds[v].astype('float64') for v in flux_vars)
        actual = ds['flux_ch4_prior'].astype('float64')
        residual = (total - actual).values
        if not np.isfinite(residual).all():
            raise ValueError('Non-finite total prior flux')
        scale = max(float(np.abs(total.values).max()), 1e-20)
        if float(np.abs(residual).max()) > 1e-6 * scale:
            raise ValueError('Prior flux differs from the sector sum')
        soil0 = ds['flux_ch4_prior_soil0'].astype('float64')
        residual = (soil0 - (actual - ds['flux_ch4_soils'])).values
        if float(np.abs(residual).max()) > 1e-6 * scale:
            raise ValueError('Prior soil0 differs from prior minus soil sink')
        print('Sector balance PASS, negative soil uptake retained')

    # --------------------------------------------------------
    # --- Original global flux integration diagnostic, Tg CH4/year
    def check_total(self, ds):
        dlat = np.deg2rad(ds.lat[1] - ds.lat[0])
        dlon = np.deg2rad(ds.lon[1] - ds.lon[0])
        lat_rad = np.deg2rad(ds.lat)
        cell_area_lat = xr.DataArray(self.R**2 * dlat * dlon * np.cos(lat_rad),
                                     coords={'lat': ds.lat}, dims=['lat'])
        cell_area = cell_area_lat.broadcast_like(ds['flux_ch4_prior'])
        earth_area = cell_area.isel(time=0).sum(dim=('lat', 'lon')).item()
        print('\nCheck total Earth area [m2]:', f'{earth_area:.5e}')
        if not np.isfinite(earth_area) or abs(earth_area - 4*np.pi*self.R**2) / (4*np.pi*self.R**2) > 1e-4:
            raise ValueError('Invalid integrated Earth area; inspect input lat/lon grid')

        time = pd.to_datetime(ds.time.values)
        sec = xr.DataArray(time.days_in_month * 24 * 3600,
                           coords={'time': ds.time}, dims=['time'])
        records = []
        for var in [v for v in ds.data_vars if v.startswith('flux_ch4_')]:
            monthly = ds[var] * cell_area * sec
            annual = monthly.sum(dim=('lat', 'lon')).groupby('time.year').sum(dim='time') / 1e9
            for year, value in zip(annual.year.values, annual.values):
                records.append({'year': int(year), 'variable': var.replace('flux_ch4_', ''),
                                'Tg_CH4_yr': float(value)})

        table = pd.DataFrame(records).pivot(index='year', columns='variable', values='Tg_CH4_yr')
        total_cols = [col for col in table if col.lower().startswith('prior')]
        table = table[[col for col in table if col not in total_cols] + total_cols]
        print(table.round(4).to_string())
        return table

    # --------------------------------------------------------
    # --- Extend original 2000-2024 fields into 1998-2026 and save NetCDF
    def extend_months(self, ds):
        start = self.yr_s
        stop = self.yr_e
        src_start, src_stop = map(int, self.years_nc)
        if (start, stop, src_start, src_stop) != (1998, 2026, 2000, 2024):
            raise ValueError('Unexpected D1 source/output years')

        source_time = pd.date_range('2000-01-01', '2024-12-01', freq='MS')
        if not np.array_equal(pd.to_datetime(ds.time.values).values, source_time.values):
            raise ValueError('Unexpected original monthly time axis')
        parts = []
        for year in [1998, 1999]:
            part = ds.sel(time=slice('2000-01', '2000-12')).copy(deep=True)
            parts.append(part.assign_coords(time=pd.date_range(f'{year}-01-01', periods=12, freq='MS')))
        parts.append(ds)
        for year in [2025, 2026]:
            part = ds.sel(time=slice('2024-01', '2024-12')).copy(deep=True)
            parts.append(part.assign_coords(time=pd.date_range(f'{year}-01-01', periods=12, freq='MS')))
        out = xr.concat(parts, dim='time')
        self.check_extended_ds(out)
        return out

    def check_extended_ds(self, out):
        expected = pd.date_range('1998-01-01', '2026-12-01', freq='MS')
        actual = pd.to_datetime(out.time.values)
        if out.sizes['time'] != 348 or not np.array_equal(actual.values, expected.values):
            raise ValueError('Wrong extended time axis; expected 348 months (1998-2026)')
        variables = [v for v in out if v.startswith('flux_ch4_')]
        for target, origin in [(1998, 2000), (1999, 2000), (2025, 2024), (2026, 2024)]:
            for var in variables:
                a = out[var].sel(time=str(target)).values
                b = out[var].sel(time=str(origin)).values
                if not np.array_equal(a, b, equal_nan=True):
                    raise ValueError(f'Edge copy differs: {var}, {target} != {origin}')
        print('Edge-copy check PASS: 1998/1999=2000; 2025/2026=2024')

    def write_nc(self, ds, prior, full):
        out = self.extend_months(ds)
        print('\nOutput months:', out.sizes['time'])
        # --- Keep output compatible with netCDF-Fortran 4.5.2 on the server.
        # --- Do not use h5netcdf: its NC_STRING attributes break nf90_get_att.
        out.attrs = {}
        for var in out.variables:
            out[var].attrs = {}
            out[var].encoding = {}

        prior.parent.mkdir(parents=True, exist_ok=True)
        if prior.exists() or full.exists():
            raise FileExistsError('Prior NetCDF outputs appeared during conversion')

        # --- Old netCDF-Fortran 4.5.2 requires NC_CHAR attributes, not NC_STRING.
        # --- NetCDF3 64-bit offset is readable by old Fortran and supports >2 GB files.
        engine = 'scipy'
        encoding = {'format': 'NETCDF3_64BIT'}
        out['flux_ch4_prior'].attrs['units'] = 'kg-CH4/m2/s'
        for var in out.data_vars:
            if var.startswith('flux_ch4_'):
                out[var].attrs['units'] = 'kg-CH4/m2/s'

        # --- Only write each output if both are free. If the second write
        # --- fails, do not leave behind a valid-looking incomplete pair.
        try:
            print('Write prior:', prior)
            out['flux_ch4_prior'].fillna(0.0).to_netcdf(prior, engine=engine, mode='w', **encoding)
            print('Write all sectors:', full)
            out.fillna(0.0).to_netcdf(full, engine=engine, mode='w', **encoding)
        except Exception:
            prior.unlink(missing_ok=True)
            full.unlink(missing_ok=True)
            raise

    # --------------------------------------------------------
    # --- NetCDF3 64-bit offset: numeric time and NC_CHAR attributes for Fortran.
    @staticmethod
    def is_fortran_nc(filename):
        with Path(filename).open('rb') as f:
            return f.read(4) == b'CDF\x02'

    # --------------------------------------------------------
    # --- Verify NetCDF monthly coverage and flux content
    def check_prior_nc(self, filename):
        if not filename.is_file() or filename.stat().st_size == 0:
            raise FileNotFoundError(filename)
        with xr.open_dataset(filename, decode_times=True) as ds:
            if 'flux_ch4_prior' not in ds:
                raise ValueError('Missing flux_ch4_prior in ' + str(filename))
            years = self.yr_e - self.yr_s + 1
            if ds.sizes.get('time') != 12*years:
                raise ValueError(f'Legacy/incompatible prior NetCDF {filename}: {ds.sizes.get("time")} months, expected {12*years} (1998-2026). Archive BOTH old prior NetCDF outputs before rebuilding D1; no files were changed.')
            if ds.sizes.get('lat') != 180 or ds.sizes.get('lon') != 360:
                raise ValueError('Wrong 1-degree grid in ' + str(filename))
            tm = pd.to_datetime(ds.time.values)
            if (tm[0].year, tm[0].month, tm[-1].year, tm[-1].month) != (self.yr_s, 1, self.yr_e, 12):
                raise ValueError('Wrong output time extent in ' + str(filename))
            expected = pd.date_range('1998-01-01', '2026-12-01', freq='MS')
            if not np.array_equal(tm.values, expected.values):
                raise ValueError('Missing, repeated or unordered monthly records')
            prior = ds['flux_ch4_prior']
            if prior.dims != ('time', 'lat', 'lon') or not bool(np.isfinite(prior).all()):
                raise ValueError('Invalid prior flux shape or non-finite values')
            self.check_extended_ds(ds)
            if 'flux_ch4_soils' in ds:
                flux_vars = [v for v in ds if v.startswith('flux_ch4_')
                             and v not in ('flux_ch4_prior', 'flux_ch4_prior_soil0')]
                self.check_sector_balance(ds, flux_vars)
            print('NetCDF check PASS:', filename.name, ds.sizes['time'], 'months')

    # --------------------------------------------------------
    # --- GT3 complete structure and header calendar (big endian Fortran sequential)
    def check_gt3(self, filename):
        import struct
        nrec = 348
        data_bytes = 128*64*4
        rec_bytes = 8 + 64*16 + 8 + data_bytes
        if not filename.is_file() or filename.stat().st_size != nrec*rec_bytes:
            raise ValueError('GT3 expected 348 complete surface records: ' + str(filename))
        with filename.open('rb') as f:
            for i in range(nrec):
                head_len = struct.unpack('>i', f.read(4))[0]
                if head_len != 1024:
                    raise ValueError(f'GT3 invalid header size at {i}')
                head = f.read(1024)
                if struct.unpack('>i', f.read(4))[0] != 1024:
                    raise ValueError(f'GT3 bad header trailer at {i}')
                words = [head[j:j+16].decode('ascii').strip() for j in range(0, 1024, 16)]
                year, month = 1998+i//12, 1+i%12
                if (words[2], words[15], words[37]) != ('QFLUX07', 'kg-CH4/m2/s', 'UR4'):
                    raise ValueError(f'GT3 invalid variable/units/format at {i}')
                if words[26][:6] != f'{year:04d}{month:02d}' or int(words[24]) != i:
                    raise ValueError(f'GT3 incorrect time header at {i}')
                if int(words[63]) != 128*64:
                    raise ValueError(f'GT3 invalid grid count at {i}')
                if struct.unpack('>i', f.read(4))[0] != data_bytes:
                    raise ValueError(f'GT3 invalid data length at {i}')
                vals = np.frombuffer(f.read(data_bytes), dtype='>f4')
                if len(vals) != 128*64 or not np.isfinite(vals).all():
                    raise ValueError(f'GT3 invalid flux values at {i}')
                if struct.unpack('>i', f.read(4))[0] != data_bytes:
                    raise ValueError(f'GT3 invalid data trailer at {i}')
            if f.read(1):
                raise ValueError('GT3 trailing bytes')
        print('GT3 check PASS:', filename.name, nrec, 'months')

    # --------------------------------------------------------
    # --- Compare all 348 monthly global integrals on the 1-degree and T42 grids.
    # --- Negative soil uptake is preserved: use net, not clipped, emissions.
    def check_flux_conservation(self, input_nc, output_gt3):
        radius = float(self.R)
        from numpy.polynomial.legendre import leggauss
        _, gw = leggauss(64)
        t42_area = (2*np.pi*radius**2/128) * gw[::-1]
        header_bytes = 1024
        data_bytes = 128*64*4
        with xr.open_dataset(input_nc, decode_times=True) as ds:
            lat = ds['lat'].values.astype(float)
            lon = ds['lon'].values.astype(float)
            if len(lat) != 180 or len(lon) != 360:
                raise ValueError('Conservation test expects 1-degree 360x180 grid')
            if not np.allclose(np.abs(np.diff(lat)), 1., atol=1e-5):
                raise ValueError('Non-1-degree input latitudes')
            if not np.allclose(np.sort(np.mod(lon, 360)), np.arange(360)+0.5, atol=1e-4) and \
               not np.allclose(np.sort(np.mod(lon, 360)), np.arange(360), atol=1e-4):
                raise ValueError('Unexpected 1-degree input longitude centers')
            area = radius**2 * (np.pi/180.) * (
                np.sin(np.deg2rad(lat+0.5)) - np.sin(np.deg2rad(lat-0.5)))
            if abs(area.sum()*360-4*np.pi*radius**2)/(4*np.pi*radius**2) > 1e-5:
                raise ValueError('Input global area integration is invalid')
            with output_gt3.open('rb') as f:
                for i in range(348):
                    f.seek(4+header_bytes+4+4, 1)
                    field = np.frombuffer(f.read(data_bytes), dtype='>f4').reshape((64, 128))
                    f.seek(4, 1)
                    source = ds['flux_ch4_prior'].isel(time=i).values.astype('float64')
                    if not np.isfinite(source).all():
                        raise ValueError('Non-finite input flux at month '+str(i))
                    old = np.sum(source * area[:, None], dtype='float64')
                    new = np.sum(field.astype('float64') * t42_area[:, None], dtype='float64')
                    scale = np.sum(np.abs(source) * area[:, None], dtype='float64')
                    if abs(old-new) > 2e-6*max(scale, 1e-30):
                        raise ValueError(f'Conservation failed {1998+i//12}-{1+i%12:02d}: '
                                         f'input={old:.12g} T42={new:.12g}')
        print('Conservation PASS: all 348 monthly global integrals, 1-degree vs T42')

    def check_prior_uncertainty_nc(self, filename):
        with xr.open_dataset(filename, decode_times=False) as ds:
            if ds.sizes.get('month') != 348 or ds.sizes.get('year') != 29:
                raise ValueError('Prior uncertainty incorrect month/year dimensions: ' + str(filename))
            if ds.sizes.get('land_region') != 54 or ds.sizes.get('datatype_count') != 2:
                raise ValueError('Prior uncertainty incorrect region/datatype dimensions')
            if 'sources_land' not in ds or 'time' not in ds:
                raise ValueError('Prior uncertainty missing variables')
            times = np.asarray(ds['time'].values, dtype=float)
            expected = np.array([year+(month-0.5)/12 for year in range(1998, 2027)
                                 for month in range(1, 13)])
            if not np.allclose(times, expected, rtol=0, atol=2e-4):
                raise ValueError('Prior uncertainty incorrect decimal years')
            if not bool(np.isfinite(ds['sources_land']).all()):
                raise ValueError('Prior uncertainty non-finite land uncertainties')
            uncertainty = ds['sources_land'].isel(datatype_count=1)
            if bool((uncertainty < 0).any()):
                raise ValueError('Prior uncertainty negative prior uncertainty')
            print('Prior uncertainty check PASS:', filename.name)

    # --------------------------------------------------------
    # --- Regional climatology format (54 regions x 12 months)
    def check_regional_table(self, filename):
        lines = filename.read_text().splitlines()
        if len(lines) != 55:
            raise ValueError('Regional table expected header plus 54 region rows')
        header = lines[0].split()
        if header[0] != 'reg#' or [int(v) for v in header[1:]] != list(range(1, 13)):
            raise ValueError('Regional table incorrect month headings')
        for reg, line in enumerate(lines[1:], 1):
            parts = line.split()
            if len(parts) != 13 or int(parts[0]) != reg:
                raise ValueError(f'Regional table incorrect region row {reg}')
            values = np.asarray([float(v) for v in parts[1:]])
            if not np.isfinite(values).all():
                raise ValueError(f'Regional table non-finite month, region {reg}')
        print('Regional table check PASS: 54 regions, 12 months')

    # --------------------------------------------------------
    # --- Write D1 namelist with matching descriptive Fortran group
    def write_case_namelist(self, operation, fields, input_files, output_files):
        names = {
            'regrid': 'd1_regrid_ch4_t42',
            'regions': 'd1_aggregate_ch4_regions',
            'uncertainty': 'd1_write_prior_uncertainty',
        }
        name = names[operation]
        cfg = self.config_dir / (name + '.nml')
        contents = '&' + name + '\n'
        contents += " hcase = '" + self.case + "'\n"
        for key, value in fields.items():
            if "'" in str(value) or '\n' in str(value):
                raise ValueError('Invalid Fortran namelist field: ' + key)
            if key in ('cfg_s_year', 'cfg_f_year'):
                contents += ' ' + key + ' = ' + str(int(value)) + '\n'
            else:
                contents += " " + key + " = '" + str(value) + "'\n"
        contents += '/\n'

        print('\nD1', operation, 'config:', cfg)
        for key, value in fields.items():
            print(' ', key, ':', value)

        for path in input_files:
            if not path.is_file():
                raise FileNotFoundError(operation + ' input missing: ' + str(path))
        for path in output_files:
            path.parent.mkdir(parents=True, exist_ok=True)
        self.config_dir.mkdir(parents=True, exist_ok=True)
        cfg.write_text(contents, encoding='utf-8')
        return cfg

    # --------------------------------------------------------
    # --- Compile when binary is missing or source is newer; then execute Fortran
    def run_exe(self, stage, config):
        import shutil
        exe_value = str(getattr(self, 'exe_d1_' + stage, '') or '').strip()
        if not exe_value:
            raise RuntimeError('Missing executable path for ' + stage + ' in _set_case.py')
        exe = Path(exe_value).expanduser().resolve()
        src = Path(self.fort_dir).expanduser().resolve() / (exe.name + '.f90')
        if not exe.is_file() or (src.is_file() and src.stat().st_mtime > exe.stat().st_mtime):
            if not src.is_file():
                raise FileNotFoundError('Missing Fortran source: ' + str(src))
            compiler = shutil.which('ifort') or shutil.which('ifx') or shutil.which('gfortran')
            nf_config = shutil.which('nf-config')
            if not compiler or not nf_config:
                raise RuntimeError('Fortran build requires ifort/ifx/gfortran and nf-config')
            flags = subprocess.check_output([nf_config, '--fflags'], text=True).split()
            libs = subprocess.check_output([nf_config, '--flibs'], text=True).split()
            exe.parent.mkdir(parents=True, exist_ok=True)
            cmd = [compiler, '-O2', str(src), '-o', str(exe)] + flags + libs
            print('Build:', ' '.join(cmd), flush=True)
            subprocess.run(cmd, check=True, cwd=str(exe.parent))
        print(stage, 'executable:', exe, flush=True)
        subprocess.run([str(exe), str(config)], check=True)

    # --------------------------------------------------------
    # --- D1 regrid: 1-degree prior NetCDF -> ACTM GT3
    def convert_flux_t42(self):
        input_file = Path(self.d1_regrid_input_dir).expanduser().resolve() / self.d1_regrid_input_nc
        output_dir = Path(self.d1_regrid_output_dir).expanduser().resolve()
        output_file = Path(self.d1_gt3_file).expanduser().resolve()
        if not input_file.is_file():
            raise FileNotFoundError(input_file)
        if output_file.is_file():
            self.check_gt3(output_file)
            self.check_flux_conservation(input_file, output_file)
            print('D1 regrid: reuse existing GT3:', output_file)
            return

        fields = {'inp_dir': str(input_file.parent) + os.sep,
                  'out_dir': str(output_dir) + os.sep,
                  'input_nc': self.d1_regrid_input_nc,
                  'input_case': self.d1_regrid_input_case,
                  'cfg_s_year': self.yr_s, 'cfg_f_year': self.yr_e}
        config = self.write_case_namelist('regrid', fields, [input_file], [output_file])
        self.run_exe('regrid', config)
        self.check_gt3(output_file)
        self.check_flux_conservation(input_file, output_file)

    # --------------------------------------------------------
    # --- Archive a stale output without overwriting previous backups
    @staticmethod
    def backup_stale(path):
        path = Path(path)
        if not path.exists():
            return
        from datetime import datetime
        suffix = datetime.now().strftime('%Y%m%d_%H%M%S')
        backup = Path(str(path) + '.backup_' + suffix)
        index = 1
        while backup.exists():
            backup = Path(str(path) + '.backup_' + suffix + '_' + str(index))
            index += 1
        path.rename(backup)
        print('Backup:', path, '->', backup)

    # --------------------------------------------------------
    # --- D1 regions: regional climatology from the current GT3
    def get_regional_flux(self):
        import hashlib
        inp = Path(self.d1_gt3_file).expanduser().resolve()
        mask = Path(self.d1_region_mask_file).expanduser().resolve()
        out = Path(self.d1_region_table_file).expanduser().resolve()
        marker = Path(str(out) + '.d1_1998-2026.sha256')

        if not inp.is_file() or not mask.is_file():
            raise FileNotFoundError('Missing D1 GT3 or regional mask: ' + str(inp) + ' / ' + str(mask))
        with inp.open('rb') as f:
            digest = hashlib.file_digest(f, 'sha256').hexdigest()

        valid = out.is_file() and marker.is_file() and marker.read_text().strip() == digest
        if valid:
            self.check_regional_table(out)
            print('D1 regions: reuse verified regional climatology:', out)
            return

        if out.exists() or marker.exists():
            print('D1 regions: old or unverified climatology; preserve and recompute for 1998-2026')
            self.backup_stale(out)
            self.backup_stale(marker)

        fields = {'flux_file': str(inp), 'mask_file': str(mask),
                  'uncr_file': str(out), 'cfg_s_year': self.yr_s, 'cfg_f_year': self.yr_e}
        config = self.write_case_namelist('regions', fields, [inp, mask], [out])
        self.run_exe('regions', config)
        self.check_nonempty(out)
        self.check_regional_table(out)
        marker.write_text(digest + '\n')

    # --------------------------------------------------------
    # --- D1 uncertainty: prior uncertainty files for p30/p50/p99
    def write_prior_uncertainty(self):
        inp = Path(self.d1_uncertainty_input_file).expanduser().resolve()
        prefix = Path(self.d1_uncertainty_output_prefix).expanduser().resolve()
        outputs = [Path(str(prefix) + str(pc) + '.nc') for pc in (30, 50, 99)]
        if not inp.is_file():
            raise FileNotFoundError('Missing D1 regional input: ' + str(inp))

        existing = [p.is_file() for p in outputs]
        if all(existing):
            try:
                for p in outputs:
                    self.check_prior_uncertainty_nc(p)
                    if p.stat().st_mtime_ns < inp.stat().st_mtime_ns:
                        raise ValueError('Prior uncertainty output predates current D1 regions: ' + str(p))
            except (ValueError, OSError) as exc:
                print('D1 uncertainty: stale or invalid outputs:', exc)
            else:
                print('D1 uncertainty: reuse verified prior uncertainty files')
                return

        if any(p.exists() for p in outputs):
            print('D1 uncertainty: preserve old uncertainty files and recompute')
            for p in outputs:
                self.backup_stale(p)

        fields = {'un_file': str(inp), 'nc_file': str(prefix),
                  'cfg_s_year': self.yr_s, 'cfg_f_year': self.yr_e}
        config = self.write_case_namelist('uncertainty', fields, [inp], [prefix])
        self.run_exe('uncertainty', config)
        for p in outputs:
            self.check_prior_uncertainty_nc(p)

    # --------------------------------------------------------
    @staticmethod
    def check_nonempty(path):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError('Missing or empty output: ' + str(path))

    # --------------------------------------------------------
    # --- Full D1 pipeline
    def run(self):
        self.get_inp_flx()
        self.convert_flux_t42()
        self.get_regional_flux()
        self.write_prior_uncertainty()
        print('\nD1 PRIOR: PASS')


# ============================================================
# === Stage interface called by _main.py
def run_stage(case='INCA', dry_run=False):
    case = str(case).upper()
    if case not in ('CYC', 'INCA'):
        raise ValueError('Unknown case: ' + case)

    print('D1 PRIOR:', case)
    print('Sequence: GCP NetCDF -> regrid T42 -> regional aggregation -> prior uncertainty')
    if dry_run:
        print('DRY-RUN: no files changed, no executables started')
        return

    # --- Support either version of SetCase case-selector during migration.
    # --- SetCase computes all paths inside __init__, so set before construction.
    old_gcp = os.environ.get('GCP_CASE')
    old_ch4 = os.environ.get('CH4_CASE')
    os.environ['GCP_CASE'] = case
    os.environ['CH4_CASE'] = case
    try:
        runner = PriorFluxes(case=case)
        runner.run()
    finally:
        if old_gcp is None:
            os.environ.pop('GCP_CASE', None)
        else:
            os.environ['GCP_CASE'] = old_gcp
        if old_ch4 is None:
            os.environ.pop('CH4_CASE', None)
        else:
            os.environ['CH4_CASE'] = old_ch4


# ============================================================
# === Standalone entry point
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=['CYC', 'INCA'], default='INCA')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    run_stage(case=args.case, dry_run=args.dry_run)
