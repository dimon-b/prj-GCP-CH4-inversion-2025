# -*- coding: utf-8 -*-
"""Check real TDI corrections, posterior GRD/GT3 budgets and seasonal changes."""

import argparse
import calendar
import struct
from pathlib import Path

import numpy as np
import pandas as pd


# =============================================================
# === Check posterior fluxes before ACTM forward simulation
class CheckFluxes:

    def __init__(self, case='INCA', icase='s057', prior='p30', unc='ctl',
                 mode='identity', root=None, predicted=None, output=None,
                 stage='sources', year_end=2024):

        self.case = case.upper()
        self.icase = icase
        self.prior = prior
        self.unc = unc
        self.mode = mode
        self.stage = stage
        self.year_end = int(year_end)
        self.root = Path(root).resolve() if root else Path(__file__).resolve().parent.parent
        self.res_dir = self.root / 'results2025' / self.case
        self.inv_dir = self.res_dir / f'{icase}_{prior}' / f'{icase}_{unc}_inv1'
        self.predicted = Path(predicted).resolve() if predicted else self.inv_dir / 'predicted.sources'
        self.prior_grd = self.root / 'fluxes/gcp2025_flux_prior/fch4_gcp2025_prior_inv1.grd'
        self.post_dir = self.res_dir / 'flux2d' / prior
        self.suffix = '_raw' if mode == 'identity' else ''
        self.post_grd = self.post_dir / f'{icase}_{unc}_inv1{self.suffix}.grd'
        self.post_gt3 = self.post_dir / f'{icase}_{unc}_inv1{self.suffix}.gt3'
        self.out_dir = Path(output).resolve() if output else self.res_dir / 'flux2d/diagnostics' / f'{icase}_{prior}_{unc}{self.suffix}'
        self.nlat = 180
        self.nlon = 360
        self.nregions = 54
        self.r_earth = 6371000.0
        self.regions = None
        self.monthly = None
        self.annual = None

        if self.case not in ('CYC', 'INCA'):
            raise ValueError('case must be CYC or INCA')
        if self.stage not in ('sources', 'full'):
            raise ValueError('stage must be sources or full')
        if self.mode not in ('identity', 'file'):
            raise ValueError('mode must be identity or file')

        self.read_sources()
        self.make_source_stats()
        if self.stage == 'full':
            self.check_grd_gt3()
        self.write_reports()

    # --- Read 54 TDI regional corrections, in Tg CH4 per 365-day year
    def read_sources(self):
        if not self.predicted.is_file():
            raise FileNotFoundError(f'Missing TDI source file: {self.predicted}')

        lines = self.predicted.read_text(encoding='utf-8').splitlines()
        records = []
        k = 0
        names = []
        while k < len(lines):
            name = lines[k].strip()
            k += 1
            if not name:
                continue
            if k >= len(lines):
                raise ValueError(f'Missing source header: {name}')
            header = lines[k].replace(',', ' ').split()
            k += 1
            if len(header) < 2:
                raise ValueError(f'Invalid source header: {name}')
            count = int(header[0])
            if count <= 0 or k + count > len(lines):
                raise ValueError(f'Invalid source count: {name}, count={count}')
            names.append(name)
            if name.startswith('land_'):
                for n in range(count):
                    row = lines[k+n].replace('D', 'E').replace('d', 'E').split()
                    if len(row) < 2:
                        raise ValueError(f'Incomplete correction row: {name}, {n}')
                    stamp = float(row[0])
                    year = int(stamp)
                    month = int(round((stamp-year)*12+0.5))
                    if not 1 <= month <= 12:
                        raise ValueError(f'Invalid year/month: {stamp}')
                    expected = year + (month-0.5)/12.0
                    if abs(stamp-expected) > 0.0002:
                        raise ValueError(f'Non-monthly decimal year: {stamp}')
                    correction = float(row[1])
                    if not np.isfinite(correction):
                        raise ValueError(f'Nonfinite source correction: {name} {year}-{month}')
                    records.append((name, year, month, correction))
            k += count

        df = pd.DataFrame(records, columns=['region', 'year', 'month', 'corr'])
        expected_names = [f'land_{n:02d}' for n in range(1, self.nregions+1)]
        if sorted(df.region.unique()) != sorted(expected_names):
            raise ValueError('TDI land regions are not exactly land_01...land_54')
        if df.duplicated(['region', 'year', 'month']).any():
            raise ValueError('Duplicate region-month in predicted.sources')
        years = sorted(df.year.unique())
        if years != list(range(min(years), max(years)+1)):
            raise ValueError('Gaps in TDI years')
        if len(df) != self.nregions * len(years) * 12:
            raise ValueError('Missing TDI region-month corrections')
        if self.year_end < min(years) or self.year_end > max(years):
            raise ValueError('year_end is outside source time range')
        df['days'] = [calendar.monthrange(int(y), int(m))[1] for y,m in zip(df.year, df.month)]
        df['corr_tg_month'] = df['corr'] * df['days'] / 365.0
        self.regions = df
        self.year_start = min(years)
        self.year_last = max(years)
        self.num_months = len(years) * 12
        print('\n=== TDI POST flux check ====================================')
        print(' Case          :', self.case, self.icase, self.prior, self.unc)
        print(' Mode          :', self.mode)
        print(' Source file   :', self.predicted)
        print(' Regions       :', self.nregions)
        print(' Months        :', self.num_months)
        print(' Analysis      :', self.year_start, '-', self.year_end)

    # --- Annual and regional budgets of the inverse correction
    def make_source_stats(self):
        df = self.regions
        months = df.groupby(['year','month'], as_index=False).agg(
            correction_tg_yr=('corr','sum'),
            correction_tg_month=('corr_tg_month','sum'),
            abs_region_corr_tg_yr=('corr', lambda v: np.abs(v).sum()),
        )
        self.monthly = months
        self.annual = months.groupby('year', as_index=False).agg(
            correction_tg=('correction_tg_month', 'sum'),
            mean_abs_region_corr_tg_yr=('abs_region_corr_tg_yr', 'mean'),
        )
        selected = self.annual.loc[self.annual.year.le(self.year_end)]
        print(' Mean correction [Tg/yr] :', f'{selected.correction_tg.mean():.3f}')
        print(' Last annual correction  :', f'{selected.iloc[-1].correction_tg:.3f}', 'Tg/yr')
        reg = df.loc[df.year.le(self.year_end)].groupby('region')['corr_tg_month'].sum()
        reg = (reg / (self.year_end-self.year_start+1)).sort_values(key=abs, ascending=False)
        print(' Strongest mean regional corrections [Tg/yr]:')
        for name, value in reg.head(10).items():
            print(f'  {name:8s}  {value:10.3f}')

    # --- Compare actual 1-degree GRD with TDI predicted source budget
    def check_grd_gt3(self):
        for file in (self.prior_grd, self.post_grd, self.post_gt3):
            if not file.is_file():
                raise FileNotFoundError(f'Missing real flux input/output: {file}')

        expected = self.num_months * self.nlat * self.nlon * 4
        for file in (self.prior_grd, self.post_grd):
            if file.stat().st_size != expected:
                raise ValueError(f'Unexpected GRD size: {file}: {file.stat().st_size}, expected {expected}')

        prior = np.memmap(self.prior_grd, dtype=np.float32, mode='r', shape=(self.num_months,self.nlat,self.nlon))
        post = np.memmap(self.post_grd, dtype=np.float32, mode='r', shape=(self.num_months,self.nlat,self.nlon))
        angles = np.linspace(-np.pi/2, np.pi/2, self.nlat+1)
        areas = self.r_earth**2 * (2*np.pi/self.nlon) * np.diff(np.sin(angles))
        years_months = self.monthly.sort_values(['year','month'])
        if len(years_months) != self.num_months:
            raise ValueError('Mismatch between GRD month count and predicted.sources')

        expected_gt3 = self.num_months * (4+1024+4+4+128*64*4+4)
        if self.post_gt3.stat().st_size != expected_gt3:
            raise ValueError(f'Unexpected GT3 size: {self.post_gt3.stat().st_size}, expected {expected_gt3}')

        lc_by_year = {}
        if self.mode == 'file':
            lc_file = self.res_dir / 'losscorr' / f'gcp_LC_{self.icase}_{self.prior}_{self.unc}_inv1.txt'
            if not lc_file.is_file():
                raise FileNotFoundError(f'Missing annual LC: {lc_file}')
            df_lc = pd.read_csv(lc_file,sep=r'\s+')
            if not {'year','LC_fact'}.issubset(df_lc.columns):
                raise ValueError('LC file must contain year and LC_fact')
            if df_lc.year.duplicated().any() or not np.isfinite(df_lc.LC_fact).all() or (df_lc.LC_fact<=0).any():
                raise ValueError('Invalid or repeated annual LC factors')
            lc_by_year = dict(zip(df_lc.year.astype(int),df_lc.LC_fact.astype(float)))
            missing = set(years_months.year)-set(lc_by_year)
            if missing:
                raise ValueError(f'Missing LC years: {sorted(missing)}')

        # --- Gaussian cell areas of ACTM T42 for checking GT3 flux budget
        from numpy.polynomial.legendre import leggauss
        _, weights = leggauss(64)
        gauss_area = 2*np.pi*self.r_earth**2 * weights/128

        records = []
        with self.post_gt3.open('rb') as f:
            for idx, row in enumerate(years_months.itertuples(index=False)):
                yr = int(row.year)
                mo = int(row.month)
                grid_prior = np.asarray(prior[idx], dtype=np.float64)
                grid_post = np.asarray(post[idx], dtype=np.float64)
                if not (np.isfinite(grid_prior).all() and np.isfinite(grid_post).all()):
                    raise ValueError(f'Nonfinite GRD value at {yr}-{mo:02d}')
                prior_global = (grid_prior.sum(axis=1) * areas).sum() * 365*86400*1e-9
                post_global = (grid_post.sum(axis=1) * areas).sum() * 365*86400*1e-9
                source_expected = float(row.correction_tg_yr)
                grd_correction = post_global-prior_global
                lc_factor = lc_by_year[yr] if self.mode == 'file' else 1.0
                actual_correction = post_global/lc_factor-prior_global
                error = actual_correction-source_expected
                tolerance = 0.02 + abs(source_expected)*0.0001
                if abs(error) > tolerance:
                    raise ValueError(f'TDI-GRD budget mismatch {yr}-{mo:02d}: expected {source_expected:.5f}, actual {actual_correction:.5f} Tg/yr')

                # --- GTOOL3 sequential unformatted records, big endian
                hlen_raw = f.read(4)
                if len(hlen_raw) != 4 or struct.unpack('>i', hlen_raw)[0] != 1024:
                    raise ValueError(f'Invalid GT3 header record {yr}-{mo:02d}')
                header = f.read(1024)
                if len(header) != 1024 or struct.unpack('>i', f.read(4))[0] != 1024:
                    raise ValueError(f'Invalid GT3 header end {yr}-{mo:02d}')
                if b'kg-CH4/m2/s' not in header:
                    raise ValueError(f'Unexpected GT3 units {yr}-{mo:02d}')
                data_len = struct.unpack('>i', f.read(4))[0]
                if data_len != 128*64*4:
                    raise ValueError(f'Invalid GT3 data record length {yr}-{mo:02d}')
                data = f.read(data_len)
                if len(data) != data_len or struct.unpack('>i', f.read(4))[0] != data_len:
                    raise ValueError(f'Invalid GT3 data record {yr}-{mo:02d}')
                t42 = np.frombuffer(data, dtype='>f4').reshape(64,128)
                if not np.isfinite(t42).all():
                    raise ValueError(f'Nonfinite GT3 values {yr}-{mo:02d}')
                t42_global = (t42.astype(float).sum(axis=1)*gauss_area).sum()*365*86400*1e-9
                mismatch = t42_global-post_global
                tolerance = 0.02+abs(post_global)*0.0001
                if abs(mismatch) > tolerance:
                    raise ValueError(f'GRD-GT3 budget mismatch {yr}-{mo:02d}: {mismatch:.5f} Tg/yr')

                records.append((yr,mo,prior_global,post_global,grd_correction,source_expected,
                                t42_global, mismatch, int((grid_prior<0).sum()),int((grid_post<0).sum())))
            if f.read(1):
                raise ValueError('GT3 file has trailing data')

        del prior, post
        check = pd.DataFrame(records, columns=[
            'year','month','prior_tg_yr','post_tg_yr','grd_corr_tg_yr','tdi_corr_tg_yr',
            't42_tg_yr','gt3_minus_grd_tg_yr','n_negative_prior','n_negative_post'])
        self.monthly = self.monthly.merge(check, on=['year','month'],validate='one_to_one')
        self.annual = self.annual.merge(
            check.assign(days=[calendar.monthrange(int(y),int(m))[1] for y,m in zip(check.year,check.month)])
                 .assign(prior_tg=lambda d:d.prior_tg_yr*d.days/365,
                         post_tg=lambda d:d.post_tg_yr*d.days/365)
                 .groupby('year',as_index=False)[['prior_tg','post_tg']].sum(),
            on='year',validate='one_to_one')
        self.annual['grd_correction_tg'] = self.annual.post_tg-self.annual.prior_tg
        sel = self.annual.loc[self.annual.year.le(self.year_end)]
        print(' Mean prior [Tg/yr]      :', f'{sel.prior_tg.mean():.3f}')
        print(' Mean posterior [Tg/yr]  :', f'{sel.post_tg.mean():.3f}')
        print(' Mean corrected [Tg/yr]  :', f'{sel.grd_correction_tg.mean():.3f}')
        print(' GRD negative grid-month :', int(self.monthly.n_negative_post.sum()))
        print(' TDI -> GRD -> GT3 budget: PASS')

    # --- Reports for INCA/CYC comparison and checking seasonality
    def write_reports(self):
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.monthly.to_csv(self.out_dir / 'monthly_budget.csv', index=False, float_format='%.10g')
        self.annual.to_csv(self.out_dir / 'annual_budget.csv', index=False, float_format='%.10g')
        selected_regions = self.regions.loc[self.regions.year.le(self.year_end)]
        selected_regions.to_csv(self.out_dir / 'regional_corrections.csv',index=False,float_format='%.10g')
        regional_clim = selected_regions.groupby(['region','month'], as_index=False)['corr'].mean()
        regional_clim.to_csv(self.out_dir / 'regional_correction_climatology.csv',index=False,float_format='%.10g')
        global_clim = self.monthly.loc[self.monthly.year.le(self.year_end)].groupby('month',as_index=False)['correction_tg_yr'].mean()
        global_clim.to_csv(self.out_dir / 'global_correction_climatology.csv',index=False,float_format='%.10g')
        seasonal_span = global_clim.correction_tg_yr.max()-global_clim.correction_tg_yr.min()
        print(' Global correction seasonal range [Tg/yr]:',f'{seasonal_span:.3f}')

        # --- Seasonal shape: average monthly flux after removing annual mean
        if self.stage == 'full':
            monthly = self.monthly.loc[self.monthly.year.le(self.year_end)].copy()
            for var in ('prior_tg_yr','post_tg_yr'):
                monthly[var+'_anomaly'] = monthly[var]-monthly.groupby('year')[var].transform('mean')
            season = monthly.groupby('month',as_index=False)[[
                'prior_tg_yr_anomaly','post_tg_yr_anomaly']].mean()
            season['change_tg_yr'] = season.post_tg_yr_anomaly-season.prior_tg_yr_anomaly
            season.to_csv(self.out_dir / 'seasonal_flux_changes.csv',index=False,float_format='%.10g')
            amp_prior = season.prior_tg_yr_anomaly.max()-season.prior_tg_yr_anomaly.min()
            amp_post = season.post_tg_yr_anomaly.max()-season.post_tg_yr_anomaly.min()
            print(' Flux seasonal range prior/post [Tg/yr]:', f'{amp_prior:.3f}',f'{amp_post:.3f}')

        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            annual = self.annual.loc[self.annual.year.le(self.year_end)]
            fig,ax=plt.subplots(figsize=(11,4.5))
            ax.plot(annual.year,annual.correction_tg,'o-',label='TDI correction')
            ax.axhline(0,color='0.5',lw=.8)
            ax.set_xlabel('Year')
            ax.set_ylabel('Global correction (Tg CH4/yr)')
            ax.grid(linestyle=':')
            ax.legend()
            fig.tight_layout()
            fig.savefig(self.out_dir / 'global_correction.png',dpi=180)
            plt.close(fig)
            fig,ax=plt.subplots(figsize=(9,4.5))
            ax.plot(global_clim.month,global_clim.correction_tg_yr,'o-',label='Global TDI correction')
            ax.axhline(0,color='0.5',lw=.8)
            ax.set_xticks(range(1,13))
            ax.set_xlabel('Month')
            ax.set_ylabel('Global monthly correction (Tg CH4/yr equivalent)')
            ax.grid(linestyle=':')
            ax.legend()
            fig.tight_layout()
            fig.savefig(self.out_dir / 'global_correction_seasonality.png',dpi=180)
            plt.close(fig)
            if self.stage=='full':
                fig,ax=plt.subplots(figsize=(9,4.5))
                ax.plot(season.month,season.prior_tg_yr_anomaly,'o-',label='Prior')
                ax.plot(season.month,season.post_tg_yr_anomaly,'s-',label='Posterior')
                ax.set_xticks(range(1,13))
                ax.set_xlabel('Month')
                ax.set_ylabel('Global emission seasonal anomaly (Tg CH4/yr)')
                ax.grid(linestyle=':')
                ax.legend()
                fig.tight_layout()
                fig.savefig(self.out_dir / 'seasonal_flux_change.png',dpi=180)
                plt.close(fig)
        except ImportError:
            print(' Matplotlib missing: CSV reports created without PNG')

        print(' Reports        :', self.out_dir)
        print(' POST flux check: PASS')


# =============================================================
# === Command line
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check posterior TDI/GRD/GT3 and seasonal flux changes')
    parser.add_argument('--case',choices=['CYC','INCA'],required=True)
    parser.add_argument('--icase',default='s057')
    parser.add_argument('--prior',default='p30')
    parser.add_argument('--unc',default='ctl')
    parser.add_argument('--mode',choices=['identity','file'],default='identity')
    parser.add_argument('--stage',choices=['sources','full'],default='sources')
    parser.add_argument('--root',default=None)
    parser.add_argument('--predicted',default=None)
    parser.add_argument('--output',default=None)
    parser.add_argument('--year-end',type=int,default=2024)
    args = parser.parse_args()
    CheckFluxes(case=args.case,icase=args.icase,prior=args.prior,unc=args.unc,
                mode=args.mode,root=args.root,predicted=args.predicted,
                output=args.output,stage=args.stage,year_end=args.year_end)
