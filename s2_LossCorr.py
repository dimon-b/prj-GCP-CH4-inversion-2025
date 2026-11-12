# -*- coding: utf-8 -*-
"""
Created:    21/11/08 17:48
Project:    GCP/GMB 2025 project [server part]
@author:    Dmitry Belikov
"""

import os
import xarray as xr
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

import _set_case


class LossCorr(_set_case.SetCase):

    def __init__(self):
        super().__init__()

        # --- input
        self.f_blos = os.path.join(self.inv_lsc_dir, 'frt_burd_')

        # --- output
        self.f_bldd = os.path.join(self.inv_lsc_dir, 'gcp_burd_add_')
        self.f_outl = os.path.join(self.inv_lsc_dir, 'gcp_LC_')

        print('\n\tInput file   : ', self.f_blos)
        print('\tOutput file  : ', self.f_bldd)
        print('\t             : ', self.f_outl)

        # --- obs
        self.site = ['SPO']
        self.site = ['SYO']

        v1 = False
        if v1:
            # --- run correction v1
            self.clc_lcorr_v1()
        else:
            # --- run correction post
            if self.hcase == 'CYC':
                self.post_model_conc = '../results2025/CYC/nc_out/MIROC4-ACTM_CH4conc_GMB_SURF_OH_Transcom.nc'
            else:
                self.post_model_conc = '../results2025/INCA/nc_out/MIROC4-ACTM_CH4conc_GMB_SURF_OH_INCA.nc'
            self.clc_lcorr_post()

    # --- get loss correction from the raw POST run
    def clc_lcorr_post(self):
        years = pd.Index(range(self.yr_s, self.yr_e + 1), name='year')
        edge_years = years.intersection([1999, 2025])
        post_years = years.difference(edge_years)

        # - station inventory
        df_st = pd.read_csv(self.val_sts)

        def get_station(site_ids):
            selected = df_st['site_id'].str.strip().isin(site_ids)
            fixed = df_st['project'].str.contains('surface|tower', case=False, na=False)
            stations = df_st.loc[selected & fixed].dropna(subset=['f_name']).drop_duplicates('f_name')

            if len(stations) != 1:
                raise ValueError(f'Expected one observation file for {site_ids} in {self.val_sts}')

            return stations.iloc[0]

        row = get_station(self.site)
        lat, lon, alt = row[['lat', 'lon', 'alt']].to_numpy(dtype=float)
        if not (-90 <= lat <= 90 and -180 <= lon <= 180 and -999 < alt <= 11000):
            raise ValueError(f'Invalid station coordinates for {self.site}')

        p_obs = 101325.0 * (1.0 - 2.25577e-5 * alt) ** 5.25588

        # - read annual burden, loss, and December reference CH4
        def get_burdloss():
            fname = self.f_blos + self.tracer + '.txt'
            cols = ['year', 'ch4_ref', 'bch4', 'loh', 'lcl', 'lod']
            df_bl = pd.read_csv(fname, sep=r'\s+', usecols=cols)
            df_bl = df_bl.apply(pd.to_numeric, errors='raise')
            df_bl = df_bl.loc[df_bl['year'].isin(years)].drop_duplicates()
            df_bl = df_bl.set_index('year', verify_integrity=True).reindex(years)
            df_bl['lch4'] = df_bl['loh'] + df_bl['lcl'] + df_bl['lod']

            if not np.isfinite(df_bl[['bch4', 'lch4']].to_numpy()).all():
                raise ValueError(f'{fname}: missing or invalid burden/loss data')
            if (df_bl[['bch4', 'lch4']] <= 0).any().any():
                raise ValueError(f'{fname}: invalid burden or loss')

            ref = df_bl.loc[edge_years, 'ch4_ref']
            valid_ref = np.isfinite(ref) & ref.gt(0)
            if not valid_ref.all():
                raise ValueError(f'{fname}: invalid December ch4_ref for {ref.index[~valid_ref].tolist()}')

            return df_bl[['ch4_ref', 'bch4', 'lch4']]

        # - read monthly observations and screen outliers
        def read_obs(station):
            fname = os.path.join(self.val_obs_dir, station['f_name'])
            print('CH4 obs from', fname)

            time_cols = ['year', 'month', 'day', 'hour', 'minute']
            df_ob = pd.read_csv(fname, sep='\t', usecols=time_cols + ['value'])
            df_ob = df_ob.apply(pd.to_numeric, errors='coerce')
            df_ob = df_ob.replace([np.inf, -np.inf], np.nan)

            valid = (
                df_ob['year'].isin(years)
                & df_ob['hour'].isin(range(24))
                & df_ob['minute'].isin(range(60))
                & df_ob['value'].gt(0.0)
            )
            df_ob = df_ob.loc[valid].copy()
            df_ob['datetime'] = pd.to_datetime(df_ob[time_cols], errors='coerce')

            by_time = df_ob.dropna(subset=['datetime']).groupby('datetime')['value'].mean()
            obs = by_time.groupby(by_time.index.to_period('M')).mean().rename('ch4_obs')
            if obs.empty:
                raise ValueError(f'{fname}: no valid observation months')

            # - fit a cubic trend and annual cycle using Huber loss
            month = obs.index.asi8
            t = (month - month.mean()) / 120.0
            phase = 2.0 * np.pi * (obs.index.month - 1) / 12.0
            design = np.column_stack([
                np.ones(len(obs)), t, t**2, t**3, np.sin(phase), np.cos(phase),
            ])
            values = obs.to_numpy(dtype=float)
            initial, _, rank, _ = np.linalg.lstsq(design, values, rcond=None)

            if rank < design.shape[1] or len(obs) <= rank:
                raise ValueError(f'{fname}: insufficient seasonal coverage for the observation fit')

            fit = least_squares(
                lambda coef: design @ coef - values,
                initial,
                loss='huber',
                f_scale=self.val_obs_diff_threshold,
            )
            if not fit.success:
                raise RuntimeError(f'Observation fit failed: {fit.message}')

            residual = obs - design @ fit.x
            bad_obs = residual.abs().gt(self.val_obs_diff_threshold)
            print(f'\tDropped observation months: {bad_obs.sum()}')
            if bad_obs.any():
                print(residual.loc[bad_obs].round(1).to_string())

            return obs.loc[~bad_obs]

        # - extract monthly POST CH4 at the selected station
        def read_model(obs):
            fname = self.post_model_conc
            print('Raw POST CH4 from', fname)

            with xr.open_dataset(fname, engine='h5netcdf') as ds:
                for dim, coord in {'x': 'longitude', 'y': 'latitude', 'plev': 'level'}.items():
                    if dim in ds.dims:
                        if coord in ds.variables:
                            ds = ds.swap_dims({dim: coord})
                        else:
                            ds = ds.rename({dim: coord})

                lat_m = ds['latitude'].to_numpy().astype(float)
                lon_m = ds['longitude'].to_numpy().astype(float)
                levels = ds['level'].to_numpy().astype(float)

                if not np.isfinite(lat_m).all() or not np.isfinite(lon_m).all():
                    raise ValueError(f'{fname}: latitude and longitude must be finite')
                if not (np.isfinite(levels) & (levels > 0) & (levels <= 1.01)).all():
                    raise ValueError(f'{fname}: expected positive sigma levels up to 1.01')

                iy = int(np.abs(lat_m - lat).argmin())
                ix = int(np.abs((lon_m - lon + 180) % 360 - 180).argmin())
                ds_p = ds[['ch4_conc', 'Ps']].isel(latitude=iy, longitude=ix).load()

                if ds_p['Ps'].dims != ('time',):
                    raise ValueError(f"{fname}: expected Ps(time) at station, found {ds_p['Ps'].dims}")

                ch4 = ds_p['ch4_conc'].transpose('time', 'level').to_numpy()  # ppb
                ps = ds_p['Ps'].to_numpy()  # Pa

                # - use the averaging interval to identify the month
                time = ds_p['time']
                bounds_name = time.attrs.get('bounds') or time.encoding.get('bounds') or 'time_bnds'
                if bounds_name in ds:
                    bounds = ds[bounds_name].transpose('time', ...).to_numpy()
                    time = xr.DataArray(
                        [start + (end - start) / 2 for start, end in bounds],
                        dims='time',
                    )

                dates = [
                    f'{year:04d}-{month:02d}'
                    for year, month in zip(time.dt.year.values, time.dt.month.values)
                ]

            df_m = pd.DataFrame(ch4, index=pd.PeriodIndex(dates, freq='M'))
            df_m['ps'] = ps
            df_m = df_m.where(np.isfinite(df_m) & df_m.gt(0)).groupby(level=0).mean()
            ps_m = df_m.pop('ps')

            # - interpolate CH4 in log pressure
            post = pd.Series(index=obs.index, dtype=float, name='ch4_ref')
            for month in obs.index.intersection(df_m.index):
                p_mod = levels * ps_m.loc[month]
                values = df_m.loc[month].to_numpy()
                valid = np.isfinite(p_mod) & (p_mod > 0) & np.isfinite(values)

                if valid.any():
                    p_mod = p_mod[valid]
                    values = values[valid]
                    order = np.argsort(p_mod)
                    post.loc[month] = np.interp(
                        np.log(p_obs), np.log(p_mod[order]), values[order],
                    )

            missing = post_years.difference(post.dropna().index.year)
            if len(missing):
                coverage = pd.DataFrame({
                    'OBS': obs.groupby(obs.index.year).size(),
                    'NC_time': df_m.groupby(df_m.index.year).size(),
                    'CH4': df_m.notna().any(axis=1).groupby(df_m.index.year).sum(),
                    'Ps': ps_m.notna().groupby(ps_m.index.year).sum(),
                    'paired': post.notna().groupby(post.index.year).sum(),
                })
                print(coverage.reindex(missing).fillna(0).astype(int).to_string())
                raise ValueError(f'{self.site}: no paired POST and observation months for {missing.tolist()}')

            return post

        # - annual means from matched POST and observation months
        obs = read_obs(row)
        if len(post_years):
            obs_post = obs.loc[obs.index.year.isin(post_years)]
            paired = pd.concat([obs_post, read_model(obs_post)], axis=1).dropna()
            by_year = paired.groupby(paired.index.year)
            df_yr = by_year.mean().reindex(years)
            n_months = by_year.size().reindex(years, fill_value=0)
        else:
            df_yr = pd.DataFrame(index=years, columns=['ch4_obs', 'ch4_ref'], dtype=float)
            n_months = pd.Series(0, index=years)

        # - frt_burd ch4_ref represents December at SPO
        if len(edge_years):
            obs_spo = obs if row['site_id'].strip() == 'SPO' else read_obs(get_station(['SPO']))
            dec_months = pd.PeriodIndex([f'{year}-12' for year in edge_years], freq='M')
            obs_dec = obs_spo.reindex(dec_months).set_axis(edge_years)
            missing_dec = obs_dec.index[obs_dec.isna()]
            if len(missing_dec):
                raise ValueError(f'SPO: missing December observations for {missing_dec.tolist()}')

            df_yr.loc[edge_years, 'ch4_obs'] = obs_dec
            n_months.loc[edge_years] = 1
            print(f'\tUsing December SPO and frt_burd ch4_ref for {edge_years.tolist()}')

        for trac in self.invcases:
            self.tracer = trac
            df_bl = get_burdloss()
            df_rs = df_yr.join(df_bl[['bch4', 'lch4']])
            df_rs.loc[edge_years, 'ch4_ref'] = df_bl.loc[edge_years, 'ch4_ref']
            df_rs = df_rs.reset_index()

            df_rs['d_ch4'] = df_rs['ch4_ref'] - df_rs['ch4_obs']
            df_rs['LC_fact'] = df_rs['ch4_obs'] / df_rs['ch4_ref']
            df_rs['LC'] = df_rs['lch4'] * (1.0 - df_rs['LC_fact'])
            df_rs['bch4_corr'] = df_rs['bch4'] * df_rs['LC_fact']

            df_rs = df_rs[[
                'year', 'ch4_obs', 'ch4_ref', 'd_ch4', 'bch4',
                'bch4_corr', 'lch4', 'LC', 'LC_fact',
            ]]
            print(df_rs.assign(n_months=n_months.to_numpy()))

            df_rs[['year', 'LC_fact']].to_csv(self.f_outl + trac + '.txt', sep='\t', index=False)
            df_rs.round(8).to_csv(self.f_bldd + trac + '.txt', sep='\t', index=False)

    # --- get loss correction coefficient
    def clc_lcorr_v1(self):
        years = pd.Index(range(self.yr_s, self.yr_e + 1), name='year')

        # - read burden/loss data
        def get_burdloss():
            fname = self.f_blos + self.tracer + '.txt'
            cols = ['year', 'ch4_ref', 'bch4', 'loh', 'lcl', 'lod']
            df_bl = pd.read_csv(fname, sep=r'\s+', usecols=cols)
            df_bl = df_bl.apply(pd.to_numeric, errors='raise')
            df_bl = df_bl.loc[df_bl['year'].isin(years)].drop_duplicates()
            df_bl = df_bl.set_index('year', verify_integrity=True).reindex(years)
            df_bl['lch4'] = df_bl['loh'] + df_bl['lcl'] + df_bl['lod']

            if not np.isfinite(df_bl.to_numpy()).all():
                raise ValueError(f'{fname}: missing or invalid burden/loss data')
            if (df_bl[['ch4_ref', 'bch4', 'lch4']] <= 0).any().any():
                raise ValueError(f'{fname}: invalid concentration, burden, or loss')

            return df_bl[['ch4_ref', 'bch4', 'lch4']]

        # - read and screen observations
        def read_obs():
            df_st = pd.read_csv(self.val_sts)
            selected = df_st['site_id'].str.strip().isin(self.site)
            fixed = df_st['project'].str.contains('surface|tower', case=False, na=False)
            files = df_st.loc[selected & fixed, 'f_name'].dropna().unique()
            if len(files) != 1:
                raise ValueError(f'Expected one observation file for {self.site} in {self.val_sts}')

            fname = os.path.join(self.val_obs_dir, files[0])
            print('\tCH4 obs from :', fname)

            time_cols = ['year', 'month', 'day', 'hour', 'minute']
            df_ob = pd.read_csv(fname, sep='\t', usecols=time_cols + ['value'])
            df_ob = df_ob.apply(pd.to_numeric, errors='coerce')
            df_ob = df_ob.replace([np.inf, -np.inf], np.nan)
            valid = (df_ob['year'].isin(years)
                     & df_ob['hour'].isin(range(24))
                     & df_ob['minute'].isin(range(60))
                     & df_ob['value'].gt(0.0))
            df_ob = df_ob.loc[valid].copy()
            df_ob['datetime'] = pd.to_datetime(df_ob[time_cols], errors='coerce')

            by_time = df_ob.dropna(subset=['datetime']).groupby('datetime')['value'].mean()
            obs = by_time.groupby(by_time.index.to_period('M')).mean().rename('ch4_obs')

            # - fit the trend and seasonal cycle, then drop outlier months
            if not obs.empty:
                month = obs.index.asi8
                t = (month - month.mean()) / 120.0
                seasonal = pd.get_dummies(obs.index.month, drop_first=True, dtype=float)
                design = np.column_stack([np.ones(len(obs)), t, t ** 2, t ** 3, seasonal])
                values = obs.to_numpy(dtype=float)
                initial, _, rank, _ = np.linalg.lstsq(design, values, rcond=None)

                if len(obs) > rank and rank == design.shape[1]:
                    fit = least_squares(lambda coef: design @ coef - values,
                                        initial,
                                        loss='huber',
                                        f_scale=self.val_obs_diff_threshold)
                    if not fit.success:
                        raise RuntimeError(f'Observation fit failed: {fit.message}')

                    residual = obs - design @ fit.x
                    bad_obs = residual.abs().gt(self.val_obs_diff_threshold)
                    print(f'\tDropped observation months: {bad_obs.sum()}')
                    print()

                    if bad_obs.any():
                        print(residual.loc[bad_obs].round(1).to_string())
                    obs = obs.loc[~bad_obs]

            # - match the December ch4_ref in frt_burd
            obs_dec = obs.loc[obs.index.month == 12].copy()
            obs_dec.index = obs_dec.index.year.rename('year')
            obs_dec = obs_dec.reindex(years)

            if self.site == ['SPO'] and 2025 in years and pd.isna(obs_dec.loc[2025]):
                obs_dec.loc[2025] = 1883.0
                print('SPO December 2025: using 1883 ppb')

            missing = obs_dec.index[obs_dec.isna()]
            if len(missing):
                raise ValueError(f'Missing December observations for years {missing.tolist()}')

            return obs_dec

        # --- obs
        obs_dec = read_obs()

        for trac in self.invcases:
            self.tracer = trac
            df_rs = get_burdloss().join(obs_dec).reset_index()

            df_rs['d_ch4'] = df_rs['ch4_ref'] - df_rs['ch4_obs']
            df_rs['LC_fact'] = df_rs['ch4_obs'] / df_rs['ch4_ref']
            df_rs['LC'] = df_rs['lch4'] * (1.0 - df_rs['LC_fact'])
            df_rs['bch4_corr'] = df_rs['bch4'] * df_rs['LC_fact']

            new_order = ['year', 'ch4_obs', 'ch4_ref', 'd_ch4', 'bch4',
                         'bch4_corr', 'lch4', 'LC', 'LC_fact', ]
            df_rs = df_rs[new_order]
            print(df_rs)

            # - to files
            df_rs[['year', 'LC_fact']].to_csv(self.f_outl + self.tracer + '.txt', sep='\t', index=False)
            df_rs.round(8).to_csv(self.f_bldd + self.tracer + '.txt', sep='\t', index=False)
