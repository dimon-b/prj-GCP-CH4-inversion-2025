# -*- coding: utf-8 -*-
"""
Finalized:  26/xx/xx
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import os
from contextlib import ExitStack

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

import _set_case


class CheckGR(_set_case.SetCase):

    def __init__(self):
        super().__init__()

        self.fig_name = 'val4_check_GR.png'
        self.check_GR()

    # - compare local Mauna Loa growth rates from the same station inventory
    def check_GR(self, **obs_settings):
        results = self.get_model_obs(self._get_sites(['MLO']), **obs_settings)
        if not results:
            print('\nNo Mauna Loa record meets the common-month requirement')
            return

        # - keep one station record and use the one with the most paired months
        row, df_m = max(results, key=lambda item: len(item[1]))
        print(f'\nMauna Loa record: {row["f_name"]} ({len(df_m)} common months)')

        cols = ['o_ch4', *self.val_model_cases]
        df = self._monthly_time_frame(df_m[cols])

        # - centered 12-month running mean, then differences of annual means
        df_smooth = df.rolling(window=12, center=True, min_periods=6).mean()
        df_yr = df_smooth.resample('YS').mean()
        df_gr = df_yr.diff()

        if not df_gr['o_ch4'].notna().any():
            print('\nNo consecutive annual means for Mauna Loa growth rates')
            return

        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(df_gr.index, df_gr['o_ch4'], color='k', lw=2.0, label='MLO Obs')

        for col, case in self.val_model_cases.items():
            residual = df_gr[col] - df_gr['o_ch4']
            bias = residual.mean()
            rmse = np.sqrt((residual ** 2).mean())
            corr = df_gr[col].corr(df_gr['o_ch4'], method='pearson')
            label = case['label']

            print(f'{label}: mean CH4 = {df_yr[col].mean():.2f} ppb, '
                  f'mean growth rate = {df_gr[col].mean():.2f} ppb/yr')
            print(f'{label}: bias = {bias:.2f} ppb/yr, RMSE = {rmse:.2f} ppb/yr, '
                  f'Pearson r = {corr:.3f}')

            ax.plot(df_gr.index, df_gr[col], color=case['color'], marker=case['marker'],
                    ls='--', lw=1.0, markersize=4, label=label)

        x_start = pd.Timestamp(year=self.years[0], month=10, day=1)
        x_end = pd.Timestamp(year=self.years[1], month=7, day=31)
        ax.set_xlim(x_start, x_end)
        ax.set_title('CH₄ Growth Rate at Mauna Loa')
        ax.set_ylabel('Growth rate (ppb/yr)')
        ax.set_xlabel('Year')
        ax.grid(linestyle=':')
        ax.legend(loc='best', fontsize=10)

        # - save
        os.makedirs(self.plt_dir, exist_ok=True)
        f_out = os.path.join(self.plt_dir, self.fig_name)
        fig.tight_layout()
        try:
            fig.savefig(f_out, format='png', dpi=600, bbox_inches='tight')
        finally:
            plt.close(fig)

    # - prepare paired monthly observations and model values for each station record
    def get_model_obs(self, df_st=None, ch4_factor=1.0, level_type='sigma',
                      ps_factor=1.0, min_months=15 * 12):

        if level_type not in ['sigma', 'Pa', 'hPa']:
            raise ValueError("level_type must be 'sigma', 'Pa' or 'hPa'")

        if min_months < 2:
            raise ValueError('min_months must be at least 2')

        target_cols = list(self.val_model_cases)
        if not target_cols:
            raise ValueError('Configure at least one model case in self.model_cases')

        if df_st is None:
            df_st = self._get_sites(self.sites)

        if df_st.empty:
            return []

        model_vars = ['ch4_conc', 'Ps'] if level_type == 'sigma' else ['ch4_conc']

        # --- read observations
        def read_obs(row):
            f_inp = os.path.join(self.val_obs_dir, row['f_name'])

            if not os.path.isfile(f_inp):
                print('\tNo observation file : ', f_inp)
                return None

            is_nies = row['type'] == 'nies'
            try:
                if is_nies:
                    df = pd.read_csv(f_inp, header=None)
                else:
                    df = pd.read_csv(f_inp, sep='\t')
            except pd.errors.EmptyDataError:
                print('\tEmpty observation file : ', f_inp)
                return None

            # - name the NIES columns and add station metadata
            if is_nies:
                if df.shape[1] != 5:
                    raise ValueError(f'Expected 5 columns in NIES surface file: {f_inp}')

                df.columns = ['year', 'month', 'day', 'hour', 'value']
                df['minute'] = 0
                df['altitude'] = pd.to_numeric(row['alt'], errors='coerce')

            # - clear missing values
            cols = ['year', 'month', 'day', 'hour', 'minute', 'value', 'altitude']

            for col in cols:
                df[col] = pd.to_numeric(df[col], errors='coerce')

            df[cols] = df[cols].replace([np.inf, -np.inf], np.nan)
            df[cols] = df[cols].mask(df[cols] <= -999)
            df.loc[df['value'] <= 0, 'value'] = np.nan

            # - check time
            df.loc[~df['hour'].between(0, 23), 'hour'] = np.nan
            df.loc[~df['minute'].between(0, 59), 'minute'] = np.nan
            df = df.dropna(subset=['year', 'month', 'day', 'hour', 'minute', 'value']).copy()
            df['datetime'] = pd.to_datetime(df[['year', 'month', 'day', 'hour', 'minute']], errors='coerce')
            df = df.dropna(subset=['datetime']).copy()

            if df.empty:
                return None

            # - monthly observation means, ppb
            df['month_id'] = df['datetime'].dt.to_period('M')
            df_m = df.groupby('month_id')[['value', 'altitude']].mean()
            df_m.rename(columns={'value': 'o_ch4'}, inplace=True)

            # - use station altitude where observation altitude is missing
            alt = pd.to_numeric(row['alt'], errors='coerce')
            if pd.notna(alt) and np.isfinite(alt) and alt > -999:
                df_m['altitude'] = df_m['altitude'].fillna(alt)
            return df_m

        # --- extract model at station
        def read_model(ds, row, df_o, col):
            lat = float(row['lat'])
            lon = float(row['lon'])

            # - nearest horizontal grid point
            lat_m = np.asarray(ds['latitude'].values, dtype=float)
            lon_m = np.asarray(ds['longitude'].values, dtype=float)
            iy = int(np.abs(lat_m - lat).argmin())
            ix = int(np.abs((lon_m - lon + 180) % 360 - 180).argmin())

            ds_p = ds[model_vars].isel(latitude=iy, longitude=ix).load()
            da = ds_p['ch4_conc'].transpose('time', 'level')

            # - concentrations in ppb
            levels = np.asarray(ds['level'].values, dtype=float)
            values = np.asarray(da.values, dtype=float) * ch4_factor
            values[(values <= 0) | ~np.isfinite(values)] = np.nan

            # - monthly model means
            years = ds_p['time'].dt.year.values
            months = ds_p['time'].dt.month.values
            dates = [f'{int(y):04d}-{int(m):02d}' for y, m in zip(years, months)]
            index = pd.PeriodIndex(dates, freq='M')

            df_m = pd.DataFrame(values, index=index, columns=range(len(levels)))
            df_m = df_m.groupby(level=0).mean()

            if level_type == 'sigma':
                if ds_p['Ps'].dims != ('time',):
                    raise ValueError(f"Expected Ps(time) at station; found {ds_p['Ps'].dims}")

                ps = np.asarray(ds_p['Ps'].values, dtype=float) * ps_factor
                ps[(ps <= 0) | ~np.isfinite(ps)] = np.nan
                ps_m = pd.Series(ps, index=index).groupby(level=0).mean()
            elif level_type == 'hPa':
                levels = levels * 100.0

            result = pd.Series(np.nan, index=df_o.index, name=col, dtype=float)

            for month in df_o.index.intersection(df_m.index):
                alt = df_o.loc[month, 'altitude']

                if not np.isfinite(alt) or alt <= -999 or alt > 11000:
                    continue

                # - approximate station pressure from standard atmosphere, Pa
                p_obs = 101325.0 * (1.0 - 2.25577e-5 * alt) ** 5.25588

                # - model pressure, Pa
                p_mod = levels * ps_m.loc[month] if level_type == 'sigma' else levels

                ch4 = df_m.loc[month].to_numpy(dtype=float)
                mask = np.isfinite(p_mod) & (p_mod > 0) & np.isfinite(ch4)

                if not mask.any():
                    continue

                p_mod = p_mod[mask]
                ch4 = ch4[mask]
                order = np.argsort(p_mod)
                p_mod = p_mod[order]
                ch4 = ch4[order]

                # - interpolate in log pressure; use boundary level outside range
                result.loc[month] = np.interp(np.log(p_obs), np.log(p_mod), ch4)
            return result

        # - open models
        datasets = {}
        results = []

        with ExitStack() as stack:
            for col, case in self.val_model_cases.items():
                f_inp = case['file']
                print(f'\nRead model {case["label"]} : {f_inp}')
                ds = stack.enter_context(xr.open_dataset(f_inp, engine='h5netcdf'))

                rename_map = {'x': 'longitude', 'y': 'latitude', 'plev': 'level'}
                for dim, coord in rename_map.items():
                    if dim not in ds.dims:
                        continue
                    if coord in ds.variables:
                        ds = ds.swap_dims({dim: coord})
                    else:
                        ds = ds.rename({dim: coord})

                missing_vars = sorted(set(model_vars) - set(ds.variables))
                if missing_vars:
                    raise ValueError(f'{f_inp}: missing model variables {", ".join(missing_vars)}')

                if set(ds['ch4_conc'].dims) != {'time', 'level', 'latitude', 'longitude'}:
                    raise ValueError(f"{f_inp}: unexpected ch4_conc dimensions {ds['ch4_conc'].dims}")

                for coord in ('time', 'latitude', 'longitude', 'level'):
                    if coord not in ds.coords or ds[coord].dims != (coord,) or ds.sizes[coord] == 0:
                        raise ValueError(f'{f_inp}: expected a non-empty 1D coordinate {coord}')

                for coord in ('latitude', 'longitude', 'level'):
                    if not np.isfinite(ds[coord].values.astype(float)).all():
                        raise ValueError(f'{f_inp}: {coord} values must be finite')

                levels = np.asarray(ds['level'].values, dtype=float)

                if (levels <= 0).any():
                    raise ValueError(f'{f_inp}: model levels must be positive')

                if level_type == 'sigma' and levels.max() > 1.01:
                    raise ValueError(f'{f_inp}: sigma levels exceed 1.01. Check pressure units or hybrid coordinates.')

                print('CH4 units : ', ds['ch4_conc'].attrs.get('units', 'not specified'))
                if level_type == 'sigma':
                    print('Ps units  : ', ds['Ps'].attrs.get('units', 'not specified'))
                print('Period    : ', ds['time'].values[0], ds['time'].values[-1])
                datasets[col] = ds

            # - process stations
            for j, row in df_st.iterrows():
                print('\nProcessing : ', j, row['site_id'], row['f_name'])
                df_o = read_obs(row)

                if df_o is None or df_o.empty:
                    continue

                df_m = df_o.copy()

                for col, ds in datasets.items():
                    df_m[col] = read_model(ds, row, df_o, col)

                # - use identical valid months for all model cases
                df_m = df_m.replace([np.inf, -np.inf], np.nan)
                df_m = df_m.dropna(subset=['o_ch4'] + target_cols)

                # - drop the whole month if the observation exceeds either model by more than the threshold
                # obs_minus_model = df_m[target_cols].rsub(df_m['o_ch4'], axis='index')
                # drop_rows = obs_minus_model.abs().gt(self.val_obs_diff_threshold).any(axis='columns')
                # df_m = df_m.loc[~drop_rows].copy()
                df_m = df_m.loc[df_m['o_ch4'] >= 1650].copy()

                print('\tCommon months : ', len(df_m))

                if len(df_m) < min_months:
                    print('\tSkipped: insufficient common months')
                    continue

                results.append((row, df_m))

        return results

    # - look up the requested site IDs in the station inventory
    def _get_sites(self, site_ids):
        df = pd.read_csv(self.val_sts, dtype={'site_id': str})
        df['site_id'] = df['site_id'].str.strip()
        df = df.loc[df['site_id'].isin(site_ids)].copy()

        # - fixed surface and tower stations
        df['project'] = df['project'].fillna('').astype(str)
        df_st = df[df['project'].str.contains('surface|tower', case=False, regex=True)].copy()

        for col in ['lat', 'lon']:
            df_st[col] = pd.to_numeric(df_st[col], errors='coerce')

        df_st = df_st[df_st['lat'].between(-90, 90) & df_st['lon'].between(-180, 180)].copy()
        df_st = df_st.dropna(subset=['f_name'])
        df_st = df_st.drop_duplicates(subset=['site_id', 'f_name'])

        print('Fixed station records              :: ', len(df_st))

        df_st = df_st.set_index('site_id')
        ordered_ids = [site_id for site_id in site_ids if site_id in df_st.index]
        return df_st.loc[ordered_ids].reset_index()

    # - restore missing calendar months as gaps and use timestamps for plotting
    def _monthly_time_frame(self, df):
        months = pd.period_range(df.index.min(), df.index.max(), freq='M')
        df = df.reindex(months)
        df.index = months.to_timestamp()
        return df
