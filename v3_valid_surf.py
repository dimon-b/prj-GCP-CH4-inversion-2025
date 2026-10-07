# -*- coding: utf-8 -*-
"""
Finalized:  26/xx/xx
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import pandas as pd
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

import _set_case


class ConcValSurf(_set_case.SetCase):

    def __init__(self):
        super().__init__()

        # self.conc3d = ['../results2025/CYC/nc_out/MIROC4-ACTM_CH4conc_GMB_SURF_OH_Transcom.nc',
        #                '../results2025/INCA/nc_out/MIROC4-ACTM_CH4conc_GMB_SURF_OH_INCA.nc']
        self.conc3d = ['../submition/v2/MIROC4-ACTM_CH4conc_GMB_SURF_OH_Transcom.nc',
                       '../submition/v2/MIROC4-ACTM_CH4conc_GMB_SURF_OH_INCA.nc']
        self.df_name = 'model_obs_stats_inv2025.csv'
        self.o_post_dir = '../../prj_apack/obs/o_post_val/'
        self.val_sts = self.o_post_dir + '_vall_all.csv'

        run = 0
        if run:
            # - prepare results
            df_res = self.get_model_obs(ch4_factor=1.0)
            os.makedirs(self.plt_dir, exist_ok=True)
            df_res.to_csv(os.path.join(self.plt_dir, self.df_name), index=False)

        # - print and plot results
        f_inp = os.path.join(self.plt_dir, self.df_name)

        try:
            df_res = pd.read_csv(f_inp)
            print(df_res)
            # - drop slope < 0.5
            df_res = df_res[~(df_res['slp_cyc'] < 0.5)].copy()
            # - drop HBA and USH
            df_res = df_res[~df_res['site'].isin(['HBA', 'USH'])].copy()
            # - drop southern sites with sd_cyc > 20
            # df_res = df_res[~((df_res['lat'] < 5) & (df_res['sd_cyc'] > 20))].copy()
            # df_res = df_res[~((df_res['lat'] < 5) & (df_res['slp_inca'] < 0.85))].copy()
        except pd.errors.EmptyDataError:
            print('\nNo prepared station results')
            return
        self.val3_sites_3panels(df_res)

    # - prepare model-observation results
    def get_model_obs(self, ch4_factor=1.0, level_type='sigma',
                      ps_factor=1.0, min_months=10 * 12):

        # --- read observations
        def read_obs(row):
            f_inp = os.path.join(self.o_post_dir, row['f_name'])

            if not os.path.isfile(f_inp):
                print('\tNo observation file : ', f_inp)
                return None

            # - read NIES surface data
            if row['type'] == 'nies':
                df = pd.read_csv(f_inp, sep=',', header=None)
                if df.shape[1] != 5:
                    raise ValueError(f'Expected 5 columns in NIES surface file: {f_inp}')

                df.columns = ['year', 'month', 'day', 'hour', 'value']
                df['minute'] = 0
                df['altitude'] = pd.to_numeric(row['alt'], errors='coerce')

            # - read converted ObsPack and WDCGG data
            else:
                df = pd.read_csv(f_inp, sep='\t')

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

            ds_p = ds[['ch4_conc', 'Ps']].isel(latitude=iy, longitude=ix).load()
            da = ds_p['ch4_conc'].transpose('time', 'level')

            if set(ds_p['Ps'].dims) != {'time'}:
                raise ValueError(f"Expected Ps(time) at station; found {ds_p['Ps'].dims}")

            # - concentrations in ppb and surface pressure in Pa
            levels = np.asarray(ds['level'].values, dtype=float)
            values = np.asarray(da.values, dtype=float) * ch4_factor
            values[(values <= 0) | ~np.isfinite(values)] = np.nan

            ps = np.asarray(ds_p['Ps'].values, dtype=float) * ps_factor
            ps[(ps <= 0) | ~np.isfinite(ps)] = np.nan

            # - monthly model means
            years = ds_p['time'].dt.year.values
            months = ds_p['time'].dt.month.values
            dates = [f'{int(y):04d}-{int(m):02d}' for y, m in zip(years, months)]
            index = pd.PeriodIndex(dates, freq='M')

            df_m = pd.DataFrame(values, index=index, columns=range(len(levels)))
            df_m = df_m.groupby(level=0).mean()
            ps_m = pd.Series(ps, index=index).groupby(level=0).mean()
            result = pd.Series(np.nan, index=df_o.index, name=col, dtype=float)

            for month in df_o.index.intersection(df_m.index):
                alt = df_o.loc[month, 'altitude']

                if not np.isfinite(alt) or alt <= -999 or alt > 11000:
                    continue

                # - approximate station pressure from standard atmosphere, Pa
                p_obs = 101325.0 * (1.0 - 2.25577e-5 * alt) ** 5.25588

                # - model pressure, Pa
                if level_type == 'sigma':
                    p_mod = levels * ps_m.loc[month]
                elif level_type == 'Pa':
                    p_mod = levels.copy()
                else:
                    p_mod = levels * 100.0

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

        # --- calculate statistics
        def get_corr_sd(df_i):
            stats = {}

            for col in target_cols:
                obs = df_i['o_ch4']
                mod = df_i[col]

                stats['sd_' + col] = np.std(mod - obs, ddof=0)

                if obs.std(ddof=1) > 0:
                    stats['slp_' + col] = obs.cov(mod) / obs.var(ddof=1)
                    stats['corr_' + col] = obs.corr(mod) if mod.std(ddof=1) > 0 else np.nan
                else:
                    stats['slp_' + col] = np.nan
                    stats['corr_' + col] = np.nan

            return stats

        # ===== =====
        if len(self.conc3d) != 2:
            raise ValueError('self.conc3d must contain the CYC and INCA files')

        if level_type not in ['sigma', 'Pa', 'hPa']:
            raise ValueError("level_type must be 'sigma', 'Pa' or 'hPa'")

        if min_months < 2:
            raise ValueError('min_months must be at least 2')

        target_cols = ['cyc', 'inca']

        # - read joined
        df = pd.read_csv(self.val_sts)
        print('Total selected sites               :: ', len(df))

        # - fixed surface and tower stations
        df['project'] = df['project'].fillna('')
        df_st = df[df['project'].str.contains('surface|tower', case=False, regex=True)].copy()

        for col in ['lat', 'lon']:
            df_st[col] = pd.to_numeric(df_st[col], errors='coerce')

        df_st = df_st[df_st['lat'].between(-90, 90) & df_st['lon'].between(-180, 180)].copy()
        df_st = df_st.drop_duplicates(subset=['site_id', 'f_name'])

        print('Fixed station records              :: ', len(df_st))

        # - open models
        datasets = []
        results = []

        try:
            for f_inp in self.conc3d:
                print('\nRead model : ', f_inp)
                ds = xr.open_dataset(f_inp, engine="h5netcdf")
                datasets.append(ds)

                rename_map = {'x': 'longitude', 'y': 'latitude', 'plev': 'level'}
                ds = ds.rename({k: v for k, v in rename_map.items() if k in ds.dims})
                datasets[-1] = ds

                if set(ds['ch4_conc'].dims) != {'time', 'level', 'latitude', 'longitude'}:
                    raise ValueError(f"Unexpected ch4_conc dimensions: {ds['ch4_conc'].dims}")

                for coord in ['latitude', 'longitude', 'level']:
                    if ds[coord].dims != (coord,):
                        raise ValueError(f'Expected a 1D coordinate: {coord}')

                levels = np.asarray(ds['level'].values, dtype=float)

                if not np.isfinite(levels).all() or (levels <= 0).any():
                    raise ValueError('Model levels must be finite and positive')

                if level_type == 'sigma' and levels.max() > 1.01:
                    raise ValueError('Levels exceed 1: check pressure units or hybrid coordinates')

                print('CH4 units : ', ds['ch4_conc'].attrs.get('units', 'not specified'))
                print('Ps units  : ', ds['Ps'].attrs.get('units', 'not specified'))
                print('Period    : ', ds['time'].values[0], ds['time'].values[-1])

            # - process stations
            for j, row in df_st.iterrows():
                print('\nProcessing : ', j, row['site_id'], row['f_name'])
                df_o = read_obs(row)

                if df_o is None or df_o.empty:
                    continue

                df_m = df_o.copy()

                for ds, col in zip(datasets, target_cols):
                    df_m[col] = read_model(ds, row, df_o, col)

                # - use identical valid months for both models
                df_m = df_m.replace([np.inf, -np.inf], np.nan)
                df_m = df_m.dropna(subset=['o_ch4'] + target_cols)
                print('Common months : ', len(df_m))

                if len(df_m) < min_months:
                    print('\tSkipped: insufficient common months')
                    continue

                stats = {'site': row['site_id'], 'lat': row['lat'], 'lon': row['lon'], 'project': row['project'],
                         'f_name': row['f_name'], 'n_months': len(df_m)}
                stats.update(get_corr_sd(df_m))
                results.append(stats)

        finally:
            for ds in datasets:
                ds.close()

        # - return prepared results
        df_res = pd.DataFrame(results)

        if not df_res.empty:
            df_res = df_res.sort_values(by=['lat', 'site']).reset_index(drop=True)

        return df_res

    # - print and plot prepared results
    def val3_sites_3panels(self, df_res):
        if df_res.empty:
            print('\nNo station records meet the common-month requirement')
            return

        target_cols = ['cyc', 'inca']
        target_labs = ['CYC', 'INCA']
        colors = ['b', 'r']
        markers = ['o', '*']

        # - print results
        print('\nTotal compared station records     :: ', len(df_res))
        print('Total distinct station IDs         :: ', df_res['site'].nunique())
        print(df_res.sort_values(by=['lat', 'site']).round(2).to_string(index=False))

        # - print median statistics
        print('\nMedian statistics')
        for col, lab in zip(target_cols, target_labs):
            sd = df_res['sd_' + col].median()
            corr = df_res['corr_' + col].median()
            slp = df_res['slp_' + col].median()
            print(f'{lab:8s} SD = {sd:.2f} ppb, correlation = {corr:.2f}, slope = {slp:.2f}')

        # - plot settings
        fig, axs = plt.subplots(1, 3, figsize=(12, 4))
        metrics = ['sd', 'corr', 'slp']
        labels = ['SD of model bias, ppb', 'Correlation coefficient', 'Slope']
        titles = ['a)', 'b)', 'c)']
        locations = ['upper left', 'lower left', 'lower left']

        # - plot statistics
        for ax, metric, ylabel, title, loc in zip(axs, metrics, labels, titles, locations):
            for col, lab, color, marker in zip(target_cols, target_labs, colors, markers):
                name = metric + '_' + col
                median = df_res[name].median()
                ax.scatter(df_res['lat'], df_res[name], color=color, marker=marker, s=45, edgecolors='black',
                           linewidths=0.5, label=f'{lab} ({median:.2f})')

            ax.set_xlabel('Latitude')
            ax.set_ylabel(ylabel)
            ax.set_title(title)
            ax.set_xlim([-90, 90])
            ax.set_xticks(np.arange(-90, 91, 30))
            ax.grid(linestyle=':')
            ax.legend(loc=loc, fontsize=10)

        axs[0].set_ylim(bottom=0)
        axs[1].set_ylim([-1.05, 1.05])
        axs[2].axhline(1.0, color='gray', linestyle='--', linewidth=0.7)

        # - reference figure limits
        axs[0].set_ylim([0, 90])
        axs[1].set_ylim([0.4, 1.1])
        axs[2].set_ylim([0.4, 1.4])

        plt.tight_layout()
        plt.savefig(self.plt_dir + 'val3_sites_3panels.png', format='png', dpi=600, bbox_inches='tight')
        plt.show()
