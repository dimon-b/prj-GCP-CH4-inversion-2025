# -*- coding: utf-8 -*-
"""
Finalized:  26/xx/xx
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import xarray as xr
import pandas as pd
import re

from pylab import *

import _set_case


class Valid3dConc(_set_case.SetCase):

    def __init__(self):
        super().__init__()

        self.conc3d = ['../results2025/CYC/nc_out/MIROC4-ACTM_CH4conc_GMB_SURF_OH_Transcom.nc',
                       '../results2025/INCA/nc_out/MIROC4-ACTM_CH4conc_GMB_SURF_OH_INCA.nc']
        self.conc3d = ['../submition/v2/MIROC4-ACTM_CH4conc_GMB_SURF_OH_Transcom.nc',
                       '../submition/v2/MIROC4-ACTM_CH4conc_GMB_SURF_OH_INCA.nc']
        # ---
        self.check_GR()

        # ---
        self.check_sites()

    # ===
    def check_sites(self):
        self.years = [1999, 2025]

        # Sites: name, lat, lon, altitude (unused here)
        self.sites = [['ch4_alt_surface-flask_1_representative', 82.4508, -62.5072, 185.0],
                      ['ch4_mlo_surface-flask_1_representative', 19.5362, -155.5763, 3397.0],
                      ['ch4_spo_surface-flask_1_representative', -89.98, -24.8, 2810.0],
                      ]

        # - Read observations
        def read_obs(site_n):
            dir = 'D:/OneDrive - 国立大学法人千葉大学/prj_apack/obs/o_orig/ObsPack/obspack_ch4_1_GLOBALVIEWplus_v7.0_2024-10-29/data/txt/'
            r_txt = dir + site_n + '.txt'
            try:
                with open(r_txt, encoding='utf8') as inputs:
                    for i, line in enumerate(inputs):
                        if '# header_lines :' in line:
                            sub_line = re.split("( )", line)
                            skipr = int(sub_line[-1])
                            break
                    df_r = pd.read_csv(r_txt, sep=r'\s+', skiprows=skipr - 1)

                # --- save to output
                df_o = df_r[['year', 'month', 'day', 'hour', 'value']].copy()
                df_o['minute'] = 0
                df_o.rename(columns={'value': 'ch4'}, inplace=True)
                df_o['ch4'] = df_o['ch4'] * 1e9
                df_o['index'] = pd.to_datetime(df_o[['year', 'month', 'day', 'hour', 'minute']])
                df_o.set_index('index', inplace=True)
                df_o.drop(columns=['year', 'month', 'day', 'hour', 'minute'], inplace=True)
                return df_o

            except IOError:
                print('\tNo File: ', r_txt)

        # - Extract model at site
        def extract_model(ds, lat, lon):
            da = ds['ch4_conc']

            # surface
            da = da.isel(level=0)

            # lon fix
            if lon < 0:
                lon = lon + 360

            da = da.sel(latitude=lat, longitude=lon, method='nearest')
            ts = da.to_series()
            return ts.to_frame(name='ch4_model')

        # Load ALL model runs once
        model_data = {}
        for f in self.conc3d:
            ds = xr.open_dataset(f, engine="h5netcdf")
            model_data[f] = ds

        # PLOTTING
        fig, axs = plt.subplots(3, 1, figsize=(6, 6), sharex=True)
        fig.subplots_adjust(hspace=0.3)

        for j, (site_n, lat, lon, lev) in enumerate(self.sites):
            ax = axs[j]

            print(f"\n=== Site: {site_n} ===")

            # --- OBS
            df_o = read_obs(site_n)
            df_oy = df_o.resample('YS').mean()

            # --- plot OBS
            df_o = df_o.loc[pd.Timestamp(str(self.years[0] - 1) + '-10-01'):
                            pd.Timestamp(str(self.years[1] + 1) + '-07-31')]
            ax.plot(df_o.index, df_o['ch4'], color='k', lw=1.0, ls='--', label='Obs')

            color_map = {'MIROC4-ACTM_CH4conc_GMB_SURF_OH_Transcom': 'r',
                         'MIROC4-ACTM_CH4conc_GMB_SURF_OH_INCA': 'lime',
                         }

            # --- MODEL (loop runs)
            for f, ds in model_data.items():
                label = f.split('/')[-1].replace('.nc', '')
                color = color_map.get(label, 'black')  # fallback

                df_m = extract_model(ds, lat, lon)
                df_my = df_m.resample('YS').mean()
                # df = df_my.join(df_oy, how='inner', rsuffix='_obs')
                # print(label)
                # print(df.head())

                # --- plot
                ax.plot(df_m.index, df_m['ch4_model'], color=color, lw=1.0, label=label)

            # --- cosmetics
            ax.set_title(site_n[4:7].upper())
            ax.set_ylabel('CH₄ (ppb)')
            ax.grid()
            ax.set_xlim(pd.Timestamp(str(self.years[0] - 1) + '-10-01'),
                        pd.Timestamp(str(self.years[1] + 1) + '-07-31'))
            # ax.set_ylim(1600, 2100)

        axs[2].legend(loc='best', fontsize=8)
        axs[-1].set_xlabel('Time')

        plt.tight_layout()
        plt.savefig(self.plt_dir + 'val1_check_3sites.png', format='png', dpi=600, bbox_inches='tight')
        plt.show()

    # ---
    def check_GR(self):

        # --- Mauna Loa coordinates
        MLO_LAT = 19.54
        MLO_LON = -155.58

        # ---
        def extract_mlo(ds):
            da = ds['ch4_conc'].isel(level=0)

            lon_vals = da.longitude.values
            if lon_vals.min() < 0:
                lon = MLO_LON
            else:
                lon = (MLO_LON + 360) % 360

            return da.sel(latitude=MLO_LAT, longitude=lon, method='nearest')

        # ---
        def yearly_growth_rate(ts):
            ts = ts.sortby('time')

            # --- smoothing (critical)
            ts_smooth = ts.rolling(time=12, center=True, min_periods=6).mean()

            # --- yearly mean
            ts_yr = ts_smooth.resample(time='YS').mean()

            # --- growth rate
            gr = ts_yr.diff('time')
            gr['time'] = ts_yr['time'][1:]

            ts_yr.name = 'yearly_mean_ch4'
            gr.name = 'yearly_growth_rate'

            return ts_yr, gr

        # ---
        def get_noaa_ch4_gr_obs_df():
            years = np.array(
                [1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000,
                 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017,
                 2018, 2019, 2020, 2021, 2022, 2023, 2024])

            gr_obs = np.array(
                [12.90, 12.18, 12.90, 11.44, 10.72, 11.17, 8.66, 14.00, 2.43, 3.85, 7.27, 3.85, 2.48, 6.33, 12.13, 2.27,
                 -1.34, -0.66, 3.18, 4.83, -4.71, 0.23, 1.84, 7.83, 6.56, 4.74, 5.04, 5.04, 5.00, 5.64, 12.74, 10.00,
                 7.12, 6.84, 8.71, 9.61, 14.85, 17.68, 12.97, 8.56, 7.23])

            gr_unc = np.array(
                [0.80, 0.89, 0.68, 0.70, 0.54, 0.52, 0.44, 0.63, 0.51, 0.57, 0.62, 0.57, 0.37, 0.68, 0.70, 0.70, 0.56,
                 0.47, 0.58, 0.67, 0.45, 0.49, 0.61, 0.60, 0.41, 0.47, 0.72, 0.64, 0.49, 0.58, 0.50, 0.77, 0.69, 0.77,
                 0.54, 0.70, 0.53, 0.36, 0.39, 0.78, 0.56])

            time_obs = pd.to_datetime(years, format='%Y')
            df = pd.DataFrame({'time': time_obs, 'gr_obs': gr_obs, 'gr_unc': gr_unc}).set_index('time')
            return df

        # =========================
        df_obs = get_noaa_ch4_gr_obs_df()
        results = {}
        for f in self.conc3d:
            print(f"\n=== Processing: {f} ===")
            ds = xr.open_dataset(f, engine="h5netcdf")

            ts = extract_mlo(ds)
            ts_yr, gr = yearly_growth_rate(ts)

            results[f] = {'ts': ts_yr, 'gr': gr}

            print("Mean CH4 (ppb):", float(ts_yr.mean().values))
            print("Mean yearly GR (ppb/yr):", float(gr.mean().values))

            # --- metrics vs NOAA
            try:
                gr_model_interp = gr.interp(time=df_obs.index)
                gr_model_series = gr_model_interp.to_series()

                bias = (gr_model_series - df_obs['gr_obs']).mean()
                rmse = np.sqrt(((gr_model_series - df_obs['gr_obs']) ** 2).mean())

                print(f"Bias vs NOAA: {bias:.2f} ppb/yr")
                print(f"RMSE vs NOAA: {rmse:.2f} ppb/yr")
            except Exception as e:
                print("Metrics skipped:", e)

        # --- PLOT
        fig, ax1 = plt.subplots(figsize=(9, 5))

        for f, v in results.items():
            label = f.split('/')[-1].replace('.nc', '')

            # ax1.plot(v['ts']['time'], v['ts'], label=f"{label} CH₄")
            ax1.plot(v['gr']['time'], v['gr'], linestyle='--', label=f"{label} GR")

        # --- NOAA obs
        ax1.plot(df_obs.index, df_obs['gr_obs'], color='black', linewidth=2, label='NOAA GR')

        ax1.fill_between(df_obs.index, df_obs['gr_obs'] - df_obs['gr_unc'],
                         df_obs['gr_obs'] + df_obs['gr_unc'],
                         color='black', alpha=0.15, label='NOAA unc')

        ax1.set_ylabel('Growth Rate (ppb/yr)')
        ax1.set_xlabel('')
        ax1.set_title('CH₄ Growth Rate at MLO (Model vs NOAA)')
        ax1.grid(linestyle=':')

        lines = ax1.get_lines()
        labels = [l.get_label() for l in lines]
        ax1.legend(lines, labels, loc='best', fontsize=12)
        ax1.set_xlim(pd.Timestamp(str(self.years[0] ) + '-10-01'),
                     pd.Timestamp(str(self.years[1] ) + '-07-31'))
        plt.tight_layout()
        plt.savefig(self.plt_dir + 'val1_check_GR.png', format='png', dpi=600, bbox_inches='tight')
        plt.show()
