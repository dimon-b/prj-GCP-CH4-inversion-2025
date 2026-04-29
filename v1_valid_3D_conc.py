# -*- coding: utf-8 -*-
"""
Finalized:  26/xx/xx
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import xarray as xr
import pandas as pd
import codecs, re

from pylab import *

import _set_case


class Valid3dConc(_set_case.SetCase):

    def __init__(self):
        super().__init__()

        self.conc3d = ['../results2025/CYC/nc_out/MIROC4-ACTM_CH4conc_GMB_SURF_OH_Transcom.nc',
                       '../results2025/INCA/nc_out/MIROC4-ACTM_CH4conc_GMB_SURF_OH_INCA.nc']
        # ---
        self.check_GR()

        # ---
        self.check_sites()

    # ===
    def check_sites(self):
        self.years = [1999, 2025]

        # Sites: name, lat, lon, altitude (unused here)
        self.sites = [
            ['ch4_alt_surface-flask_1_representative', 82.4508, -62.5072, 185.0],
            ['ch4_mlo_surface-flask_1_representative', 19.5362, -155.5763, 3397.0],
            ['ch4_spo_surface-flask_1_representative', -89.98, -24.8, 2810.0],
        ]

        # - Read observations
        def read_obs(site_n):
            dir = 'D:/OneDrive - 国立大学法人千葉大学/prj_apack/obs/o_orig/ObsPack/obspack_ch4_1_GLOBALVIEWplus_v7.0_2024-10-29/data/txt/'
            r_txt = dir + site_n + '.txt'
            try:
                with codecs.open(r_txt, encoding='utf8') as inputs:
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
        # plt.show()
        plt.savefig(self.plt_dir + 'check_3sites.png', format='png', dpi=600, bbox_inches='tight')


    # ---
    def check_GR(self):
        # --- Mauna Loa
        MLO_LAT = 19.54
        MLO_LON = -155.58  # will convert to 0–360

        # -
        def extract_mlo(ds):
            da = ds['ch4_conc']
            da = da.isel(level=0)

            lon = (MLO_LON + 360) % 360  # → ~204.42
            da_mlo = da.sel(latitude=MLO_LAT, longitude=lon, method='nearest')

            return da_mlo

        # -
        def yearly_growth_rate(ts):
            # --- ensure datetime (already ok, but safe)
            ts['time'] = xr.decode_cf(ts.to_dataset(name='tmp')).time

            # --- annual mean
            ts_yr = ts.resample(time='YS').mean()  # start-of-year

            # --- growth rate (year-to-year)
            gr = ts_yr.diff('time')

            gr.name = 'yearly_growth_rate'
            return ts_yr, gr

        # =========================
        results = {}

        for f in self.conc3d:
            print(f"\n=== Processing: {f} ===")
            ds = xr.open_dataset(f, engine="h5netcdf")

            # --- extract MLO time series
            ts = extract_mlo(ds)

            # --- compute yearly GR
            ts_yr, gr = yearly_growth_rate(ts)

            results[f] = {'ts': ts_yr, 'gr': gr}

            print("Mean CH4 (ppb):", float(ts_yr.mean().values))
            print("Mean yearly GR (ppb/yr):", float(gr.mean().values))

        # === PLOT
        fig, ax1 = plt.subplots(figsize=(8, 5))
        ax2 = ax1.twinx()

        for f, v in results.items():
            label = f.split('/')[-1].replace('.nc', '')

            ts = v['ts']
            gr = v['gr']

            # --- CH4 (solid)
            ax1.plot(ts['time'], ts, label=f"{label} CH₄")
            # --- GR (dashed)
            ax2.plot(gr['time'], gr, linestyle='--', label=f"{label} GR")

        ax1.set_ylabel('CH₄ (ppb)')
        ax2.set_ylabel('Yearly Growth Rate (ppb/yr)')
        ax1.set_xlabel('Time')
        ax1.set_title('Yearly CH₄ and Growth Rate for MLO')

        ax1.set_ylim(1750, 1990)
        ax1.set_yticks(np.arange(1750, 1990.1, 30))
        ax2.set_ylim(-15, 25)
        ax1.grid(linestyle=':')

        lines = ax1.get_lines() + ax2.get_lines()
        labels = [l.get_label() for l in lines]
        ax1.legend(lines, labels, loc='best', fontsize=8)

        plt.tight_layout()
        plt.savefig(self.plt_dir + 'check_GR.png', format='png', dpi=600, bbox_inches='tight')
        # plt.show()

        self.results = results
