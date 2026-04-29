# -*- coding: utf-8 -*-
"""
Finalized:  26/xx/xx
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import os
import netCDF4
import pandas as pd
from netCDF4 import date2num
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

import _set_case


class ValidComp(_set_case.SetCase):

    def __init__(self):
        super().__init__()

        # -
        # self.run_1()
        # -
        self.run_2()

    # ---
    def run_2(self):
        dirs = [
            ['../trout_MIROC/CH4_t42l67_CYC_H_251014/', '../trout_MIROC/CH4_t42l67_CYC_H_post_260323/'],
            ['../trout_MIROC/CH4_t42l67_INCA_H_260306/', '../trout_MIROC/CH4_t42l67_INCA_H_post_260324/']
        ]

        f_name = 'ch4_hip_aircraft-insitu_59_allvalid_out.csv'

        # --- read reference
        df_ref = pd.read_csv(os.path.join(dirs[0][0], f_name))
        df_ref['dateTime'] = pd.to_datetime(df_ref['dateTime'])

        bins = np.arange(df_ref['altitude'].min(), df_ref['altitude'].max() + 500, 500)

        def bin_profile(df, var):
            df = df.copy()
            df['alt_bin'] = pd.cut(df['altitude'], bins)
            grouped = df.groupby('alt_bin', observed=True)[var]

            mean = grouped.mean()
            rmse = np.sqrt(grouped.apply(lambda x: np.mean((x - x.mean()) ** 2)))

            centers = [interval.mid for interval in mean.index]
            return centers, mean.values, rmse.values

        # --- get years
        years = sorted(df_ref['dateTime'].dt.year.unique())

        # --- create subplots
        n = len(years)
        fig, axes = plt.subplots(1, n, figsize=(5 * n, 6), sharey=True)

        if n == 1:
            axes = [axes]

        # --- loop years
        for ax, year in zip(axes, years):

            ch4_o_plotted = False

            for d_apr, d_pst in dirs:

                df_apr = pd.read_csv(os.path.join(d_apr, f_name))
                df_pst = pd.read_csv(os.path.join(d_pst, f_name))

                df_apr['dateTime'] = pd.to_datetime(df_apr['dateTime'])
                df_pst['dateTime'] = pd.to_datetime(df_pst['dateTime'])

                # --- filter by year
                df_apr_y = df_apr[df_apr['dateTime'].dt.year == year]
                df_pst_y = df_pst[df_pst['dateTime'].dt.year == year]

                if df_apr_y.empty:
                    continue

                name_apr = os.path.basename(os.path.normpath(d_apr)).split('_')[2]
                name_pst = name_apr + '_post'

                # --- OBS once per subplot
                if not ch4_o_plotted:
                    z_o, ch4_o, err_o = bin_profile(df_apr_y, 'ch4_o')
                    ax.errorbar(ch4_o, z_o, xerr=err_o,
                                color='black', linewidth=2, label='OBS', fmt='-')
                    ch4_o_plotted = True

                # --- model
                z_apr, ch4_apr, err_apr = bin_profile(df_apr_y, 'ch4_a')
                z_pst, ch4_pst, err_pst = bin_profile(df_pst_y, 'ch4_a')

                ax.errorbar(ch4_apr, z_apr, xerr=err_apr,
                            linestyle='-', label=name_apr, fmt='o', alpha=0.7)

                ax.errorbar(ch4_pst, z_pst, xerr=err_pst,
                            linestyle='--', label=name_pst, fmt='s', alpha=0.7)

            ax.set_title(f'{year}')
            ax.grid()

        axes[0].set_ylabel('Altitude (m)')
        for ax in axes:
            ax.set_xlabel('CH₄ (ppb)')

        axes[-1].legend(loc='best')

        plt.suptitle('CH₄ profiles by year (binned, 500 m, RMSE)')
        plt.tight_layout()
        # plt.show()
        plt.savefig(self.plt_dir + 'check_hippo.png', format='png', dpi=600, bbox_inches='tight')

    # ---
    def run_1(self):

        # ---
        def valis_comp(apr_csv_, pst_csv_, id_):

            # ---
            def read_comp_csv():
                def get_1file(ifile_: str) -> pd.DataFrame:
                    df = pd.read_csv(ifile_)  # , sep=r'\s+', skiprows=0)
                    return df

                ifile = apr_csv_
                df_apr = get_1file(ifile)
                len_apr = len(df_apr)
                # print(df_apr.head())

                ifile = pst_csv_
                df_pst = get_1file(ifile)
                len_pst = len(df_pst)
                # print(df_pst.head())

                return df_apr, df_pst

            # ---
            def check_1_comp(apr, pst):
                # use common altitude range for both datasets
                zmin = min(apr['altitude(m)'].min(), pst['altitude(m)'].min())
                zmax = max(apr['altitude(m)'].max(), pst['altitude(m)'].max())

                bins = np.arange(zmin, zmax + 100, 100)

                def compute_rmse(df):
                    df = df.copy()
                    df['alt_bin'] = pd.cut(df['altitude(m)'], bins)
                    grouped = df.groupby('alt_bin', observed=True)['ch4_a']

                    # RMSE around bin mean (i.e. variability)
                    rmse = np.sqrt(grouped.apply(lambda x: np.mean((x - x.mean()) ** 2)))

                    # bin centers
                    centers = np.array([interval.mid for interval in rmse.index])

                    return centers, rmse.values

                # compute for both
                z_apr, rmse_apr = compute_rmse(apr)
                z_pst, rmse_pst = compute_rmse(pst)

                # plot
                plt.plot(rmse_apr, z_apr, label='APR', marker='o')
                plt.plot(rmse_pst, z_pst, label='PST', marker='s')

                plt.xlabel('RMSE CH₄ (ppb)')
                plt.ylabel('Altitude (m)')
                plt.title('RMSE of CH₄ vs Altitude')
                plt.grid()
                plt.legend()
                plt.show()
                # exit()

            # ===========================================================
            print(f'\n\t\tRun read_comp_txt for: {apr_csv_, pst_csv_, id_}')
            for j1, invc in enumerate(self.invcases[:]):
                a_conc, p_conc = read_comp_csv()
                check_1_comp(a_conc, p_conc)

        # ---
        def get_csv_files(dir_):
            return [os.path.join(dir_, f) for f in os.listdir(dir_) if f.lower().endswith('.csv')]

        # ===========================================================
        self.invcases = ['inv1']
        dirs = [
            ['../trout_MIROC/CH4_t42l67_CYC_H_251014/', '../trout_MIROC/CH4_t42l67_CYC_H_post_260323/'],
            ['../trout_MIROC/CH4_t42l67_INCA_H_260306/', '../trout_MIROC/CH4_t42l67_INCA_H_post_260324/']
        ]

        for ddir in dirs:
            print('')
            print("\n List of apr csv files :")
            apr_csv_files = get_csv_files(ddir[0])
            apr_csv_files = [f for f in apr_csv_files if 'station' in f]
            print(",\n".join(f'\t\t{f}' for f in apr_csv_files))

            print("\n List of pst csv files :")
            pst_csv_files = get_csv_files(ddir[1])
            pst_csv_files = [f for f in pst_csv_files if 'station' in f]
            print(",\n".join(f'\t\t{f}' for f in pst_csv_files))

            # --- safer pairing (recommended)
            apr_csv_files.sort()
            pst_csv_files.sort()

            for apr_csv, pst_csv in zip(apr_csv_files, pst_csv_files):
                filename = os.path.basename(str(apr_csv))
                cmp_id = filename.split("_")[0]

                valis_comp(apr_csv, pst_csv, cmp_id)
