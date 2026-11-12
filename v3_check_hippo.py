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
import matplotlib.colors as mcolors

import _set_case


class CheckHippo(_set_case.SetCase):

    def __init__(self):
        super().__init__()

        # # - separate profiles
        # self.val2_hippo_profiles()
        #
        # # - yearly profiles
        # self.val2_hippo_profiles_yr()
        #
        # # - rmse profiles
        # self.val2_hippo_rmse()
        #
        # # - scatter
        # self.val2_hippo_scatter()
        # self.val2_hippo_scatter_post()

    # ---
    def val2_hippo_scatter_post(self):

        paths = ['../trout_MIROC/CH4_t42l67_CYC_H_post_260323/',
                 '../trout_MIROC/CH4_t42l67_INCA_H_post_260324/'
                 ]
        labels = ['Post CYC', 'Post INCA']
        f_name = 'ch4_hip_aircraft-insitu_59_allvalid_out.csv'

        # - campaign periods, exclusive stop dates
        campaigns = [['HIPPO-1 Jan 2009', '2009-01-01', '2009-02-01'],
                     ['HIPPO-2 Nov 2009', '2009-11-01', '2009-12-01'],
                     ['HIPPO-3 Mar–Apr 2010', '2010-03-01', '2010-05-01'],
                     ['HIPPO-4 Jun–Jul 2011', '2011-06-01', '2011-08-01'],
                     ['HIPPO-5 Aug–Sep 2011', '2011-08-01', '2011-10-01']
                     ]

        # - read files once
        dfs = []

        for path, label in zip(paths, labels):
            f_inp = os.path.join(path, f_name)
            df = pd.read_csv(f_inp)
            df['dateTime'] = pd.to_datetime(df['dateTime'], errors='coerce', utc=True).dt.tz_localize(None)

            cols = ['latitude', 'altitude', 'ch4_o', 'ch4_a']
            for col in cols:
                df[col] = pd.to_numeric(df[col], errors='coerce')

            df[cols] = df[cols].replace([np.inf, -np.inf], np.nan)
            df = df.dropna(subset=['dateTime'] + cols).copy()
            df = df[
                df['latitude'].between(-90, 90) & (df['altitude'] >= 0) & (df['ch4_o'] > 0) & (df['ch4_a'] > 0)].copy()

            # - altitude in km, concentrations already in ppb
            df['alt_km'] = df['altitude'] / 1000.0
            df['diff'] = df['ch4_a'] - df['ch4_o']

            # - retain plotted campaign periods
            mask = pd.Series(False, index=df.index)
            for title, start, stop in campaigns:
                mask = mask | ((df['dateTime'] >= start) & (df['dateTime'] < stop))

            df = df[mask].copy()
            dfs.append(df)
            print(label, ' plotted observations : ', len(df))

        if all(df.empty for df in dfs):
            print('No valid observations in the selected campaign periods')
            return

        # - common altitude range
        zmax = max(df['alt_km'].max() for df in dfs if not df.empty)
        zmax = max(1.0, np.ceil(zmax))

        # - common discrete colour scale
        bounds = np.arange(-60, 61, 10)
        cmap = plt.get_cmap('PuOr_r', len(bounds) - 1)
        norm = mcolors.BoundaryNorm(bounds, cmap.N)

        # - create subplots
        fig, axes = plt.subplots(2, 5, figsize=(13, 6), sharex=True, sharey=True)
        fig.subplots_adjust(left=0.07, right=0.94, bottom=0.20, top=0.87, wspace=0.12, hspace=0.24)

        # - loop model cases and campaigns
        for i, (df, label) in enumerate(zip(dfs, labels)):
            for j, (title, start, stop) in enumerate(campaigns):
                ax = axes[i, j]
                df_y = df[(df['dateTime'] >= start) & (df['dateTime'] < stop)].copy()

                ax.set_title(f'{chr(97 + i * 5 + j)})', fontsize=12)
                ax.set_xlim([-90, 90])
                ax.set_ylim([0, zmax])
                ax.set_yticks(np.arange(0, 16, 3))
                ax.set_xticks(np.arange(-90, 91, 30))
                ax.tick_params(axis='both', labelsize=10)
                ax.set_axisbelow(True)
                ax.grid(linestyle=':', linewidth=0.5)

                if i == 0:
                    ax.text(0.5, 1.16, title, transform=ax.transAxes, ha='center', va='bottom', fontsize=12)

                if j == 4:
                    ax.text(1.08, 0.5, label, transform=ax.transAxes, rotation=270, ha='center', va='center',
                            fontsize=12)

                if df_y.empty:
                    ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center', va='center')
                    continue

                # - scatter of model minus observations
                ax.scatter(df_y['latitude'], df_y['alt_km'], c=df_y['diff'], cmap=cmap, norm=norm, marker='s', s=9,
                           edgecolors='none', rasterized=True)

                # - statistics for plotted observations
                bias = df_y['diff'].mean()
                sd = df_y['diff'].std(ddof=0)
                ax.text(0.03, 0.03, f'Bias = {bias:.2f}\nSD = {sd:.2f}', transform=ax.transAxes, ha='left', va='bottom',
                        fontsize=10, bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))
                print(f'{label:12s} {title:25s} N = {len(df_y):7d}, Bias = {bias:.2f}, SD = {sd:.2f}')

        # - common axis labels
        # fig.text(0.5, 0.13, 'Latitude', ha='center', fontsize=13)
        fig.text(0.025, 0.54, 'Altitude (km)', va='center', rotation=90, fontsize=13)

        # - common colour bar
        cax = fig.add_axes([0.24, 0.095, 0.52, 0.025])
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        cb = fig.colorbar(sm, cax=cax, orientation='horizontal', ticks=bounds, extend='both')
        cb.set_label('Model − observations CH₄ (ppb)', fontsize=12)
        cb.ax.tick_params(labelsize=10)

        # - save plot
        os.makedirs(self.plt_dir, exist_ok=True)
        f_out = os.path.join(self.plt_dir, 'val2_hippo_scatter_post.png')
        fig.savefig(f_out, format='png', dpi=600, bbox_inches='tight')
        print('\nPlot saved to : ', f_out)

        plt.show()
        plt.close(fig)

    # ---
    def val2_hippo_scatter(self):

        dirs = [
            ['../trout_MIROC/CH4_t42l67_CYC_H_251014/', '../trout_MIROC/CH4_t42l67_CYC_H_post_260323/'],
            ['../trout_MIROC/CH4_t42l67_INCA_H_260306/', '../trout_MIROC/CH4_t42l67_INCA_H_post_260324/']
        ]

        f_name = 'ch4_hip_aircraft-insitu_59_allvalid_out.csv'

        # - row order
        paths = [dirs[0][0], dirs[1][0], dirs[0][1], dirs[1][1]]
        labels = ['Prior CYC', 'Prior INCA', 'Post CYC', 'Post INCA']

        # - campaign periods, exclusive stop dates
        campaigns = [
            ['HIPPO-1 Jan 2009', '2009-01-01', '2009-02-01'],
            ['HIPPO-2 Nov 2009', '2009-11-01', '2009-12-01'],
            ['HIPPO-3 Mar–Apr 2010', '2010-03-01', '2010-05-01'],
            ['HIPPO-4 Jun–Jul 2011', '2011-06-01', '2011-08-01'],
            ['HIPPO-5 Aug–Sep 2011', '2011-08-01', '2011-10-01']
        ]

        # - read files once
        dfs = []

        for path, label in zip(paths, labels):
            f_inp = os.path.join(path, f_name)
            df = pd.read_csv(f_inp)
            df['dateTime'] = pd.to_datetime(df['dateTime'], errors='coerce', utc=True).dt.tz_localize(None)

            cols = ['latitude', 'altitude', 'ch4_o', 'ch4_a']
            for col in cols:
                df[col] = pd.to_numeric(df[col], errors='coerce')

            df[cols] = df[cols].replace([np.inf, -np.inf], np.nan)
            df = df.dropna(subset=['dateTime'] + cols).copy()
            df = df[
                df['latitude'].between(-90, 90) & (df['altitude'] >= 0) & (df['ch4_o'] > 0) & (df['ch4_a'] > 0)].copy()

            # - altitude in km, concentrations already in ppb
            df['alt_km'] = df['altitude'] / 1000.0
            df['diff'] = df['ch4_a'] - df['ch4_o']

            # - retain plotted campaign periods
            mask = pd.Series(False, index=df.index)
            for title, start, stop in campaigns:
                mask = mask | ((df['dateTime'] >= start) & (df['dateTime'] < stop))

            df = df[mask].copy()
            dfs.append(df)
            print(label, ' plotted observations : ', len(df))

        if all(df.empty for df in dfs):
            print('No valid observations in the selected campaign periods')
            return

        # - common altitude range
        zmax = max(df['alt_km'].max() for df in dfs if not df.empty)
        zmax = max(1.0, np.ceil(zmax))

        # - common discrete colour scale
        bounds = np.arange(-60, 61, 10)
        cmap = plt.get_cmap('PuOr_r', len(bounds) - 1)
        norm = mcolors.BoundaryNorm(bounds, cmap.N)

        # - create subplots
        fig, axes = plt.subplots(4, 5, figsize=(18, 13), sharex=True, sharey=True)
        fig.subplots_adjust(left=0.07, right=0.94, bottom=0.15, top=0.92, wspace=0.12, hspace=0.24)

        # - loop model cases and campaigns
        for i, (df, label) in enumerate(zip(dfs, labels)):
            for j, (title, start, stop) in enumerate(campaigns):
                ax = axes[i, j]
                df_y = df[(df['dateTime'] >= start) & (df['dateTime'] < stop)].copy()

                ax.set_title(f'{chr(97 + i * 5 + j)})', fontsize=12)
                ax.set_xlim([-90, 90])
                ax.set_ylim([0, zmax])
                ax.set_xticks(np.arange(-90, 91, 30))
                ax.tick_params(axis='both', labelsize=10)
                ax.set_axisbelow(True)
                ax.grid(linestyle=':', linewidth=0.5)

                if i == 0:
                    ax.text(0.5, 1.16, title, transform=ax.transAxes, ha='center', va='bottom', fontsize=12)

                if j == 4:
                    ax.text(1.08, 0.5, label, transform=ax.transAxes, rotation=270, ha='center', va='center',
                            fontsize=12)

                if df_y.empty:
                    ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center', va='center')
                    continue

                # - scatter of model minus observations
                ax.scatter(df_y['latitude'], df_y['alt_km'], c=df_y['diff'], cmap=cmap, norm=norm, marker='s', s=9,
                           edgecolors='none', rasterized=True)

                # - statistics for plotted observations
                bias = df_y['diff'].mean()
                sd = df_y['diff'].std(ddof=0)
                ax.text(0.03, 0.03, f'Bias = {bias:.2f}\nSD = {sd:.2f}', transform=ax.transAxes, ha='left', va='bottom',
                        fontsize=10, bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))

                print(f'{label:12s} {title:25s} N = {len(df_y):7d}, Bias = {bias:.2f}, SD = {sd:.2f}')

        # - common axis labels
        fig.text(0.5, 0.11, 'Latitude', ha='center', fontsize=13)
        fig.text(0.025, 0.54, 'Altitude (km)', va='center', rotation=90, fontsize=13)

        # - common colour bar
        cax = fig.add_axes([0.24, 0.055, 0.52, 0.018])
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        cb = fig.colorbar(sm, cax=cax, orientation='horizontal', ticks=bounds, extend='both')
        cb.set_label('Model − observations CH₄ (ppb)', fontsize=12)
        cb.ax.tick_params(labelsize=10)

        # - save plot
        os.makedirs(self.plt_dir, exist_ok=True)
        f_out = os.path.join(self.plt_dir, 'val2_hippo_scatter_all.png')
        fig.savefig(f_out, format='png', dpi=600, bbox_inches='tight')
        print('\nPlot saved to : ', f_out)

        plt.show()
        plt.close(fig)


    # ---
    def val2_hippo_rmse(self):
        import os
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D

        paths = [
            '../trout_MIROC/CH4_t42l67_CYC_H_post_260323/',
            '../trout_MIROC/CH4_t42l67_INCA_H_post_260324/'
        ]
        labels = ['Post CYC', 'Post INCA']
        styles = ['-', '--']
        colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
        f_name = 'ch4_hip_aircraft-insitu_59_allvalid_out.csv'

        # - campaign periods, exclusive stop dates
        campaigns = [
            ['HIPPO-1 Jan 2009', '2009-01-01', '2009-02-01'],
            ['HIPPO-2 Nov 2009', '2009-11-01', '2009-12-01'],
            ['HIPPO-3 Mar–Apr 2010', '2010-03-01', '2010-05-01'],
            ['HIPPO-4 Jun–Jul 2011', '2011-06-01', '2011-08-01'],
            ['HIPPO-5 Aug–Sep 2011', '2011-08-01', '2011-10-01']
        ]

        # - altitude bins, 500 m
        bins = np.arange(0, 15001, 500)
        centers = (bins[:-1] + bins[1:]) / 2000.0

        # --- RMSE profile
        def bin_profile(df):
            df = df.copy()
            df['alt_bin'] = pd.cut(df['altitude'], bins, labels=False, include_lowest=True)
            df['diff2'] = df['diff'] ** 2
            mean_diff2 = df.groupby('alt_bin')['diff2'].mean()
            rmse = np.sqrt(mean_diff2.reindex(range(len(centers))))
            return centers, rmse.to_numpy()

        # - read files once
        dfs = []

        for path, label in zip(paths, labels):
            df = pd.read_csv(os.path.join(path, f_name))
            df['dateTime'] = pd.to_datetime(df['dateTime'], errors='coerce', utc=True).dt.tz_localize(None)

            cols = ['latitude', 'altitude', 'ch4_o', 'ch4_a']
            for col in cols:
                df[col] = pd.to_numeric(df[col], errors='coerce')

            df[cols] = df[cols].replace([np.inf, -np.inf], np.nan)
            df = df.dropna(subset=['dateTime'] + cols).copy()
            df = df[df['latitude'].between(-90, 90) & df['altitude'].between(0, 15000) & (df['ch4_o'] > 0) & (df['ch4_a'] > 0)].copy()

            # - model minus observations, ppb
            df['diff'] = df['ch4_a'] - df['ch4_o']
            dfs.append(df)

        # - one panel
        fig, ax = plt.subplots(figsize=(8, 6))
        campaign_handles = []
        plotted_cases = set()

        # - campaign colour; model line style
        for (title, start, stop), color in zip(campaigns, colors):
            plotted = False

            for i, (df, label, style) in enumerate(zip(dfs, labels, styles)):
                df_y = df[(df['dateTime'] >= start) & (df['dateTime'] < stop)].copy()

                if df_y.empty:
                    continue

                z, rmse = bin_profile(df_y)
                ax.plot(rmse, z, color=color, linestyle=style, linewidth=2)
                plotted = True
                plotted_cases.add(i)

                rmse_all = np.sqrt(np.mean(df_y['diff'] ** 2))
                print(f'{label:12s} {title:25s} N = {len(df_y):7d}, RMSE = {rmse_all:.2f} ppb')

            if plotted:
                campaign_handles.append(Line2D([0], [0], color=color, linewidth=2, label=title))

        if not campaign_handles:
            print('No valid observations in the selected campaign periods')
            plt.close(fig)
            return

        # - axes
        ax.set_xlabel('RMSE of model − observations CH₄ (ppb)', fontsize=12)
        ax.set_ylabel('Altitude (km)', fontsize=12)
        ax.set_xlim(left=0)
        ax.set_ylim([0, 15])
        ax.set_yticks(np.arange(0, 16, 3))
        ax.tick_params(axis='both', labelsize=11)
        ax.grid(linestyle=':', linewidth=0.5)

        # - legend in lower right
        model_handles = [Line2D([0], [0], color='black', linestyle=styles[i], linewidth=2, label=labels[i]) for i in sorted(plotted_cases)]
        ax.legend(handles=campaign_handles + model_handles, loc='lower right', ncol=1, fontsize=10, frameon=False)

        # - save plot
        # fig.subplots_adjust(left=0.13, right=0.97, top=0.97, bottom=0.30)
        # os.makedirs(self.plt_dir, exist_ok=True)
        f_out = os.path.join(self.plt_dir, 'val2_hippo_rmse_profiles_post.png')
        fig.savefig(f_out, format='png', dpi=600, bbox_inches='tight')
        print('\nPlot saved to : ', f_out)

        plt.show()
        plt.close(fig)

    # ---
    def val2_hippo_profiles_yr(self):
        dirs = [
            ['../trout_MIROC/CH4_t42l67_CYC_H_251014/', '../trout_MIROC/CH4_t42l67_CYC_H_post_260323/'],
            ['../trout_MIROC/CH4_t42l67_INCA_H_260306/', '../trout_MIROC/CH4_t42l67_INCA_H_post_260324/']
        ]

        f_name = 'ch4_hip_aircraft-insitu_59_allvalid_out.csv'

        # --- read reference
        df_ref = pd.read_csv(os.path.join(dirs[0][0], f_name))
        print(df_ref)
        # exit()
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

                # ax.errorbar(ch4_apr, z_apr, xerr=err_apr, linestyle='-', label=name_apr, fmt='o', alpha=0.7)
                ax.errorbar(ch4_pst, z_pst, xerr=err_pst, linestyle='--', label=name_pst, fmt='s', alpha=0.7)

            ax.set_title(f'{year}')
            ax.grid()

        axes[0].set_ylabel('Altitude (m)')
        for ax in axes:
            ax.set_xlabel('CH₄ (ppb)')

        axes[-1].legend(loc='best')

        plt.suptitle('CH₄ profiles by year (binned, 500 m, RMSE)')
        plt.tight_layout()
        plt.savefig(self.plt_dir + 'val2_hippo_profiles_yr.png', format='png', dpi=600, bbox_inches='tight')
        plt.show()

    # ---
    def val2_hippo_profiles(self):

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
