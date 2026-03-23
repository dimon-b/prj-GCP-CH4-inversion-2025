# -*- coding: utf-8 -*-
"""
Finalized:  26/xx/xx
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pylab import *

import _set_case


class ConcValSurf(_set_case.SetCase):

    def __init__(self):
        super().__init__()
        '''Used here dataset prepared by apack project'''

        self.target_cols = ['apr_cyc', 'apr_inc', 'pst_cyc', 'pst_inc']
        self.target_labs = ['Prior CYC', 'Prior INCA', 'Post CYC', 'Post INCA']

        # self.target_cols = ['apr_cyc', 'pst_cyc']
        # self.target_labs = ['Prior CYC', 'Post CYC']

        # self.target_cols = ['apr_inc',  'pst_inc']
        # self.target_labs = [ 'Prior INC',  'Post INC']

        df_j = self.get_4mod_df()
        # df_j = df_j[df_j['year'] < 2009]
        print(df_j)
        df_r = self.get_corr_sd(df_j)
        # self.prt_table(df_r)
        self.plt_corr_sd(df_r)

        self.plot_ts(df_j, 'ALT')
        self.plot_ts(df_j, 'SPO')
        self.plot_ts(df_j, 'BHD')
        self.plot_ts(df_j, 'SMO')
        self.plot_ts(df_j, 'DEM')
        self.plot_ts(df_j, 'MNM')
        self.plot_ts(df_j, 'BRW')
        # self.plt_site_map(df_r)

        # with pd.option_context('display.max_rows', None):
        #     print(df_r.sort_values(by='lat', ascending=False))

    # --- prt_table
    def prt_table(self, df_):

        #
        print('\t Number of sites concidered = ', len(df_['lat'].sort_values().unique()), ' Total = ', len(df_))

        #
        with pd.option_context('display.max_rows', None):
            # print(df_.sort_values(by='site', ascending=True))
            df2 = df_.sort_values(by='site', ascending=True).reset_index(drop=True).round(2)
            df2.index = df2.index + 1
            print(df2)

    # --- plt_corr_sd
    def plt_corr_sd(self, df_):

        print('\t Number of sites concidered = ', len(df_['lat'].sort_values().unique()),
              ' Total = ', len(df_))

        fig, axs = plt.subplots(1, 3, figsize=(12, 4))

        colors = ['lime', 'b', 'r', 'yellow']
        markers = ['s', 'o', '*', '^']

        metrics = {
            'sd': {
                'cols': [f'sd_{c}' for c in self.target_cols],
                'ylabel': 'SD of model bias, ppb',
                'ylim': [0, 70],
                'title': 'a)',
                'legend_loc': 'upper left'
            },
            'corr': {
                'cols': [f'corr_{c}' for c in self.target_cols],
                'ylabel': 'Correlation coefficient',
                'ylim': [0.5, 1.1],
                'title': 'b)',
                'legend_loc': 'lower left'
            },
            'slp': {
                'cols': [f'slp_{c}' for c in self.target_cols],
                'ylabel': 'Slope',
                'ylim': [0.5, 1.4],
                'title': 'c)',
                'legend_loc': 'lower left'
            }
        }

        # ---------- a) SD ----------
        ax = axs[0]
        for col, lab, c, m in zip(metrics['sd']['cols'],
                                  self.target_labs,
                                  colors,
                                  markers):
            mmn = f" ({df_[col].median().round(2)})"
            ax.scatter(df_['lat'], df_[col],
                       color=c, marker=m,
                       edgecolors='black', linewidths=0.5,
                       label=lab + mmn)

        ax.set_xlabel('Latitude')
        ax.set_ylabel(metrics['sd']['ylabel'])
        ax.set_title(metrics['sd']['title'])
        ax.grid(linestyle=':')
        ax.set_xlim([-90, 90])
        ax.set_xticks(np.arange(-90, 90.01, 30))
        # ax.set_ylim(metrics['sd']['ylim'])
        ax.legend(loc=metrics['sd']['legend_loc'], fontsize=self.fs - 3)

        # ---------- Bias reduction ----------
        print('\nBias reduction, %')
        mmn = [df_[c].median() for c in metrics['sd']['cols']]
        print(((mmn[0] - mmn[1]) / ((mmn[0] + mmn[1]) / 2) * 100).round(2))
        # print(((mmn[1] - mmn[3]) / ((mmn[1] + mmn[3]) / 2) * 100).round(2))

        print('\nNumber of values < 30, %')
        for col in metrics['sd']['cols']:
            print(round(((df_[col] < 30).sum() / len(df_)) * 100, 2))

        # ---------- b) Correlation ----------
        ax = axs[1]
        for col, lab, c, m in zip(metrics['corr']['cols'],
                                  self.target_labs,
                                  colors,
                                  markers):
            mmn = f" ({df_[col].median().round(2)})"
            ax.scatter(df_['lat'], df_[col],
                       color=c, marker=m,
                       edgecolors='black', linewidths=0.5,
                       label=lab + mmn)

        ax.set_xlabel('Latitude')
        ax.set_ylabel(metrics['corr']['ylabel'])
        ax.set_title(metrics['corr']['title'])
        ax.grid(linestyle=':')
        ax.set_xlim([-90, 90])
        ax.set_xticks(np.arange(-90, 90.01, 30))
        # ax.set_ylim(metrics['corr']['ylim'])
        ax.legend(loc=metrics['corr']['legend_loc'], fontsize=self.fs - 3)

        print('\nMin nsmallest(10) correlation')
        for col in metrics['corr']['cols']:
            print(df_[col].nsmallest(10).mean().round(2))

        print('\nNumber of values > 0.8, %')
        for col in metrics['corr']['cols']:
            print(round(((df_[col] > 0.8).sum() / len(df_)) * 100, 2))

        # ---------- c) Slope ----------
        ax = axs[2]
        for col, lab, c, m in zip(metrics['slp']['cols'],
                                  self.target_labs,
                                  colors,
                                  markers):
            mmn = f" ({df_[col].median().round(2)})"
            ax.scatter(df_['lat'], df_[col],
                       color=c, marker=m,
                       edgecolors='black', linewidths=0.5,
                       label=lab + mmn)

        ax.set_xlabel('Latitude')
        ax.set_ylabel(metrics['slp']['ylabel'])
        ax.set_title(metrics['slp']['title'])
        ax.grid(linestyle=':')
        ax.set_xlim([-90, 90])
        ax.set_xticks(np.arange(-90, 90.01, 30))
        # ax.set_ylim(metrics['slp']['ylim'])
        ax.legend(loc=metrics['slp']['legend_loc'], fontsize=self.fs - 3)

        plt.tight_layout()
        plt.show()

    def plt_corr_sd_0(self, df_):

        #
        print('\t Number of sites concidered = ', len(df_['lat'].sort_values().unique()), ' Total = ', len(df_))

        # Plot
        fig, axs = plt.subplots(1, 3, figsize=(12, 4))

        # a) Standard Deviation
        ax = axs[0]
        mmn = ' (' + str(df_['sd_pri_cyc'].median().round(2)) + ')'
        ax.scatter(df_['lat'], df_['sd_pri_cyc'], color='lime', marker='s', edgecolors='black', linewidths=0.5,
                   label='Prior CYC' + mmn)
        mmn = ' (' + str(df_['sd_pri_iav'].median().round(2)) + ')'
        ax.scatter(df_['lat'], df_['sd_pri_iav'], color='b', marker='o', edgecolors='black', linewidths=0.5,
                   label='Prior IAV' + mmn)
        mmn = ' (' + str(df_['sd_pst_cyc'].median().round(2)) + ')'
        ax.scatter(df_['lat'], df_['sd_pst_cyc'], color='r', marker='*', edgecolors='black', linewidths=0.5,
                   label='Post CYC' + mmn)
        mmn = ' (' + str(df_['sd_pst_iav'].median().round(2)) + ')'
        ax.scatter(df_['lat'], df_['sd_pst_iav'], color='yellow', marker='^', edgecolors='black', linewidths=0.5,
                   label='Post IAV' + mmn)

        ax.set_xlabel('Latitude')
        ax.set_ylabel('SD of model bias, ppb')
        ax.set_title('a)')
        ax.grid(linestyle=':')
        ax.set_xlim([-90, 90])
        ax.set_xticks(np.arange(-90, 90.01, 30))
        ax.set_ylim([0, 70])
        ax.legend(loc="upper left", ncol=1, fontsize=self.fs - 3)

        # - Bias reduction
        print('\nBias reduction, %')
        mmn1 = df_['sd_pri_cyc'].median()
        mmn2 = df_['sd_pri_iav'].median()
        mmn3 = df_['sd_pst_cyc'].median()
        mmn4 = df_['sd_pst_iav'].median()
        print(((mmn1 - mmn3) / ((mmn1 + mmn3) / 2) * 100).round(2))
        print(((mmn2 - mmn4) / ((mmn2 + mmn4) / 2) * 100).round(2))

        #
        print('\nNumber of values < 30, %')
        print(round(((df_["sd_pri_cyc"] < 30).sum() / len(df_)) * 100, 2))
        print(round(((df_["sd_pri_iav"] < 30).sum() / len(df_)) * 100, 2))
        print(round(((df_["sd_pst_cyc"] < 30).sum() / len(df_)) * 100, 2))
        print(round(((df_["sd_pst_iav"] < 30).sum() / len(df_)) * 100, 2))

        # b) Correlation
        ax = axs[1]
        mmn = ' (' + str(df_['corr_pri_cyc'].median().round(2)) + ')'
        ax.scatter(df_['lat'], df_['corr_pri_cyc'], color='lime', marker='s', edgecolors='black', linewidths=0.5,
                   label='Prior CYC' + mmn)
        mmn = ' (' + str(df_['corr_pri_iav'].median().round(2)) + ')'
        ax.scatter(df_['lat'], df_['corr_pri_iav'], color='b', marker='o', edgecolors='black', linewidths=0.5,
                   label='Prior IAV' + mmn)
        mmn = ' (' + str(df_['corr_pst_cyc'].median().round(2)) + ')'
        ax.scatter(df_['lat'], df_['corr_pst_cyc'], color='r', marker='*', edgecolors='black', linewidths=0.5,
                   label='Post CYC' + mmn)
        mmn = ' (' + str(df_['corr_pst_iav'].median().round(2)) + ')'
        ax.scatter(df_['lat'], df_['corr_pst_iav'], color='yellow', marker='^', edgecolors='black', linewidths=0.5,
                   label='Post IAV' + mmn)
        ax.set_xlabel('Latitude')
        ax.set_ylabel('Correlation coefficient')
        ax.set_title('b)')
        ax.grid(linestyle=':')
        ax.set_xlim([-90, 90])
        ax.set_xticks(np.arange(-90, 90.01, 30))
        ax.set_ylim([0.5, 1.1])
        ax.legend(loc="lower left", ncol=1, fontsize=self.fs - 3)

        # - Min corr
        print('\nMin nsmallest(10) correlation')
        print(df_['corr_pri_cyc'].nsmallest(10).mean().round(2))
        print(df_['corr_pri_iav'].nsmallest(10).mean().round(2))
        print(df_['corr_pst_cyc'].nsmallest(10).mean().round(2))
        print(df_['corr_pst_iav'].nsmallest(10).mean().round(2))

        # -
        print('\nNumber of values > 0.8, %')
        print(round(((df_["corr_pri_cyc"] > 0.8).sum() / len(df_)) * 100, 2))
        print(round(((df_["corr_pri_iav"] > 0.8).sum() / len(df_)) * 100, 2))
        print(round(((df_["corr_pst_cyc"] > 0.8).sum() / len(df_)) * 100, 2))
        print(round(((df_["corr_pst_iav"] > 0.8).sum() / len(df_)) * 100, 2))

        # c) Slope
        ax = axs[2]
        mmn = ' (' + str(df_['slp_pri_cyc'].median().round(2)) + ')'
        ax.scatter(df_['lat'], df_['slp_pri_cyc'], color='lime', marker='s', edgecolors='black', linewidths=0.5,
                   label='Prior CYC' + mmn)
        mmn = ' (' + str(df_['slp_pri_iav'].median().round(2)) + ')'
        ax.scatter(df_['lat'], df_['slp_pri_iav'], color='b', marker='o', edgecolors='black', linewidths=0.5,
                   label='Prior IAV' + mmn)
        mmn = ' (' + str(df_['slp_pst_cyc'].median().round(2)) + ')'
        ax.scatter(df_['lat'], df_['slp_pst_cyc'], color='r', marker='*', edgecolors='black', linewidths=0.5,
                   label='Post CYC' + mmn)
        mmn = ' (' + str(df_['slp_pst_iav'].median().round(2)) + ')'
        ax.scatter(df_['lat'], df_['slp_pst_iav'], color='yellow', marker='^', edgecolors='black', linewidths=0.5,
                   label='Post IAV' + mmn)
        ax.set_xlabel('Latitude')
        ax.set_ylabel('Slope')
        ax.set_title('c)')
        ax.grid(linestyle=':')
        ax.set_xlim([-90, 90])
        ax.set_xticks(np.arange(-90, 90.01, 30))
        ax.set_ylim([0.5, 1.4])
        ax.legend(loc="lower left", ncol=1, fontsize=self.fs - 3)

        plt.tight_layout()
        plt.show()
        # plt.savefig(self.cmp_dir + 'F_2' + '.png', format='png', dpi=600, bbox_inches='tight')

    # --- get_corr_sd
    def get_corr_sd(self, df_):
        correlation_results = []
        grouped = df_.groupby(['site', 'lat', 'lon'])

        list_s = []
        for (site, lat, lon), group in grouped:
            if site in ['SGP']:
                continue
            if lat < 85:
                lat = lat + 2
            stats = {'site': site,
                     'lat': lat,
                     'lon': lon, }
            for col in self.target_cols:
                # Drop rows with NaNs in either column
                valid = group[['o_ch4', col]].dropna()

                len_o = len(valid['o_ch4'])
                if (not valid.empty) & (len_o > 12 * 8):
                    if col == self.target_cols[0]:
                        # print(site, len_o)
                        list_s.append(site)
                    # Pearson correlation
                    corr = valid['o_ch4'].corr(valid[col])
                    # Standard deviation of difference
                    sd_diff = np.std(valid[col] - valid['o_ch4'])
                    # sd_diff = np.mean(valid[col] - valid['o_ch4'])
                    slope = valid[col].std() / valid['o_ch4'].std() * corr
                else:
                    corr = np.nan
                    sd_diff = np.nan
                    slope = np.nan

                if (corr > 0.045) & (sd_diff < 100) & (slope > 0.0) & (lat < 990):
                # if (corr > -0.99) & (sd_diff < 180) & (slope > -0.99) & (lat < 990):
                    stats[f'corr_{col}'] = corr
                    stats[f'sd_{col}'] = sd_diff
                    stats[f'slp_{col}'] = slope
            correlation_results.append(stats)

        # Create DataFrame from results
        df_res = pd.DataFrame(correlation_results)
        df_res.dropna(inplace=True)
        # print(df_res.sort_values(by='corr_pri_cyc', ascending=False))
        # print(df_res[(df_res['slp_pst_cyc'] > df_res['slp_pst_iav']) & (df_res['slp_pst_cyc'] < 1)])

        # Make table
        filtered = df_[df_['site'].isin(list_s)]
        df_1 = filtered.groupby('site', as_index=False).first()
        df_1 = df_1.drop(['year', 'month', 'freq', 'o_ch4'] + self.target_cols, axis=1)
        df_1.rename(
            columns={'site': 'Code', 'type': 'Provider', 'project': 'Type', 'lat': 'Latitude', 'lon': 'Longitude',
                     'alt': 'Altitude'}, inplace=True)

        return df_res

    # --- get_4mod_df
    def get_4mod_df(self):
        invc = 'inv1'
        frms = []
        first_df = None
        for i, (mod, dir_) in enumerate(self.conc_cases):
            print(mod, dir_)
            path = dir_ + '/_surf_mn_' + invc + '.csv'
            df_r = pd.read_csv(path)
            df_r.rename(columns={'m_ch4': mod}, inplace=True)

            if i == 0:
                # Keep first DataFrame with all columns
                first_df = df_r
                # Identify common columns (all except the renamed one)
                common_cols = [col for col in df_r.columns if col != mod]
            else:
                # For subsequent DataFrames, keep only the renamed column
                frms.append(df_r[[mod]])

        # Concatenate first DataFrame with just the new columns from others
        df_b = pd.concat([first_df] + frms, axis=1)

        # clean
        df_b.loc[df_b['o_ch4'] <= -0.0, 'o_ch4'] = np.nan
        df_b.dropna(inplace=True)
        print(df_b)
        return df_b

    # --- plot_ts
    def plot_ts(self, df_, site_):
        # Filter for specific site
        df_site = df_[df_["site"] == site_].copy()

        # Create a datetime column for plotting
        df_site["date"] = pd.to_datetime(df_site[["year", "month"]].assign(day=1))

        # Plotting
        plt.figure(figsize=(10, 6))
        plt.plot(df_site["date"], df_site["o_ch4"], label="Observed (o_ch4)")
        plt.plot(df_site["date"], df_site["apr_cyc"], label="Prior (CYC)")
        plt.plot(df_site["date"], df_site["apr_inc"], label="Prior (AIV)")
        plt.plot(df_site["date"], df_site["pst_cyc"], label="Posterior (CYC)")
        plt.plot(df_site["date"], df_site["pst_inc"], label="Posterior (AIV)")

        plt.title(f"CH₄ Time Series at {site_}")
        plt.xlabel("Date")
        plt.ylabel("CH₄ (ppb)")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    # --- plt_site_map
    def plt_site_map(self, df_):
        # Plotting
        plt.figure(figsize=(8, 4))
        prjct = ccrs.PlateCarree()
        # prjct = ccrs.Robinson()
        ax = plt.axes(projection=prjct)

        # Plot the sites
        ax.scatter(df_['lon'], df_['lat'], color='red', edgecolor='b', s=20, marker='o',
                   transform=ccrs.PlateCarree())

        # Add site labels
        for i, row in df_.iterrows():
            site = row['site']
            xx, yy = 1, 2
            if site in ['BHD', 'YON', 'AMY', 'SHM', 'BIS', 'WLG']:
                xx, yy = -12, 1
            if site in ['MLO', 'UTA', 'RYO', 'LLN', 'VGN', 'EGB', 'LMP', 'JFJ', 'LAU', ]:
                xx, yy = -5, -7
            if site in ['ALT', 'HUN', 'MDN', 'MQA']:
                xx, yy = 1, -5
            if site not in ['PRS', 'CMN', 'TRN', 'SSL', 'ZSF', 'PUY', 'OPE', 'HEI']:
                ax.text(row['lon'] + xx, row['lat'] + yy, row['site'], fontsize=self.fs - 7,
                        transform=ccrs.PlateCarree())

        # self.mpl = [0, 30, 30, 60]
        ax.set_extent(self.mpl)
        ax.set_xticks(np.arange(self.mpl[0], self.mpl[1] + 0.01, 60), crs=ccrs.PlateCarree())
        ax.set_yticks(np.arange(self.mpl[2], self.mpl[3] + 0.01, 30), crs=ccrs.PlateCarree())
        ax.tick_params(axis='both', labelsize=self.fs - 3)
        ax.xaxis.set_label_text('')
        ax.yaxis.set_label_text('')
        ax.grid(linestyle=':', lw=0.5)

        # # Add gridlines with labeled ticks
        # gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), linestyle=':', linewidth=0.5, color='gray')
        # gl.xlocator = plt.FixedLocator(np.arange(self.mpl[0], self.mpl[1] + 0.01, 60))
        # gl.ylocator = plt.FixedLocator(np.arange(self.mpl[2], self.mpl[3] + 0.01, 30))
        # gl.xlabel_style = {'size': 10}
        # gl.ylabel_style = {'size': 10}
        # gl.top_labels = False
        # gl.right_labels = False

        # Add map features
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.OCEAN)
        ax.add_feature(cfeature.COASTLINE, edgecolor='grey')

        # plt.tight_layout()
        # plt.show()
        plt.savefig(self.cmp_dir + 'F_S4' + '.png', format='png', dpi=1200, bbox_inches='tight')
        plt.close()
