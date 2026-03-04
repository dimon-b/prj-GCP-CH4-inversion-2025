# -*- coding: utf-8 -*-
"""
Created:    22/05/23 17:48
Project:    GCP-21 get statistics about GT, Loss, Loss correction
@author:    Dmitry Belikov
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import _set_case


class LossCorrPlot(_set_case.SetCase):

    def __init__(self):
        super().__init__()

        self.fs = 10
        self.f_cgt = self.inv_pst_dir + 'check_GT_'
        self.f_bldd = self.inv_lsc_dir + 'gcp_burd_add_'

        # --- get loss and burden plots
        self.run_GTLLC()

    def run_GTLLC(self):
        # -
        def one_plot():
            plt.rc('font', family='serif', size=self.fs)
            fig, ax = plt.subplots(figsize=(6, 7), nrows=5, ncols=1, sharex=True)
            print('\t run_GTLLC for', unp, unx)

            # --- read GT file
            df_gt = pd.read_csv(self.f_cgt + unp + '_' + unx + '.txt', sep=r'\s+', skiprows=0)
            print(df_gt)
            df_gt = df_gt[df_gt['year'] < 2025]

            # --- plots
            for j, inv in enumerate(self.invcases[:]):
                df_x = df_gt.loc[df_gt['invc'] == inv]

                # --- combined fluxes
                ax[0].plot(df_x['year'], df_x['prior'], color=colors[j], ls=':', label=inv + ' prior')
                ax[0].plot(df_x['year'], df_x['post'], color=colors[j + 1], ls='--', label=inv + ' post')
                ax[0].plot(df_x['year'], df_x['post*lc_f'], color=colors[j + 2], ls='-', label=inv + ' post×LC')
                ax[1].plot(df_x['year'], df_x['flxc'], color=colors[j],  ls='-', marker='o', ms=0, label=inv + ' LC')

                # --- read additional burden file
                df_brdd = pd.read_csv(self.f_bldd + inv + '.txt', sep='\t', skiprows=0)
                print(df_brdd)

                ax[2].plot(df_brdd['year'], df_brdd['d_ch4'], color=colors[j], label=inv)
                ax[3].plot(df_brdd['year'], df_brdd['LC_fact'], color=colors[j], label=inv)
                ax[4].plot(df_brdd['year'], df_brdd['ch4_ref'], color=colors[j], label='ch4_ref')
                ax[4].plot(df_brdd['year'], df_brdd['ch4_obs'], color=colors[j + 2], ls='--', label='ch4_obs')

            # --- formatting
            ax[0].set_ylabel('CH$_4$ flux', fontsize=8)
            ax[0].set_ylim([450, 650])
            ax[0].legend(loc='upper left', ncol=3, fontsize=8)
            ax[1].set_ylabel('Loss correction flux', fontsize=8)
            ax[2].set_ylabel(r'$\Delta$CH$_4$ at ref. site (SPO)', fontsize=8)
            ax[3].set_ylabel('LC_fact', fontsize=8)
            ax[4].legend(loc='upper left', ncol=3, fontsize=8)

            for ax1 in ax:
                ax1.grid(linestyle=':')
                ax1.set_xlabel('')

            major_ticks = np.arange(self.years[0], self.years[1] + 1, 5)
            minor_ticks = np.arange(self.years[0], self.years[1] + 1, 1)

            ax[-1].set_xticks(major_ticks)
            ax[-1].set_xticks(minor_ticks, minor=True)
            ax[-1].set_xlim([self.years[0], self.years[1] + 1])
            ax[-1].set_xlabel('Year')

            plt.tight_layout()
            # plt.show()
            plt.savefig(self.plt_dir + 'chk_GT_' + self.hcase + '_' + unp + '_' + unx + '.png',
                        format='png', dpi=600, bbox_inches='tight')

        # =======================================
        cmap = plt.get_cmap('jet')
        colors = [cmap(i) for i in np.linspace(0, 1, len(self.invcases) + 2)]

        for unp in self.unpcases[:1]:
            for unx in self.unxcases[:1]:
                one_plot()
