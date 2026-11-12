# -*- coding: utf-8 -*-
"""
Created:    22/05/23 17:48
Project:    GCP-21 get statistics about GT, Loss, Loss correction
@author:    Dmitry Belikov
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import _set_case


class LossCorrPlot(_set_case.SetCase):

    def __init__(self):
        super().__init__()

        self.fs = 10
        self.f_cgt = os.path.join(self.inv_pst_dir, 'check_GT_')
        self.f_bldd = os.path.join(self.inv_lsc_dir, 'gcp_burd_add_')

        # --- get loss and burden plots
        self.run_GTLLC()

    def run_GTLLC(self):
        single_case = len(self.invcases) == 1
        colors = plt.get_cmap('tab10').colors

        # --- plot one uncertainty combination
        def one_plot():
            print('\t run_GTLLC for', unp, unx)

            # - read GT file
            f_gt = f'{self.f_cgt}{unp}_{unx}.txt'
            df_gt = pd.read_csv(f_gt, sep=r'\s+')
            df_gt = df_gt.loc[df_gt['year'] < 2025].copy()

            # - support GT files without the LC_dif column
            if 'LC_dif' not in df_gt.columns:
                df_gt['LC_dif'] = df_gt['post*lc_f'] - df_gt['post']
            print('\n\tFrom ', self.f_cgt)
            print(df_gt)

            plt.rc('font', family='serif', size=self.fs)
            fig, ax = plt.subplots(figsize=(6, 7), nrows=5, ncols=1, sharex=True)

            for j, inv in enumerate(self.invcases):
                df_x = df_gt.loc[df_gt['invc'].eq(inv)].sort_values('year')

                if single_case:
                    prior_color, post_color, lc_color, obs_color = '0.35', 'tab:blue', 'tab:red', 'black'
                else:
                    prior_color = post_color = lc_color = obs_color = colors[j % len(colors)]

                # - global fluxes and the change due to loss correction
                ax[0].plot(df_x['year'], df_x['prior'], color=prior_color, ls=':', label=f'{inv} prior')
                ax[0].plot(df_x['year'], df_x['post'], color=post_color, ls='--', label=f'{inv} post')
                ax[0].plot(df_x['year'], df_x['post*lc_f'], color=lc_color, ls='-', label=f'{inv} post×LC')
                ax[1].plot(df_x['year'], df_x['LC_dif'], color=lc_color)

                # - loss correction diagnostics
                f_brdd = f'{self.f_bldd}{inv}.txt'
                df_brdd = pd.read_csv(f_brdd, sep='\t').sort_values('year')
                print('\n\tFrom ', self.f_bldd)
                print(df_brdd)

                year = df_brdd['year']
                ax[2].plot(year, df_brdd['d_ch4'], color=post_color)
                ax[3].plot(year, df_brdd['LC_fact'], color=lc_color)
                ax[4].plot(year, df_brdd['ch4_ref'], color=post_color, label=f'{inv} model')
                ax[4].plot(year, df_brdd['ch4_obs'], color=obs_color, ls='--', label=f'{inv} obs')

            # - formatting
            ax[0].set_ylabel('CH$_4$ flux\n(Tg yr$^{-1}$)', fontsize=8)
            ax[0].set_ylim([450, 650])
            ax[0].legend(loc='upper left', ncol=3, fontsize=8)
            ax[1].set_ylabel('LC flux change\n(Tg yr$^{-1}$)', fontsize=8)
            ax[2].set_ylabel(r'$\Delta$CH$_4$ at SPO (ppb)', fontsize=8)
            ax[3].set_ylabel('LC factor', fontsize=8)
            ax[4].set_ylabel('CH$_4$ at SPO (ppb)', fontsize=8)
            ax[4].legend(loc='upper left', ncol=3, fontsize=8)

            for ax1 in ax:
                ax1.grid(linestyle=':')

            major_ticks = np.arange(self.years[0]+1, self.years[1] + 1, 4)
            minor_ticks = np.arange(self.years[0], self.years[1] + 1)
            ax[-1].set_xticks(major_ticks)
            ax[-1].set_xticks(minor_ticks, minor=True)
            ax[-1].set_xlim(self.years[0], self.years[1] + 1)
            ax[-1].set_xlabel('Year')

            os.makedirs(self.plt_dir, exist_ok=True)
            f_out = os.path.join(self.plt_dir, f'chk_GT_{self.hcase}_{unp}_{unx}.png')
            fig.tight_layout()
            fig.savefig(f_out, format='png', dpi=600, bbox_inches='tight')
            plt.close(fig)

        # --- first configured uncertainty combination
        for unp in self.unpcases[:1]:
            for unx in self.unxcases[:1]:
                one_plot()
