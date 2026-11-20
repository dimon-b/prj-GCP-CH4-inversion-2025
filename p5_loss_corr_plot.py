# -*- coding: utf-8 -*-
"""Diagnostic plots for global CH4 flux, burden, and loss correction."""

import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import _set_case
from p6_fortran_tools import check_inputs, check_outputs, full_path


# =============================================================
# === Plot posterior flux and loss correction
class LossCorrPlot(_set_case.SetCase):

    def __init__(self, setup=None, prior='p30', uncertainty='ctl', site='SYO',
                 dry_run=False, overwrite=False):
        if setup is None:
            super().__init__()
        else:
            self.__dict__.update(vars(setup))

        self.fs = 10
        self.f_cgt = os.path.join(self.inv_pst_dir, f'check_GT_{self.icase}_')
        self.f_bldd = os.path.join(self.inv_lsc_dir, f'gcp_burd_add_{self.icase}_{prior}_{uncertainty}_')
        self.site = site.upper()
        self.run_GTLLC(prior=prior, uncertainty=uncertainty,
                       dry_run=dry_run, overwrite=overwrite)

    # --- run_GTLLC
    def run_GTLLC(self, prior='p30', uncertainty='ctl', dry_run=False, overwrite=False):
        f_gt = f'{self.f_cgt}{prior}_{uncertainty}.txt'
        f_burds = [f'{self.f_bldd}{inv}.txt' for inv in self.invcases]
        f_out = os.path.join(self.plt_dir, f'chk_GT_{self.hcase}_{self.icase}_{prior}_{uncertainty}.png')
        print('\n    POST plot inputs   :', f_gt, f_burds)
        print('    POST plot output   :', f_out)
        if dry_run:
            print('    DRY-RUN (no image created)')
            return

        check_inputs([f_gt] + f_burds)
        check_outputs([f_out], overwrite=overwrite)
        full_path(self.plt_dir).mkdir(parents=True, exist_ok=True)

        df_gt = pd.read_csv(f_gt, sep=r'\s+')
        required = ['year', 'prior', 'post', 'post*lc_f', 'invc']
        if any(c not in df_gt.columns for c in required):
            raise ValueError(f'Missing check_GT columns: {required}')
        df_gt = df_gt.loc[df_gt['year'].between(self.yr_s, min(self.yr_e, 2024))].copy()
        if df_gt.empty:
            raise ValueError('No annual global flux data before 2025')
        if 'LC_dif' not in df_gt.columns:
            df_gt['LC_dif'] = df_gt['post*lc_f'] - df_gt['post']

        fig, ax = plt.subplots(figsize=(8, 9), nrows=5, sharex=True)
        colors = plt.get_cmap('tab10').colors
        for j, inv in enumerate(self.invcases):
            df_x = df_gt.loc[df_gt['invc'].astype(str).eq(inv)].sort_values('year')
            if df_x.empty:
                raise ValueError(f'No global flux data for {inv}: {f_gt}')
            color = colors[j % len(colors)]
            ax[0].plot(df_x['year'], df_x['prior'], color='0.4', ls=':', label=f'{inv} prior')
            ax[0].plot(df_x['year'], df_x['post'], color=color, ls='--', label=f'{inv} post')
            ax[0].plot(df_x['year'], df_x['post*lc_f'], color=color, label=f'{inv} post x LC')
            ax[1].plot(df_x['year'], df_x['LC_dif'], color=color, label=inv)

            df_b = pd.read_csv(f_burds[j], sep=r'\s+').sort_values('year')
            df_b = df_b.loc[df_b['year'].between(self.yr_s, min(self.yr_e, 2024))]
            for column in ['d_ch4', 'LC_fact', 'ch4_ref', 'ch4_obs']:
                if column not in df_b.columns:
                    raise ValueError(f'Missing {column} in {f_burds[j]}')
            ax[2].plot(df_b['year'], df_b['d_ch4'], color=color)
            ax[3].plot(df_b['year'], df_b['LC_fact'], color=color)
            ax[4].plot(df_b['year'], df_b['ch4_ref'], color=color, label=f'{inv} ACTM')
            ax[4].plot(df_b['year'], df_b['ch4_obs'], color=color, ls=':', label=f'{inv} OBS')

        ax[0].set_ylabel('Global CH4 flux\n(Tg CH4/yr)')
        ax[1].set_ylabel('LC flux change\n(Tg CH4/yr)')
        ax[2].set_ylabel('ACTM - OBS (ppb)')
        ax[3].set_ylabel('LC factor')
        ax[4].set_ylabel('CH4 concentration (ppb)')
        ax[4].set_xlabel('Year')
        ax[0].legend(fontsize=8, ncol=2)
        ax[4].legend(fontsize=8)
        ax[3].axhline(1.0, ls=':', color='0.4')
        for a in ax:
            a.grid(linestyle=':', alpha=0.7)
        ax[-1].set_xlim(self.yr_s, min(self.yr_e, 2024))
        fig.suptitle(f'{self.hcase} {self.icase} {prior}/{uncertainty}: LC at {self.site} (1999/2025 SPO anchors)', fontsize=11)
        fig.tight_layout()
        fig.savefig(f_out, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print('    POST plot finished :', f_out)
