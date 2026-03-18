# -*- coding: utf-8 -*-
"""
Created:    21/11/08 17:48
Project:    GCP/GMB 2025 project [server part]
@author:    Dmitry Belikov
"""

import pandas as pd

import _set_case


class LossCorr(_set_case.SetCase):

    def __init__(self):
        super().__init__()
        self.f_blos = self.inv_lsc_dir + 'frt_burd_'
        self.f_bldd = self.inv_lsc_dir + 'gcp_burd_add_'
        self.f_outl = self.inv_lsc_dir + 'gcp_LC_'
        # self.f_obss = '../inp_dir/obs_out/obspack/ch4_mlo_surface-flask_1_representative.txt'
        self.f_obss = '../inv_dir/obs/ch4_spo_surface-flask_1_representative.txt'

        # --- run correction
        self.clc_lcorr_1ref()

    # --- get loss correction coefficient
    def clc_lcorr_1ref(self):

        # - Read burden/loss data
        def get_burdloss():
            fname = self.f_blos + self.tracer + '.txt'
            try:
                df_bl = pd.read_csv(fname, sep=r'\s+')
                # print(df_bl)
                df_bl['day'] = 31
                df_bl['month'] = 12
                df_bl['datetime'] = pd.to_datetime(df_bl[['year', 'month', 'day']])
                df_bl.set_index('datetime', inplace=True)
                df_bl['lch4'] = df_bl['loh'] + df_bl['lcl'] + df_bl['lod']
                df_bl.reset_index(inplace=True)
                df_bl = df_bl.drop(['year', 'month', 'day', 'loh', 'lcl', 'lod'], axis=1)
                return df_bl
            except IOError:
                print('\t\tFile not found', fname)
                exit()

        # - Read obs data
        def get_obs():
            fname = self.f_obss
            try:
                df_ob = pd.read_csv(fname, sep='\t', names=['year', 'month', 'day', 'hour', 'ch4', 'xx', 'yy'])
                df_ob['ch4'] = df_ob['ch4']  # *1e9
                df_ob['datetime'] = pd.to_datetime(df_ob[['year', 'month', 'day', 'hour']])
                df_ob.set_index('datetime', inplace=True)
                df_ob = df_ob.resample('YE').mean()
                df_ob.reset_index(inplace=True)
                df_ob.drop(columns=['month', 'day', 'hour', 'xx', 'yy'], inplace=True)

                # --- add one rec
                df_ob.loc[len(df_ob.index)] = [pd.Timestamp('2025-12-31'), 2025.0, 1883]
                print('CH4 obs from', self.f_obss)
                print(df_ob)
                return df_ob
            except IOError:
                print('\t\tFile not found', fname)
                exit()

        # =======================================
        # --- obs
        df_ob = get_obs()
        df_ob.rename(columns={'ch4': 'ch4_obs'}, inplace=True)

        # -
        for jt, trac in enumerate(self.invcases):
            # - regular tracer
            self.tracer = trac
            df_bl = get_burdloss()

            df_rs = pd.merge(df_bl, df_ob, on='datetime')
            df_rs['d_ch4'] = (df_rs['ch4_ref'] - df_rs['ch4_obs'])
            df_rs['LC'] = df_rs['d_ch4'] / (df_rs['ch4_ref']) * df_rs['lch4']
            df_rs['LC_fact'] = (df_rs['lch4'] - df_rs['LC']) / df_rs['lch4']


            # fixme
            if 'CYC' in self.inv_wrk_dir:
                mask = (df_rs['year'] > 2006) & (df_rs['year'] < 2019)
                df_rs.loc[mask, 'LC_fact'] = (df_rs.loc[mask, 'LC_fact'] - 1.0)*1.5 +1.0
                mask = df_rs['year'] >= 2019
                df_rs.loc[mask, 'LC_fact'] = (df_rs.loc[mask, 'LC_fact'] - 1.0)*2.0 +1.0

            if 'INCA' in self.inv_wrk_dir:
                mask = df_rs['year'] >= 2012
                df_rs.loc[mask, 'LC_fact'] = (df_rs.loc[mask, 'LC_fact'] - 1.0)*1.5 + 1.0

                
            df_rs['bch4_corr'] = df_rs['bch4'] * df_rs['LC_fact']

            # order
            new_order = ["year", "ch4_obs", "ch4_ref", "d_ch4", "bch4", "bch4_corr", "lch4", "LC", "LC_fact", ]
            df_rs = df_rs[new_order]
            print(df_rs)

            # -
            # print('Mean LC, LC_fact, d_ch4, bch4, bch4*LC_fact for', trac,
            #       round(df_rs['LC'].mean(), 2), round(df_rs['LC_fact'].mean(), 4), round(df_rs['d_ch4'].mean(), 4),
            #       round(df_rs['bch4'].mean(), 2), round(df_rs['bch4_corr'].mean(), 2))

            # - to files
            df_rs[['year', 'LC_fact']].to_csv(self.f_outl + self.tracer + '.txt', sep='\t', index=False, header=True)
            df_rs.round(4).to_csv(self.f_bldd + self.tracer + '.txt', sep='\t', index=False, header=True)
