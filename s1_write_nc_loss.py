# -*- coding: utf-8 -*-
"""
Created:    22/01/26 19:17
Project:    Write total fCH4 *.nc for GCP2021 as originally created by Naveen Negi
@author:    Dmitry Belikov
"""
from datetime import datetime

import netCDF4
import numpy as np
import pandas as pd
from netCDF4 import date2num
import matplotlib.pyplot as plt

# import b_sub_flux
import _set_case


class WriteNcLoss(_set_case.SetCase):

    def __init__(self):
        super().__init__()
        self.invcases = ['inv1', 'inv7', 'vCao']
        self.dirs = ['../runs/pr_gcp_brw/', '../runs/pr_gcp_brw_OH/']  # res_check

        for dir in self.dirs:
            self.wrt_loss_nc(dir)

    def wrt_loss_nc(self, dir):
        print(f'\tRun wrt_loss_nc')
        print(f'\tRun wrt_loss_nc for: {dir}')

        for j1, invc in enumerate(self.invcases[:]):
            a_los, p_los = self.read_loss_txt(dir, invc)
            self.write_1nc(invc, a_los, p_los, dir)

    def read_loss_txt(self, dir, invc):
        def get_1file(ifile):
            df = pd.read_csv(ifile, delim_whitespace=True, skiprows=0)
            print(df)
            return df

        ifile = dir + 'losscorr' + '/' + 'gcp_loss_apr_' + invc + '.txt'
        print(f'\tRead ifile: {ifile}')
        df_apr = get_1file(ifile)
        df_apr.columns = ['year', 'month',
                          'OH_tot_loss_prior', 'Cl_tot_loss_prior', 'O1D_tot_loss_prior',
                          'OH_SH_loss_prior', 'Cl_SH_loss_prior', 'O1D_SH_loss_prior',
                          'OH_NH_loss_prior', 'Cl_NH_loss_prior', 'O1D_NH_loss_prior',
                          'OH_strat_loss_prior', 'Cl_strat_loss_prior', 'O1D_strat_loss_prior']
        print(df_apr)
        df_apr_1 = df_apr.copy()
        df_apr_1['OH_tot_loss_ck'] = df_apr_1['OH_SH_loss_prior'] + df_apr_1['OH_NH_loss_prior']
        df_apr_1[['OH_tot_loss_prior', 'OH_SH_loss_prior', 'OH_NH_loss_prior', 'OH_tot_loss_ck']].plot()
        print(df_apr_1[['OH_tot_loss_prior', 'OH_SH_loss_prior', 'OH_NH_loss_prior', 'OH_tot_loss_ck']])
        plt.show()
        df_apr_1 = df_apr.copy()
        df_apr_1['Cl_tot_loss_ck'] = df_apr_1['Cl_SH_loss_prior'] + df_apr_1['Cl_NH_loss_prior']
        df_apr_1[['Cl_tot_loss_prior', 'Cl_SH_loss_prior', 'Cl_NH_loss_prior', 'Cl_tot_loss_ck']].plot()
        print(df_apr_1[['Cl_tot_loss_prior', 'Cl_SH_loss_prior', 'Cl_NH_loss_prior', 'Cl_tot_loss_ck']])
        plt.show()
        df_apr_1 = df_apr.copy()
        df_apr_1['O1D_tot_loss_ck'] = df_apr_1['O1D_SH_loss_prior'] + df_apr_1['O1D_NH_loss_prior']
        df_apr_1[['O1D_tot_loss_prior', 'O1D_SH_loss_prior', 'O1D_NH_loss_prior', 'O1D_tot_loss_ck']].plot()
        print(df_apr_1[['O1D_tot_loss_prior', 'O1D_SH_loss_prior', 'O1D_NH_loss_prior', 'O1D_tot_loss_ck']])
        plt.show()

        ifile = dir + 'res_check' + '/' + 'gcp_loss_pst_' + invc + '.txt'
        print(f'\tRead ifile: {ifile}')
        df_pst = get_1file(ifile)
        df_pst.columns = ['year', 'month',
                          'OH_tot_loss_post', 'Cl_tot_loss_post', 'O1D_tot_loss_post',
                          'OH_SH_loss_post', 'Cl_SH_loss_post', 'O1D_SH_loss_post',
                          'OH_NH_loss_post', 'Cl_NH_loss_post', 'O1D_NH_loss_post',
                          'OH_strat_loss_post', 'Cl_strat_loss_post', 'O1D_strat_loss_post']

        return df_apr, df_pst

    # ===
    def write_1nc(self, invc, apr, pst, dir):

        fout = dir + 'nc_out' + '/' + '/MIROC4-ACTM_sink_' + invc + 'SURF.nc'
        nc = netCDF4.Dataset(fout, 'w', format='NETCDF4')
        nc.description = 'The prior and posteor CH4 loss resulted from the surface based inversion for GCP-CH4, 2021. The results are produced at JAMSTEC, Japan. For details please contact at prabir@jamstec.go.jp and d.belikov@chiba-u.jp.'
        nc.Institution = 'Japan Agency for Marine-Earth Science and Technology (JAMSTEC)'
        nc.Contact = 'd.belikov@chiba-u.jp and prabir@jamstec.go.jp'
        nc.Disclaimer = 'This data is created for GCP-CH4 2021.'
        nc.History = 'Created on January 2023 by Dmitry Belikov'
        nc.Additional_Info = 'MIROC4.0-based Atmospheric Chemistry-Transport Model (MIROC4-ACTM) is used for forward simulations (Prabir et al., 2018) and Time-independent Bayesian inversion is used for optimising flux (Patra et al., 2016)'

        # --- set dimensions
        time_dim = nc.createDimension('time', None)

        # --- set variables
        time = nc.createVariable('time', 'f4', ('time',))
        time.units = 'hours since 1970-01-01 00:00:00'
        time.calendar = 'gregorian'
        OH_tot_loss_prior = nc.createVariable('OH_tot_loss_prior', 'f4', ('time'))
        OH_SH_loss_prior = nc.createVariable('OH_SH_loss_prior', 'f4', ('time'))
        OH_NH_loss_prior = nc.createVariable('OH_NH_loss_prior', 'f4', ('time'))
        OH_strat_loss_prior = nc.createVariable('OH_strat_loss_prior', 'f4', ('time'))
        Cl_tot_loss_prior = nc.createVariable('Cl_tot_loss_prior', 'f4', ('time'))
        Cl_SH_loss_prior = nc.createVariable('Cl_SH_loss_prior', 'f4', ('time'))
        Cl_NH_loss_prior = nc.createVariable('Cl_NH_loss_prior', 'f4', ('time'))
        Cl_strat_loss_prior = nc.createVariable('Cl_strat_loss_prior', 'f4', ('time'))
        O1D_tot_loss_prior = nc.createVariable('O1D_tot_loss_prior', 'f4', ('time'))
        O1D_SH_loss_prior = nc.createVariable('O1D_SH_loss_prior', 'f4', ('time'))
        O1D_NH_loss_prior = nc.createVariable('O1D_NH_loss_prior', 'f4', ('time'))
        O1D_strat_loss_prior = nc.createVariable('O1D_strat_loss_prior', 'f4', ('time'))

        OH_tot_loss_post = nc.createVariable('OH_tot_loss_post', 'f4', ('time'))
        OH_SH_loss_post = nc.createVariable('OH_SH_loss_post', 'f4', ('time'))
        OH_NH_loss_post = nc.createVariable('OH_NH_loss_post', 'f4', ('time'))
        OH_strat_loss_post = nc.createVariable('OH_strat_loss_post', 'f4', ('time'))
        Cl_tot_loss_post = nc.createVariable('Cl_tot_loss_post', 'f4', ('time'))
        Cl_SH_loss_post = nc.createVariable('Cl_SH_loss_post', 'f4', ('time'))
        Cl_NH_loss_post = nc.createVariable('Cl_NH_loss_post', 'f4', ('time'))
        Cl_strat_loss_post = nc.createVariable('Cl_strat_loss_post', 'f4', ('time'))
        O1D_tot_loss_post = nc.createVariable('O1D_tot_loss_post', 'f4', ('time'))
        O1D_SH_loss_post = nc.createVariable('O1D_SH_loss_post', 'f4', ('time'))
        O1D_NH_loss_post = nc.createVariable('O1D_NH_loss_post', 'f4', ('time'))
        O1D_strat_loss_post = nc.createVariable('O1D_strat_loss_post', 'f4', ('time'))

        # --- set data
        OH_tot_loss_prior[:] = apr['OH_tot_loss_prior'].values
        OH_SH_loss_prior[:] = apr['OH_SH_loss_prior'].values
        OH_NH_loss_prior[:] = apr['OH_NH_loss_prior'].values
        OH_strat_loss_prior[:] = apr['OH_strat_loss_prior'].values
        Cl_tot_loss_prior[:] = apr['Cl_tot_loss_prior'].values
        Cl_SH_loss_prior[:] = apr['Cl_SH_loss_prior'].values
        Cl_NH_loss_prior[:] = apr['Cl_NH_loss_prior'].values
        Cl_strat_loss_prior[:] = apr['Cl_strat_loss_prior'].values
        O1D_tot_loss_prior[:] = apr['O1D_tot_loss_prior'].values
        O1D_SH_loss_prior[:] = apr['O1D_SH_loss_prior'].values
        O1D_NH_loss_prior[:] = apr['O1D_NH_loss_prior'].values
        O1D_strat_loss_prior[:] = apr['OH_strat_loss_prior'].values

        OH_tot_loss_post[:] = pst['OH_tot_loss_post']
        OH_SH_loss_post[:] = pst['OH_SH_loss_post']
        OH_NH_loss_post[:] = pst['OH_NH_loss_post']
        OH_strat_loss_post[:] = pst['OH_strat_loss_post']
        Cl_tot_loss_post[:] = pst['Cl_tot_loss_post']
        Cl_SH_loss_post[:] = pst['Cl_SH_loss_post']
        Cl_NH_loss_post[:] = pst['Cl_NH_loss_post']
        Cl_strat_loss_post[:] = pst['Cl_strat_loss_post']
        O1D_tot_loss_post[:] = pst['O1D_tot_loss_post']
        O1D_SH_loss_post[:] = pst['O1D_SH_loss_post']
        O1D_NH_loss_post[:] = pst['O1D_NH_loss_post']
        O1D_strat_loss_post[:] = pst['OH_strat_loss_post']

        # --- date-time
        dd = []
        for ym in np.arange(self.years[0], self.years[1] + 1, 1):
            for mm in np.arange(1, 13, 1):
                dd.append(datetime(int(ym), int(mm), int(self.mdays[mm - 1])))
        dt = date2num(dd, units=time.units, calendar=time.calendar)
        time[:] = dt[:]

        # --- units
        OH_tot_loss_prior.units = 'g-CH4/m2/month'
        OH_SH_loss_prior.units = 'g-CH4/m2/month'
        OH_NH_loss_prior.units = 'g-CH4/m2/month'
        OH_strat_loss_prior.units = 'g-CH4/m2/month'
        Cl_tot_loss_prior.units = 'g-CH4/m2/month'
        Cl_SH_loss_prior.units = 'g-CH4/m2/month'
        Cl_NH_loss_prior.units = 'g-CH4/m2/month'
        Cl_strat_loss_prior.units = 'g-CH4/m2/month'
        O1D_tot_loss_prior.units = 'g-CH4/m2/month'
        O1D_SH_loss_prior.units = 'g-CH4/m2/month'
        O1D_NH_loss_prior.units = 'g-CH4/m2/month'
        O1D_strat_loss_prior.units = 'g-CH4/m2/month'

        OH_tot_loss_post.units = 'g-CH4/m2/month'
        OH_SH_loss_post.units = 'g-CH4/m2/month'
        OH_NH_loss_post.units = 'g-CH4/m2/month'
        OH_strat_loss_post.units = 'g-CH4/m2/month'
        Cl_tot_loss_post.units = 'g-CH4/m2/month'
        Cl_SH_loss_post.units = 'g-CH4/m2/month'
        Cl_NH_loss_post.units = 'g-CH4/m2/month'
        Cl_strat_loss_post.units = 'g-CH4/m2/month'
        O1D_tot_loss_post.units = 'g-CH4/m2/month'
        O1D_SH_loss_post.units = 'g-CH4/m2/month'
        O1D_NH_loss_post.units = 'g-CH4/m2/month'
        O1D_strat_loss_post.units = 'g-CH4/m2/month'


        print('*** SUCCESS writing for ', invc)
        nc.close()
