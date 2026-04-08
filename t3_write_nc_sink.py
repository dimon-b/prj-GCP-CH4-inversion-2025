# -*- coding: utf-8 -*-
"""
Created:    22/01/26 19:17
Project:    Write total fCH4 *.nc for GCP2021 as originally created by Naveen Negi
@author:    Dmitry Belikov
"""
from datetime import datetime
import netCDF4
import pandas as pd
from netCDF4 import date2num
from pathlib import Path

# import b_sub_flux
import _set_case


class WriteNcLoss(_set_case.SetCase):

    def __init__(self):
        super().__init__()
        self.invcases = ['inv1']  # , 'inv7', 'vCao']
        self.dirs = ['../results2025/CYC/', '../results2025/INCA/']

        # --- run FORTRAN to get loss txt
        # ./c_gcpv3_f/run t3_ch4_sink.f90

        # --- run here: txt -> nc
        for ddir in self.dirs:
            self.wrt_loss_nc(ddir)

    # ---
    def wrt_loss_nc(self, dir_):

        # ---
        def read_loss_txt():
            def get_1file(ifile_: str) -> pd.DataFrame:
                df = pd.read_csv(ifile_, sep=r'\s+', skiprows=0)
                return df

            ifile = dir_ + 'losscorr/' + 'sink_prior_' + invc + '.txt'
            print(f'\t\tRead ifile: {ifile}')
            df_apr = get_1file(ifile)
            # print(df_apr.head())
            df_apr.columns = ['year', 'month',
                              'OH_tot_loss_prior', 'Cl_tot_loss_prior', 'O1D_tot_loss_prior',
                              'OH_SH_loss_prior', 'Cl_SH_loss_prior', 'O1D_SH_loss_prior',
                              'OH_NH_loss_prior', 'Cl_NH_loss_prior', 'O1D_NH_loss_prior',
                              'OH_strat_loss_prior', 'Cl_strat_loss_prior', 'O1D_strat_loss_prior',
                              'OH_tropo_loss_prior', 'Cl_tropo_loss_prior', 'O1D_tropo_loss_prior',
                              ]
            # print(df_apr)
            df_apr_1 = df_apr.copy()
            df_apr_1['OH_tot_loss_ck'] = df_apr_1['OH_SH_loss_prior'] + df_apr_1['OH_NH_loss_prior']
            df_apr_1[['OH_tot_loss_prior', 'OH_SH_loss_prior', 'OH_NH_loss_prior', 'OH_tot_loss_ck']].plot()
            # print(df_apr_1[['OH_tot_loss_prior', 'OH_SH_loss_prior', 'OH_NH_loss_prior', 'OH_tot_loss_ck']])
            # plt.show()
            df_apr_1 = df_apr.copy()
            df_apr_1['Cl_tot_loss_ck'] = df_apr_1['Cl_SH_loss_prior'] + df_apr_1['Cl_NH_loss_prior']
            df_apr_1[['Cl_tot_loss_prior', 'Cl_SH_loss_prior', 'Cl_NH_loss_prior', 'Cl_tot_loss_ck']].plot()
            # print(df_apr_1[['Cl_tot_loss_prior', 'Cl_SH_loss_prior', 'Cl_NH_loss_prior', 'Cl_tot_loss_ck']])
            # plt.show()
            df_apr_1 = df_apr.copy()
            df_apr_1['O1D_tot_loss_ck'] = df_apr_1['O1D_SH_loss_prior'] + df_apr_1['O1D_NH_loss_prior']
            df_apr_1[['O1D_tot_loss_prior', 'O1D_SH_loss_prior', 'O1D_NH_loss_prior', 'O1D_tot_loss_ck']].plot()
            # print(df_apr_1[['O1D_tot_loss_prior', 'O1D_SH_loss_prior', 'O1D_NH_loss_prior', 'O1D_tot_loss_ck']])
            # plt.show()

            ifile = dir_ + 'losscorr/' + 'sink_post_' + invc + '.txt'
            print(f'\t\tRead ifile: {ifile}')
            df_pst = get_1file(ifile)
            df_pst.columns = ['year', 'month',
                              'OH_tot_loss_post', 'Cl_tot_loss_post', 'O1D_tot_loss_post',
                              'OH_SH_loss_post', 'Cl_SH_loss_post', 'O1D_SH_loss_post',
                              'OH_NH_loss_post', 'Cl_NH_loss_post', 'O1D_NH_loss_post',
                              'OH_strat_loss_post', 'Cl_strat_loss_post', 'O1D_strat_loss_post',
                              'OH_tropo_loss_post', 'Cl_tropo_loss_post', 'O1D_tropo_loss_post',
                              ]

            def cat_time(df_):
                # df_['date'] = pd.to_datetime(dict(year=df_.year, month=df_.month, day=1))
                df_['date'] = pd.to_datetime(df_[['year', 'month']].assign(day=1))
                df_ = df_.set_index('date')
                df_ = df_.loc['2000-01':'2024-12']
                return df_

            return cat_time(df_apr), cat_time(df_pst)

        # ---
        def write_1nc_sink(apr, pst):
            file_out = ''
            if 'CYC' in dir_:
                file_out = Path('../results2025/CYC/nc_out/') / 'MIROC4-ACTM_sink_GMB_SURF_OH_Transcom.nc'
            elif 'INCA' in dir_:
                file_out = Path('../results2025/INCA/nc_out/') / 'MIROC4-ACTM_sink_GMB_SURF_OH_INCA.nc'

            nc = netCDF4.Dataset(file_out, 'w', format='NETCDF4')
            nc.description = 'Net prior and posteor CH4 emissions resulted from the surface based inversion for GCP-CH4, 2025. The results are produced at JAMSTEC, Japan. For details please contact at prabir@jamstec.go.jp and d.belikov@chiba-u.jp.'
            nc.Institution = 'Japan Agency for Marine-Earth Science and Technology (JAMSTEC) and Chiba University'
            nc.Contact = 'd.belikov@chiba-u.jp and prabir@jamstec.go.jp'
            nc.Disclaimer = 'This data is created for GCP-CH4 2025.'
            nc.History = 'April 2026'
            nc.Additional_Info = 'MIROC4.0-based Atmospheric Chemistry-Transport Model (MIROC4-ACTM) is used for forward simulations (Prabir et al., 2018) and Time-independent Bayesian inversion is used for optimising flux (Patra et al., 2016)'

            # --- set dimensions
            nc.createDimension('time', None)

            # --- set variables
            time = nc.createVariable('time', 'f4', ('time',))
            time.units = 'hours since 1970-01-01 00:00:00'
            time.calendar = 'Gregorian'
            OH_tot_loss_prior = nc.createVariable('OH_tot_loss_prior', 'f4', 'time')
            OH_SH_loss_prior = nc.createVariable('OH_SH_loss_prior', 'f4', 'time')
            OH_NH_loss_prior = nc.createVariable('OH_NH_loss_prior', 'f4', 'time')
            OH_strat_loss_prior = nc.createVariable('OH_strat_loss_prior', 'f4', 'time')
            Cl_tot_loss_prior = nc.createVariable('Cl_tot_loss_prior', 'f4', 'time')
            Cl_SH_loss_prior = nc.createVariable('Cl_SH_loss_prior', 'f4', 'time')
            Cl_NH_loss_prior = nc.createVariable('Cl_NH_loss_prior', 'f4', 'time')
            Cl_strat_loss_prior = nc.createVariable('Cl_strat_loss_prior', 'f4', 'time')
            O1D_tot_loss_prior = nc.createVariable('O1D_tot_loss_prior', 'f4', 'time')
            O1D_SH_loss_prior = nc.createVariable('O1D_SH_loss_prior', 'f4', 'time')
            O1D_NH_loss_prior = nc.createVariable('O1D_NH_loss_prior', 'f4', 'time')
            O1D_strat_loss_prior = nc.createVariable('O1D_strat_loss_prior', 'f4', 'time')

            OH_tot_loss_post = nc.createVariable('OH_tot_loss_post', 'f4', 'time')
            OH_SH_loss_post = nc.createVariable('OH_SH_loss_post', 'f4', 'time')
            OH_NH_loss_post = nc.createVariable('OH_NH_loss_post', 'f4', 'time')
            OH_strat_loss_post = nc.createVariable('OH_strat_loss_post', 'f4', 'time')
            Cl_tot_loss_post = nc.createVariable('Cl_tot_loss_post', 'f4', 'time')
            Cl_SH_loss_post = nc.createVariable('Cl_SH_loss_post', 'f4', 'time')
            Cl_NH_loss_post = nc.createVariable('Cl_NH_loss_post', 'f4', 'time')
            Cl_strat_loss_post = nc.createVariable('Cl_strat_loss_post', 'f4', 'time')
            O1D_tot_loss_post = nc.createVariable('O1D_tot_loss_post', 'f4', 'time')
            O1D_SH_loss_post = nc.createVariable('O1D_SH_loss_post', 'f4', 'time')
            O1D_NH_loss_post = nc.createVariable('O1D_NH_loss_post', 'f4', 'time')
            O1D_strat_loss_post = nc.createVariable('O1D_strat_loss_post', 'f4', 'time')

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
            dd = [datetime(year, month, self.mdays[month - 1])
                  for year in range(self.years_nc[0], self.years_nc[1] + 1)
                  for month in range(1, 13)]
            time[:] = date2num(dd, units=time.units, calendar=time.calendar)

            # --- units
            loss_vars = [
                # PRIOR
                OH_tot_loss_prior, OH_SH_loss_prior, OH_NH_loss_prior, OH_strat_loss_prior,
                Cl_tot_loss_prior, Cl_SH_loss_prior, Cl_NH_loss_prior, Cl_strat_loss_prior,
                O1D_tot_loss_prior, O1D_SH_loss_prior, O1D_NH_loss_prior, O1D_strat_loss_prior,

                # POST
                OH_tot_loss_post, OH_SH_loss_post, OH_NH_loss_post, OH_strat_loss_post,
                Cl_tot_loss_post, Cl_SH_loss_post, Cl_NH_loss_post, Cl_strat_loss_post,
                O1D_tot_loss_post, O1D_SH_loss_post, O1D_NH_loss_post, O1D_strat_loss_post]

            for var in loss_vars:
                var.units = 'Tg-CH4/month'

            print('\t\tSUCCESS writing sink to >>> ', file_out)
            nc.close()

        # ===========================================================
        print(f'\n\t\tRun wrt_loss_nc for: {dir_}')
        for j1, invc in enumerate(self.invcases[:]):
            a_los, p_los = read_loss_txt()
            write_1nc_sink(a_los, p_los)

