# -*- coding: utf-8 -*-
"""
Created:    22/01/26 19:17
Project:    Write total fCH4 *.nc for GCP2021 as originally created by Naveen Negi
@author:    Dmitry Belikov
"""
import os
from datetime import datetime
import netCDF4
import pandas as pd
from netCDF4 import date2num
from pathlib import Path

# import b_sub_flux
import _set_case


class WriteNcComp(_set_case.SetCase):

    def __init__(self):
        super().__init__()
        self.invcases = ['inv1']
        dirs = [['../trout_MIROC/CH4_t42l67_CYC_H_251014/', '../trout_MIROC/CH4_t42l67_CYC_H_post_260323/'],
                ['../trout_MIROC/CH4_t42l67_INCA_H_260306/', '../trout_MIROC/CH4_t42l67_INCA_H_post_260324/']
                ]

        # --- run
        # /S/data01/G5070/y0715/ccode/prj_apack/c_apack/b3_get_Model_hr_GCP.py

        # --- run here: txt -> nc
        for ddir in dirs:
            # --- get list of csv file
            apr_csv_files = self.get_csv_files(ddir[0])
            print('')
            print("List of apr csv files :")
            print(",\n".join(f'\t\t{f}' for f in apr_csv_files))
            pst_csv_files = self.get_csv_files(ddir[1])
            print("List of pst csv files :")
            print(",\n".join(f'\t\t{f}' for f in pst_csv_files))

            # ---
            for apr_csv, pst_csv in zip(apr_csv_files, pst_csv_files):
                if 'station' in apr_csv:
                    filename = os.path.basename(str(apr_csv))
                    cmp_id = filename.split("_")[0]
                    self.wrt_comp_nc(apr_csv, pst_csv, cmp_id)

    # ---
    def wrt_comp_nc(self, apr_csv_, pst_csv_, id_):

        # ---
        def read_comp_txt():
            def get_1file(ifile_: str) -> pd.DataFrame:
                df = pd.read_csv(ifile_)#, sep=r'\s+', skiprows=0)
                return df

            ifile = apr_csv_
            print(f'\t\tRead ifile: {ifile}')
            df_apr = get_1file(ifile)
            print(df_apr.head())

            ifile = pst_csv_
            print(f'\t\tRead ifile: {ifile}')
            df_pst = get_1file(ifile)
            print(df_pst.head())
            # exit()

            return df_apr, df_pst

        # ---
        def write_1nc_comp(apr, pst):
            file_out = ''
            if 'CYC' in apr_csv_:
                file_out = Path('../results2025/CYC/nc_out/') / f'MIROC4-ACTM_{id_}_comp_GMB_SURF_OH_Transcom.nc'
            elif 'INCA' in apr_csv_:
                file_out = Path('../results2025/INCA/nc_out/') / f'MIROC4-ACTM_{id_}_comp_GMB_SURF_OH_INCA.nc'

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
            time = nc.createVariable('time', 'f8', ('time',))
            time.units = 'hours since 1970-01-01 00:00:00'
            time.calendar = 'Gregorian'
            comp_prior = nc.createVariable('comp_prior', 'f4', 'time')
            comp_post = nc.createVariable('comp_post', 'f4', 'time')

            # --- set data
            comp_prior[:] = apr['ch4_a'].values
            comp_post[:] = pst['ch4_a'].values

            # --- date-time
            apr['dateTime'] = pd.to_datetime(apr['dateTime'])
            dd = apr['dateTime'].dt.to_pydatetime()
            # print(apr)
            # print(dd)
            time[:] = date2num(dd, units=time.units, calendar=time.calendar)
            # print(time[:])
            # exit()

            # --- units
            loss_vars = [comp_prior, comp_post]

            for var in loss_vars:
                var.units = 'ppb'

            print('\t\t\tSUCCESS writing comp to >>> ', file_out)
            nc.close()

        # ===========================================================
        print(f'\n\t\tRun read_comp_txt for: {apr_csv_, pst_csv_, id_}')
        for j1, invc in enumerate(self.invcases[:]):
            a_conc, p_conc = read_comp_txt()
            write_1nc_comp(a_conc, p_conc)

    # ---
    @staticmethod
    def get_csv_files(dir_):
        return [os.path.join(dir_, f) for f in os.listdir(dir_) if f.lower().endswith('.csv')]
