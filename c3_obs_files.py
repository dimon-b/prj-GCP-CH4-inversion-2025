# -*- coding: utf-8 -*-
"""
Created/Corrected:    25/09/09 09:09
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import os
import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset, stringtochar

import _set_case


class ObsFilesTxt(_set_case.SetCase):
    def __init__(self):
        super().__init__()

        print('Read Input Fluxes from dir: ')
        print('\t\t', self.inp_dir)

        # ---
        self.get_obs_files_txt()

    # --- get_obs_model_nc
    def get_obs_files_txt(self):
        # - get_obs_sites
        def get_obs_sites():
            df = pd.read_csv(self.sites_f, sep='\\s+', skiprows=2, header=None, usecols=range(6))
            df.columns = ['idx', 'sno', 'Sitename', 'Lat_o', 'Lon_o', 'Alt_o']
            return df

        # --- check_obspack
        def check_obspack(dir, f_name_):
            # - read file
            def get_obspack_ts():
                df_r = pd.read_csv(dir + filename, sep='\\s+', comment='#')
                df = df_r[['year', 'month', 'day', 'hour', 'value']].copy()
                df.rename(columns={'value': 'o_ch4'}, inplace=True)
                df.loc[:, 'o_ch4'] = df['o_ch4'] * 1e9
                return df

            found = False
            df_list = []

            # - check exact filename first
            for root, dirs, files in os.walk(dir):
                for filename in files:
                    if f_name_.lower() in filename.lower():
                        found = True
                        print(f"\tFound file: {filename}")
                        df = get_obspack_ts()
                        return found, df, None

            # - if no exact match, check first 7 chars
            if not found:
                for root, dirs, files in os.walk(dir):
                    for filename in files:
                        if f_name_[:11].lower() in filename.lower():
                            found = True
                            print(f"\tFound partial file: {filename}")
                            df = get_obspack_ts()
                            return found, df, filename

            if not found:
                print(f"\tNo file found for pattern: {f_name_}")
                return False, None

        # =======================================
        # - obs
        df_obs = get_obs_sites()
        print(df_obs.head())

        # - loop sites
        for index, row in df_obs.iloc[:94].iterrows():
            idx, sno, f_name, lat, lon, alt = (row.iloc[0], row.iloc[1], row.iloc[2],
                                               row.iloc[3], row.iloc[4], row.iloc[5])
            # if sno == 9999:
            if sno >-1:
                print(idx, sno, f_name)
                found = False

                # check_obspack
                found, df, f_name_n = check_obspack(self.obspack_dir, f_name)
                if found:
                    if f_name_n == None:
                        print(f"\tGet file date from Obspack: {f_name}")
                        df.to_csv(self.obsout_dir + f_name + '.txt', sep='\t', index=False, header=False)
                    else:
                        print(f"\t\tGet partial file date from Obspack: {f_name}")
                        df.to_csv(self.obsout_dir + f_name + '.txt', sep='\t', index=False, header=False)

                if not found:
                    print(f"\tNo file found for pattern: {f_name}")
                    exit()
