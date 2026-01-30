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


class ObsModNc(_set_case.SetCase):
    def __init__(self):
        super().__init__()

        print('Read Input Fluxes from dir: ')
        print('\t\t', self.flx_inp_dir)

        # ---
        self.get_obs_model_nc()

    # --- get_obs_model_nc
    def get_obs_model_nc(self):
        # - get_obs_sites
        def get_obs_sites():
            df = pd.read_csv(self.sites_f, sep='\\s+', skiprows=2, header=None, usecols=range(6))

            # Assign column names manually (based on your data structure)
            df.columns = ['idx', 'sno', 'Sitename',
                          'Lat_o', 'Lon_o', 'Alt_o']
            return df

        # --- get_obs_ts
        def get_obspack_ts(f_name_):
            df_r = pd.read_csv(self.obspack_dir + f_name_ + '.txt', sep='\\s+', comment='#')
            df = df_r[['year', 'month', 'day', 'hour', 'minute', 'value', 'latitude', 'longitude', 'altitude']].copy()
            df.rename(columns={'value': 'o_ch4'}, inplace=True)
            df.loc[:, 'o_ch4'] = df['o_ch4'] * 1e9
            df['index'] = pd.to_datetime(df[['year', 'month', 'day', 'hour', 'minute']])
            df.set_index('index', inplace=True)
            df_om = df.resample('ME').mean()
            return df_om[['o_ch4']]

        # --- get_opost_ts
        def get_opost_ts(f_name_):
            found = False
            for filename in os.listdir(self.opost_dir):
                if f_name_[:7].lower() in filename:
                    found = True
                    print(f"\tFound file: {filename}")
                    df = pd.read_csv(os.path.join(self.opost_dir, filename), sep='\t')
                    df.columns = ['year', 'month', 'day', 'hour', 'minute', 'lat', 'lon', 'alt', 'o_ch4']
                    df.loc[df['o_ch4'] <= -0.0, 'o_ch4'] = np.nan
                    df['index'] = pd.to_datetime(df[['year', 'month', 'day', 'hour', 'minute']])
                    df.set_index('index', inplace=True)
                    df_om = df.resample('ME').mean()
                    df_om.loc[df_om['o_ch4'] <= 1400.0, 'o_ch4'] = np.nan
                    print(f"\t\t Min/Max", df_om['o_ch4'].min().round(2), df_om['o_ch4'].max().round(2))
                    return df_om[['o_ch4']]

            if not found:
                print(f"\tNo file found for pattern: {f_name_[:7]}")
                exit()
                # return pd.DataFrame(columns=['o_ch4'])

        # --- get_model_ts
        def get_model_ts(var_, ds_z_, lat_, lon_, alt_, syr_, inq_):
            # get level = level_val
            z_profile = ds_z_['Z'].isel(time=0).sel(y=lat_, x=lon_, method='nearest')
            level_idx = np.abs(z_profile - alt_).argmin().item()
            level_val = z_profile['level'].values[level_idx]
            # print(alt_, level_val)

            df_a = var_.sel(y=lat_, x=lon_, level=level_val, method='nearest').to_dataframe()
            df_a.reset_index(inplace=True)
            fixed_date = syr_ + '-01-31'
            fixed_date = pd.to_datetime(fixed_date)
            # print(df_a['time'])
            df_a['index'] = pd.to_datetime(fixed_date) + pd.to_timedelta(df_a['time'], unit='D')
            # df_a['index'] = df_a['time'].astype(int).apply(lambda m: fixed_date + pd.DateOffset(months=m))
            df_a.set_index('index', inplace=True)
            df_a.drop(columns=['x', 'y', 'time', 'level'], inplace=True)
            df_a.rename(columns={inq_: 'm_ch4'}, inplace=True)
            return df_a

        # --- get join_df
        def join_df(df_o_, df_m_):
            # print(len(df_m_))
            # print(len(df_o_))
            df_j = df_o_.join(df_m_[['m_ch4']], how='right')
            # df_j.dropna(inplace=True)
            return df_j

        # --- write_obsrvCH4_nc
        def write_obsrvCH4_nc():

            filename = self.inv_dir + 'obsrvCH4_test.nc'
            station_count = len(nc_sites)
            tlen = len(nc_time)
            data_type = 2
            namelen = 45

            # Create a new NetCDF file
            ncfile = Dataset(filename, mode='w', format='NETCDF4')

            # Define dimensions
            ncfile.createDimension('station_count', station_count)
            ncfile.createDimension('tlen', tlen)
            ncfile.createDimension('data_type', data_type)
            ncfile.createDimension('namelen', namelen)

            # Define variables
            site_name = ncfile.createVariable('site_name', 'S1', ('station_count', 'namelen'))
            station_num = ncfile.createVariable('station_count', 'i2', ('station_count',))
            station_num.long_name = "station number"

            tlen_var = ncfile.createVariable('tlen', 'f4', ('tlen',))
            tracer_data = ncfile.createVariable('tracer_data', 'f4', ('tlen', 'station_count', 'data_type'))
            tracer_data.long_name = "monthly mean tracer DMFs"
            tracer_data.units = "ppm"

            # Fill with data
            # site_name[:] = stringtochar(np.array(nc_sites, f'S{namelen}'))
            site_name[:] = stringtochar(np.array([s.ljust(namelen) for s in nc_sites], f'S{namelen}'))
            station_num[:] = nc_idx
            tlen_var[:] = nc_time
            tracer_data[:, :, :] = nc_tracer_data

            ncfile.close()
            print(f"✅ NetCDF file '{filename}' successfully created.")

        # =======================================
        # - obs
        df_obs = get_obs_sites()
        print(df_obs.head())

        # - model
        speriod = str(self.yr_s) + '_' + str(self.yr_e)
        syr = str(self.yr_s)
        path = self.mod_dir + '/' + 'ch4_c1' + '_' + speriod + '.nc'
        ds = xr.open_dataset(path, decode_times=False)
        path = self.mod_dir + '/' + 'z' + '_' + speriod + '.nc'
        ds_z = xr.open_dataset(path, decode_times=False)
        inq = self.invq[0]
        var = ds[inq]

        # - nc actual data
        nc_sites = df_obs.Sitename
        # print(nc_sites)
        nc_idx = df_obs.idx
        # nc_time = (np.arange(self.yr_s, self.yr_e + 1, 1 / 12) + 1 / (2 * 12)).round(3)
        nc_time = (np.arange(self.yr_s-2, self.yr_e-2 + 1, 1 / 12) + 1 / (2 * 12)).round(3)
        # print(nc_time)

        # - preallocate nc array
        nc_tracer_data = np.zeros((len(nc_time), len(df_obs), 2), dtype=np.float32)
        # print(nc_tracer_data.shape)

        # - loop sites
        for index, row in df_obs.iloc[:].iterrows():
            idx, sno, f_name, lat, lon, alt = (row.iloc[0], row.iloc[1], row.iloc[2],
                                               row.iloc[3], row.iloc[4], row.iloc[5])
            print(idx, sno, f_name)
            if sno == 9999:
                # df_o = get_obspack_ts(f_name)
                df_o = get_opost_ts(f_name)
                # print(df_o.head())
                df_m = get_model_ts(var, ds_z, lat, lon, alt, syr, inq)
                # print(df_m.head())
                df_j = join_df(df_o, df_m)
                # print(len(df_j))
                # print(df_j.head())
                # print(df_j.tail())
                nc_tracer_data[:, idx - 1, 0] = df_j['m_ch4'].values/50. #- df_j['o_ch4'].values
                nc_tracer_data[:, idx - 1, 1] = df_j['m_ch4'].values/100. #- df_j['o_ch4'].values

        # --- write_obsrvCH4_nc
        write_obsrvCH4_nc()
