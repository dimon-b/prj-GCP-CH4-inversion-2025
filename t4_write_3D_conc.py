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
import matplotlib.pyplot as plt
import xarray as xr

import _set_case


class Write3dConc(_set_case.SetCase):

    def __init__(self):
        super().__init__()
        self.invcases = ['inv1']  # , 'inv7', 'vCao']
        self.dirs = ['D:/dbase/ACTM/CH4_t42l67_CYC_M_post_260302/', 'D:/dbase/ACTM/CH4_t42l67_INCA_M_post_260311/']

        for dir in self.dirs:
            self.wrt_3D_conc_nc(dir)

    # ---
    def wrt_3D_conc_nc(self, dir):
        print(f'\tRun wrt_3D_conc_nc for: {dir}')

        for j1, invc in enumerate(self.invcases[:]):
            self.build_ch4_dataset(dir, invc)

    # ---
    def build_ch4_dataset(self, dir, invc_):

        ch4_file = dir + 'ch4_' + invc_ + '_1999_2025.nc'
        ps_file = dir + 'ps_1999_2025.nc'
        t_file = dir + 't_1999_2025.nc'

        ds_ch4 = xr.open_dataset(ch4_file)
        ds_ps = xr.open_dataset(ps_file)
        ds_t = xr.open_dataset(t_file)

        print(ds_ch4)
        print(ds_ps)
        print(ds_t)
        # exit()

        # align
        ds_ch4, ds_ps, ds_t = xr.align(ds_ch4, ds_ps, ds_t, join='inner')

        # merge
        ds = xr.Dataset({'ch4_conc': ds_ch4['Q03'],
                         'Ps': ds_ps['PS'] * 1e2,
                         'temperature': ds_t['T']
                         })

        # names & attributes
        rename_map = {'x': 'longitude', 'y': 'latitude', 'plev': 'level'}
        ds = ds.rename({k: v for k, v in rename_map.items() if k in ds_ch4.dims})
        ds['level'].attrs.update({'long_name': 'Sigma-pressure level'})
        ds['ch4_conc'].attrs.update({'long_name': 'Methane concentration'})
        ds['Ps'].attrs.update({'long_name': 'Surface pressure'})
        ds['temperature'].attrs.update({'long_name': 'Air temperature'})

        ds['longitude'].attrs['actual_range'] = [float(ds.longitude.min()), float(ds.longitude.max())]
        ds['latitude'].attrs['actual_range'] = [float(ds.latitude.min()), float(ds.latitude.max())]
        ds['level'].attrs['actual_range'] = [float(ds.level.min()), float(ds.level.max())]
        ds['ch4_conc'].attrs['actual_range'] = [float(ds.ch4_conc.min()), float(ds.ch4_conc.max())]
        ds['Ps'].attrs['actual_range'] = [float(ds.Ps.min()), float(ds.Ps.max())]
        ds['temperature'].attrs['actual_range'] = [float(ds.temperature.min()), float(ds.temperature.max())]
        print(ds)

        fout = self.flx_ncd_dir + '/MIROC4-ACTM_sink_GMB_SURF_OH_Transcom.nc'
        if 'CYC' in dir:
            fout = dir + '/MIROC4-ACTM_CH4conc_OH_Transcom_SURF.nc'
        elif 'INCA' in dir:
            fout = dir + '/MIROC4-ACTM_CH4conc_OH_INCA_SURF.nc'

        encoding = {v: {'_FillValue': -999.0} for v in ds.data_vars}
        ds.to_netcdf(fout, encoding=encoding)
