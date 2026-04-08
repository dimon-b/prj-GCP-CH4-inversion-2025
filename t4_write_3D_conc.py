# -*- coding: utf-8 -*-
"""
Created:    22/01/26 19:17
Project:    Write total fCH4 *.nc for GCP2021 as originally created by Naveen Negi
@author:    Dmitry Belikov
"""
import netCDF4
import numpy as np
import xarray as xr
from pathlib import Path

import _set_case


class Write3dConc(_set_case.SetCase):

    def __init__(self):
        super().__init__()
        self.invcases = ['inv1']  # , 'inv7', 'vCao']
        self.dirs = ['D:/OneDrive - 国立大学法人千葉大学/prj_GCP_v25/trout_MIROC/CH4_t42l67_CYC_M_post_260302/',
                     'D:/OneDrive - 国立大学法人千葉大学/prj_GCP_v25/trout_MIROC/CH4_t42l67_INCA_M_post_260311/']

        # - get nc files for each case (CYC, INCA) and run (apr, post)
        for dir_ in self.dirs:
            self.wrt_3d_conc_nc(dir_)

        # - combine files for (CYC, INCA)
        self.write_1nc_conc(self.dirs[0], 'MIROC4-ACTM_CH4conc_GMB_SURF_OH_Transcom.nc')
        self.write_1nc_conc(self.dirs[1], 'MIROC4-ACTM_CH4conc_GMB_SURF_OH_INCA.nc')

    # ===
    def write_1nc_conc(self, dir_, f_out):
        file_inp = Path(dir_) / "MIROC4-ACTM_CH4conc.nc"
        ds_ = xr.open_dataset(file_inp, engine="h5netcdf")
        ds = ds_.sel(time=slice(f"{self.years_nc[0]}-01-01", f"{self.years_nc[1]}-12-31"))

        file_out = ''
        if 'CYC' in dir_:
            file_out = Path('../results2025/CYC/nc_out/') / f_out
        elif 'INCA' in dir_:
            file_out = Path('../results2025/INCA/nc_out/') / f_out

        nc = netCDF4.Dataset(file_out, 'w', format='NETCDF4')
        nc.description = 'Net prior and posteor CH4 emissions resulted from the surface based inversion for GCP-CH4, 2025. The results are produced at JAMSTEC, Japan. For details please contact at prabir@jamstec.go.jp and d.belikov@chiba-u.jp.'
        nc.Institution = 'Japan Agency for Marine-Earth Science and Technology (JAMSTEC) and Chiba University'
        nc.Contact = 'd.belikov@chiba-u.jp and prabir@jamstec.go.jp'
        nc.Disclaimer = 'This data is created for GCP-CH4 2025.'
        nc.History = 'April 2026'
        nc.Additional_Info = 'MIROC4.0-based Atmospheric Chemistry-Transport Model (MIROC4-ACTM) is used for forward simulations (Prabir et al., 2018) and Time-independent Bayesian inversion is used for optimising flux (Patra et al., 2016)'

        # - Create dimensions
        for dim_name, dim_size in ds.sizes.items():
            nc.createDimension(str(dim_name), int(dim_size))

        # - Create variables
        for var_name, da in ds.variables.items():
            name = str(var_name)
            dims = tuple(str(d) for d in da.dims)

            dtype = da.dtype

            # Handle datetime64 safely
            if np.issubdtype(dtype, np.datetime64):
                dtype = np.dtype("float64")  # netCDF-safe float64 for time

            # Force native endianness to avoid warning
            dtype = np.dtype(dtype).newbyteorder("=")

            # --- create variable ---
            nc_var = nc.createVariable(name, dtype, dims, zlib=True)

            # Copy attributes
            for attr_name in da.attrs:
                setattr(nc_var, attr_name, da.attrs[attr_name])

            # Write data
            if np.issubdtype(da.dtype, np.datetime64):
                # Convert time to numeric (important!)
                time_units = 'days since 1970-01-01 00:00:00'
                calendar = 'standard'
                nc_var.units = time_units
                nc_var.calendar = calendar
                nc_var[:] = netCDF4.date2num(
                    da.values.astype('datetime64[ns]').astype('datetime64[ms]').astype(object),
                    units=time_units,
                    calendar=calendar
                )
            else:
                nc_var[:] = da.values

        print('\t\t SUCCESS writing 3d conc to >>> ', file_out)
        nc.close()

    # ---
    def wrt_3d_conc_nc(self, dir_):
        # ---
        def build_ch4_dataset(ddir_, invc_):
            ddir_path = Path(ddir_)
            ch4_file = ddir_path / f"ch4_{invc_}_1999_2025.nc"
            ps_file = ddir_path / "ps_1999_2025.nc"
            t_file = ddir_path / "t_1999_2025.nc"

            ds_ch4 = xr.open_dataset(ch4_file, engine="h5netcdf")
            ds_ps = xr.open_dataset(ps_file, engine="h5netcdf")
            ds_t = xr.open_dataset(t_file, engine="h5netcdf")

            # print(ds_ch4)
            # print(ds_ps)
            # print(ds_t)

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
            # print(ds)

            fout = Path(ddir_) / "MIROC4-ACTM_CH4conc.nc"
            encoding = {v: {'_FillValue': -999.0} for v in ds.data_vars}
            ds.to_netcdf(fout, encoding=encoding, engine="h5netcdf")

        # ===========================================================
        print(f'\t\t Run wrt_3D_conc_nc for: {dir_}')
        for j1, invc in enumerate(self.invcases[:]):
            build_ch4_dataset(dir_, invc)
