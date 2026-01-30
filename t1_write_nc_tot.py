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
import xarray as xr

import _set_case


class WriteNcTot(_set_case.SetCase):

    def __init__(self):
        super().__init__()
        self.invcases = ['inv1']
        self.flxcases = ['gcp2021_v2_soil0_inv1']

        self.unp = 'p30'
        self.unx = 'ctl'
        self.wrt_flux_tot()

    # --- write flux categories
    def wrt_flux_tot(self):
        # - get_post_flx_grd
        def get_post_flx_grd(unp, unx, invc):
            # - read bin
            def read_flx_bin(bfile):
                size_2d = self.d1_nlon * self.d1_nlat * self.nmonth
                with open(bfile, 'rb') as f:
                    fl2d = np.fromfile(f, dtype='<f4', count=size_2d * self.nyear_nc
                                       ).reshape(self.nyear_nc, self.nmonth, self.d1_nlat, self.d1_nlon)

                print(f'\t\t\t post flux shape from grd file: {fl2d.shape}')
                return fl2d

            # - flux shape correction [yr, mn, :, :] -> [time, :, :]
            def cor_shape(flx):
                flx1 = np.zeros([self.nyear_nc * self.nmonth, self.d1_nlat, self.d1_nlon])
                for yr in np.arange(self.nyear_nc):
                    for mn in np.arange(self.nmonth):
                        ss = yr * self.nmonth + mn
                        flx1[ss, :, :] = flx[yr, mn, :, :]
                print(f'\t\t\t post flux shape from after correction: {flx1.shape}, '
                      f'time length: {self.nyear_nc * self.nmonth}, years: {self.nyear_nc}, months: {self.nmonth}')
                return flx1

            # ===========================================================
            # - post
            i_file = self.inv_pst_dir + unp + '/' + self.icase + '_' + unx + '_' + invc + '.grd'
            print('\t\tPosterior flux reading: ', invc, i_file)
            flx = read_flx_bin(i_file)
            flx1 = cor_shape(flx)
            return flx1

        # - get_prior_full_flx_nc
        def get_prior_full_flx_nc():
            # --- get_inp_nc
            def get_prior_full_nc():
                path = self.flx_apr_dir + self.flx_apr_nc_full
                ds = xr.open_dataset(path, decode_times=False)
                return ds

            ds = get_prior_full_nc()
            return ds

        # - get_joined_apr2post
        def get_joined_apr2post(apr, pst):
            # - soil from apr
            apr_soil = apr["flux_ch4_soils"]

            # - extend pst
            pst_ext = np.concatenate([pst[:12], pst, pst[-12:], ], axis=0, )

            # - pst_ext into a DataArray
            pst_da = xr.DataArray(pst_ext, dims=apr_soil.dims, coords=apr_soil.coords, name="pst")
            assert pst_da.shape == apr_soil.shape
            assert pst_da.time.equals(apr_soil.time)

            pst_flux = xr.Dataset({
                "fch4_soils": apr["flux_ch4_soils"],
                "fch4_prior": apr["flux_ch4_prior"],
                "fch4_prior_soil0": apr["flux_ch4_prior_soil0"],
                "fch4_post": pst_da,
                "fch4_post_soil0": pst_da - apr["flux_ch4_soils"],
                "fch4_corr": pst_da - apr["flux_ch4_prior"],
            })
            if not np.issubdtype(pst_flux.time.dtype, np.datetime64):
                pst_flux['time'] = xr.decode_cf(pst_flux[['time']]).time

            return pst_flux

        # - check_total
        def check_total(ds):

            # - Grid-cell area [m2]
            dlat = np.deg2rad(ds.lat[1] - ds.lat[0])
            dlon = np.deg2rad(ds.lon[1] - ds.lon[0])
            lat_rad = np.deg2rad(ds.lat)
            cell_area_2d = (self.R ** 2 * dlat * dlon * np.cos(lat_rad))
            cell_area = cell_area_2d.broadcast_like(ds["fch4_prior"])

            # - sanity check: total Earth surface area
            earth_area = cell_area.isel(time=0).sum(dim=("lat", "lon"))
            print("")
            print(f"Check total Earth area = {earth_area.item():.2e} m2")
            print("Expected ≈ 5.10e14 m2")
            print("")

            # - Seconds per month (time-aware, xarray-safe)
            time = pd.to_datetime(ds.time.values)
            seconds_per_month = xr.DataArray(time.days_in_month * 24 * 3600,
                                             coords={"time": ds.time},
                                             dims=("time",),
                                             name="seconds_per_month",
                                             )

            # - Flux variables to integrate
            FLUX_VARS = [v for v in ds.data_vars if v.startswith("fch4_")]

            # - Integration
            records = []
            for var in FLUX_VARS:
                # kg m-2 s-1 -> kg per month per grid cell
                monthly_flux = ds[var] * cell_area * seconds_per_month
                # global monthly sum [kg]
                global_monthly = monthly_flux.sum(dim=("lat", "lon"))
                # annual totals [Tg yr-1]
                global_annual = (global_monthly.groupby("time.year").sum(dim="time") / 1e9)
                for year, value in zip(global_annual.year.values, global_annual.values):
                    records.append(
                        {"year": int(year), "variable": var.replace("fch4_", ""), "Tg_CH4_yr": float(value), })

            # - Output table
            df = pd.DataFrame(records)
            table = df.pivot(index="year", columns="variable", values="Tg_CH4_yr", )

            # - Put totals at the end
            total_cols = [c for c in table.columns if c.lower().startswith("total")]
            other_cols = [c for c in table.columns if c not in total_cols]
            table = table[other_cols + total_cols]

            print(table.round(4))

        # - write_1nc
        def write_1nc(invc, apr, pst):

            fout = self.inv_ncd_dir + '/MIROC4-ACTM_totflux_' + invc + 'SURF.nc'
            nc = netCDF4.Dataset(fout, 'w', format='NETCDF4')
            nc.description = 'Net prior and posteor CH4 emissions resulted from the surface based inversion for GCP-CH4, 2025. The results are produced at JAMSTEC, Japan. For details please contact at prabir@jamstec.go.jp and d.belikov@chiba-u.jp.'
            nc.Institution = 'Japan Agency for Marine-Earth Science and Technology (JAMSTEC)'
            nc.Contact = 'd.belikov@chiba-u.jp and prabir@jamstec.go.jp'
            nc.Disclaimer = 'This data is created for GCP-CH4 2025.'
            nc.History = 'Created on January 2026 by Dmitry Belikov'
            nc.Additional_Info = 'MIROC4.0-based Atmospheric Chemistry-Transport Model (MIROC4-ACTM) is used for forward simulations (Prabir et al., 2018) and Time-independent Bayesian inversion is used for optimising flux (Patra et al., 2016)'

            # --- set dimensions
            time_dim = nc.createDimension('time', None)
            lat_dim = nc.createDimension('lat', self.d1_nlat)
            lon_dim = nc.createDimension('lon', self.d1_nlon)
            # carea = nc.createDimension('carea', self.d1_nlat)

            # --- set variables
            time = nc.createVariable('time', 'f4', ('time',))
            time.units = 'hours since 1970-01-01 00:00:00'
            time.calendar = 'gregorian'
            lat = nc.createVariable('lat', 'f4', ('lat',))
            lat.units = 'degrees north from -90'
            lat.long_name = 'latitudes'
            lon = nc.createVariable('lon', 'f4', ('lon',))
            lon.units = 'degrees east from -180'
            lon.long_name = 'longitudes'
            cellarea = nc.createVariable('carea', 'f4', ('lat',))
            prior = nc.createVariable('fch4_prior', 'f4', ('time', 'lat', 'lon'))
            post = nc.createVariable('fch4_post', 'f4', ('time', 'lat', 'lon'))

            # --- set data
            lat[:] = self.d1_lats
            lon[:] = self.d1_lons
            cellarea[:] = self.garia_d1
            prior[:, :, :] = apr[:, :, :]
            post[:, :, :] = pst[:, :, :]

            # --- check G-total
            nms = ['apr', 'pst']
            vrs = [apr, pst]

            # --- date-time
            dd = []
            for ym in np.arange(self.years_nc[0], self.years_nc[1] + 1, 1):
                for mm in np.arange(1, 13, 1):
                    dd.append(datetime(int(ym), int(mm), int(self.mdays[mm - 1])))
            dt = date2num(dd, units=time.units, calendar=time.calendar)
            time[:] = dt[:]

            # --- units
            cellarea.units = 'm2'
            prior.units = 'g-CH4/m2/month'
            post.units = 'g-CH4/m2/month'
            prior.valid_range = np.array((np.min(prior), np.max(prior)))
            post.valid_range = np.array((np.min(post), np.max(post)))

            print('*** SUCCESS writing for ', invc)
            nc.close()

        # ===========================================================
        print(f'\tRun wrt_flux_cat')
        print(f'\tApr flux from: {self.inv_apr_dir}')
        print(f'\tPst flux from: {self.inv_pst_dir}')
        print(f'\tOutput *.nc file length: {self.nyear_nc} years for {self.years_nc} \n')

        for j1, invc in enumerate(self.invcases):
            flx_prior = get_prior_full_flx_nc()
            flx_post = get_post_flx_grd(self.unp, self.unx, invc)
            flx_join = get_joined_apr2post(flx_prior, flx_post)
            check_total(flx_join)
            # todo
            write_1nc()
