# -*- coding: utf-8 -*-
"""
Created:    22/01/26 19:17
Project:    Write total fCH4 *.nc for GCP2021 as originally created by Naveen Negi
@author:    Dmitry Belikov
"""
from datetime import datetime
import calendar
import netCDF4
import numpy as np
import pandas as pd
from netCDF4 import date2num

import _set_case


class WriteNcTot(_set_case.SetCase):

    def __init__(self):
        super().__init__()
        self.invcases = ['inv1']
        self.flxcases = ['fch4_gcp2025_soil0_inv1']

        self.unp = 'p30'
        self.unx = 'ux4'
        self.wrt_flux_tot()

    # --- wrt_flux_tot
    def wrt_flux_tot(self):
        # --- write_1nc
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
            Gtotal_tot(nms, vrs, invc)

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

        # --- Gtotal_tot
        def Gtotal_tot(nms, vrs, invc):
            print(f'\tGtotal_tot: ', invc, '  len: ', len(nms), len(vrs))

            def get_gt_1yr(yr, var):
                gt_1yr = 0
                for mn in range(0, self.nmonth):
                    for j in range(0, self.d1_nlat):
                        gt_1yr = gt_1yr + sum(var[yr * 12 + mn, j, :] * self.garia_d1[j] * 1e-12)
                return gt_1yr

            llst = ['year'] + nms + ['diff']
            df_stt = pd.DataFrame(columns=llst)

            for yr in range(0, self.nyear_nc):
                ln1 = []
                for y in vrs:
                    gt_1yr = get_gt_1yr(yr, y)
                    ln1.append(round(gt_1yr, 2))
                ln1 = [self.years_nc[0] + yr] + ln1 + [round(ln1[1] - ln1[0], 2)]
                df_stt.loc[len(df_stt), :] = ln1
            print(df_stt)
            df_stt.to_csv(self.inv_ncd_dir + 'budget_ch4_gcp25_tot_' + invc + '.txt', sep='\t', index=False)

        # ===========================================================
        print(f'\tRun wrt_flux')
        print(f'\tApr flux from: {self.inv_apr_dir}')
        print(f'\tPst flux from: {self.inv_pst_dir}')

        self.apr_flx_set, self.pst_flx_set = self.get_gcp_flx(self.unp, self.unx, 1)
        for j1, invc in enumerate(self.invcases[0:]):
            apr = self.apr_flx_set[j1]
            pst = self.pst_flx_set[j1]
            # print(apr.shape)
            write_1nc(invc, apr, pst)

    # --- get_gcp_flx
    def get_gcp_flx(self, unp, unx, shape_cor):
        # - read bin
        def read_flx_bin(bfile, skp_yr):
            size_2d = self.d1_nlon * self.d1_nlat * self.nmonth
            with open(bfile, 'rb') as f:
                f.seek(size_2d * skp_yr * 4)
                fl2d = np.fromfile(f, dtype='<f4', count=size_2d * self.nyear_nc
                                   ).reshape(self.nyear_nc, self.nmonth, self.d1_nlat, self.d1_nlon)

            print(fl2d.shape)
            return fl2d

        # - VISIT flux correction
        def cor_visit(flx):
            flx1 = np.full_like(flx, 0)
            for yr in range(0, self.nyear):
                for mn in range(0, self.nmonth):
                    ndys = calendar.monthrange(self.years_nc[0] + yr, mn + 1)[1]
                    flx1[yr, mn, :, :] = flx[yr, mn, :, :] / ndys / 86400.0
            return flx1

        # - flux shape correction [yr, mn, :, :] -> [time, :, :] and scale to monthly
        def cor_shape(flx):
            flx1 = np.zeros([self.nyear_nc * self.nmonth, self.d1_nlat, self.d1_nlon])
            for yr in np.arange(self.nyear_nc):
                for mn in np.arange(self.nmonth):
                    ss = yr * self.nmonth + mn
                    # Kg-CH4/m2/sec ==>> g-CH4/m2/month -GCP recom
                    flx1[ss, :, :] = flx[yr, mn, :, :] * 1000 * 86400 * self.ndays[mn]
            return flx1

        # - flux shape correction: [yr_s : yr_f] -> [yr_s + 1 : yr_f + 1]
        def cor_shape_2(flx):
            flx1 = np.zeros([self.nyear_nc * self.nmonth, self.d1_nlat, self.d1_nlon])
            ss1 = 0
            ss2 = (self.nyear_nc - 1) * self.nmonth
            flx1[ss1:ss2, :, :] = flx[ss1 + 12:ss2 + 12, :, :]
            flx1[ss2:ss2 + 12, :, :] = flx[ss2:ss2 + 12, :, :]
            return flx1

        # - flux value correction
        def cor_val(flx):
            flx1 = flx.copy()
            for mn in np.arange(self.nmonth):
                # Kg-CH4/m2/sec ==>> g-CH4/m2/month
                flx1[:, mn, :, :] = flx[:, mn, :, :] * 1000 * 86400 * self.ndays[mn]
            return flx1

        # ===========================================================
        # - apr
        apr_flx_set = []
        for j1, fl in enumerate(self.flxcases):
            i_file = self.inv_apr_dir + fl + '.grd'
            skp_y = 1
            print('\tPrior     flux reading: ', fl, i_file, skp_y)
            flx = read_flx_bin(i_file, skp_y)
            if shape_cor:
                flx1 = cor_shape(flx)
                flx1 = cor_shape_2(flx1)
            else:
                flx1 = cor_val(flx)
            apr_flx_set.append(flx1)

        # - post
        pst_flx_set = []
        for fl in self.invcases:
            i_file = self.inv_pst_dir + unp + '/' + self.icase + '_' + unx + '_' + fl + '.grd'
            skp_y = 1
            print('\tPosterior flux reading: ', fl, i_file, skp_y)
            flx = read_flx_bin(i_file, skp_y)
            if shape_cor:
                flx1 = cor_shape(flx)
            else:
                flx1 = cor_val(flx)
            pst_flx_set.append(flx1)

        print('\tPrior/Post flux list shape: ', np.array(apr_flx_set).shape, np.array(pst_flx_set).shape)
        return apr_flx_set, pst_flx_set
