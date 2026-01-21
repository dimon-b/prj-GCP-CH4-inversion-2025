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

# import b_sub_flux
import _set_case


class WriteNcCat(_set_case.SetCase):

    def __init__(self):
        super().__init__()
        self.invcases = ['inv1', 'inv7', 'vCao', 'viWH']
        # self.invcindx = [0, 4, 5]
        self.flxcases = ['gcp2021_v2_soil0_inv1', 'gcp2021_v2_soil0_inv7', 'gcp2021_v2_soil0_vCao']
        self.invcases = ['viWH']
        #self.invcindx = [5]
        self.flxcases = ['gcp2021_v2_soil0_viWH']

        self.unp = 'p30'
        self.unx = 'ux4'
        self.wrt_flux_cat()

    # --- write flux categories
    def wrt_flux_cat(self, ):
        print(f'\tRun wrt_flux_cat')
        print(f'\tApr flux from: {self.aprf_dir}')
        print(f'\tPst flux from: {self.pstf_dir}')

        # --- total of apr and post
        self.apr_flx_set, self.pst_flx_set = b_sub_flux.get_gcp_flx(self, self.unp, self.unx, 1)
        for j1, inv in enumerate(self.invcases):
            print('\tConversion now:', j1, inv)
            self.inv = inv
            # self.idx = self.invcindx[j1]
            apr = self.apr_flx_set[j1]
            pst = self.pst_flx_set[j1]
            self.write_1nc_cat(apr, pst)

    # --- Global total cat
    def cat_GTotal(self, nms, vrs, invc, cs):

        print(f'\tGtotal for flux: ', invc, '  category len: ', len(nms), len(vrs))

        def get_gt_1yr(yr, var):
            gt_1yr = 0
            for mn in range(0, self.nmonth):
                for j in range(0, self.d1_nlat):
                    gt_1yr = gt_1yr + sum(var[yr*12 + mn, j, :]*self.garia_d1[j]*1e-12)
            return gt_1yr

        # llst = ['year'] + nms + ['sum_apr']
        llst = ['year'] + nms + ['sum_apr/pst'] + ['diff']
        df_stt = pd.DataFrame(columns=llst)

        for yr in range(0, self.nyear):
            ln1 = []
            for y in vrs:
                gt_1yr = get_gt_1yr(yr, y)
                ln1.append(round(gt_1yr, 2))
            # ln1 = [self.years[0] + yr] + ln1 + [round(sum(ln1[2:]), 2)]
            ln1 = [self.years[0] + yr] + ln1 + [round(sum(ln1[2:]), 2)] + [round(ln1[0] - sum(ln1[2:]), 2)]
            df_stt.loc[len(df_stt), :] = ln1
        print(df_stt)
        df_stt.to_csv(self.nc_dir + 'budget_ch4_gcp21_cat_' + invc + '_' + cs + '.txt', sep='\t', index=False)

    # --- Global total cat
    def ctg_scale(self, vrs_inp, apr, pst):
        vrs_out = []
        apr1 = np.where(apr == 0, 1.0, apr)
        scl = np.array(pst)/np.array(apr1)
        scl = np.where(scl > 10, 10.0, scl)
        scl = np.where(scl < -10, -10.0, scl)
        for var in vrs_inp[:]:
            apr_flx_cor = var * scl
            vrs_out.append(apr_flx_cor)
        return vrs_out

    # --- Write the netcdf file
    def write_1nc_cat(self, apr, pst):

        # --- apr cat
        ctg_apr, fl_sol, ctg_nms = b_sub_flux.get_flx_categor(self, self.inv)
        #ctg_apr = [wet_apr, bb_apr, biofuel_apr, oilgas_apr, coal_apr, agri_apr, waste_apr, geol_apr,
        #           termite_apr, oce_apr, soil_apr]
        ctg_pst = self.ctg_scale(ctg_apr, apr, pst)

        # --- check G-total for category names and vars
        ctg_nms2 = ['apr', 'pst'] + ctg_nms
        ctg_apr2 = [apr, pst] + ctg_apr
        self.cat_GTotal(ctg_nms2, ctg_apr2, self.inv, 'ap')
        ctg_pst2 = [apr, pst] + ctg_pst
        self.cat_GTotal(ctg_nms2, ctg_pst2, self.inv, 'ps')
        # exit()

        fout = self.nc_dir + 'MIROC4-ACTM_12cat_' + self.inv + 'SURF.nc'
        nc = netCDF4.Dataset(fout, 'w', format='NETCDF4')
        nc.description = 'Monthly category wise flux for GCP-CH4, 2021. The results are produced at JAMSTEC, Japan. For details please contact at prabir@jamstec.go.jp and d.belikov@chiba-u.jp.'
        nc.Institution = 'Japan Agency for Marine-Earth Science and Technology (JAMSTEC)'
        nc.Contact = 'd.belikov@chiba-u.jp and prabir@jamstec.go.jp'
        nc.Disclaimer = 'This data is created for GCP-CH4 2021.'
        nc.History = 'Created on January 2022 by Dmitry Belikov'
        nc.Additional_Info = 'MIROC4.0-based Atmospheric Chemistry-Transport Model (MIROC4-ACTM) is used for forward simulations (Prabir et al., 2018) and Time-independent Bayesian inversion is used for optimising flux (Patra et al., 2016)'

        # --- set dimensions
        time_dim = nc.createDimension('time', None)
        lat_dim = nc.createDimension('lat', self.d1_nlat)
        lon_dim = nc.createDimension('lon', self.d1_nlon)

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
        cellarea = nc.createVariable('carea', 'f4', ('lat', 'lon'))

        fch4_tot_prior = nc.createVariable('fch4_tot_prior', 'f4', ('time', 'lat', 'lon'))
        fch4_tot_post = nc.createVariable('fch4_tot_post', 'f4', ('time', 'lat', 'lon'))
        fch4_wet_prior = nc.createVariable('fch4_wet_prior', 'f4', ('time', 'lat', 'lon'))
        fch4_wet_post = nc.createVariable('fch4_wet_post', 'f4', ('time', 'lat', 'lon'))
        fch4_bb_prior = nc.createVariable('fch4_bb_prior', 'f4', ('time', 'lat', 'lon'))
        fch4_bb_post = nc.createVariable('fch4_bb_post', 'f4', ('time', 'lat', 'lon'))
        fch4_biofuel_prior = nc.createVariable('fch4_biofuel_prior', 'f4', ('time', 'lat', 'lon'))
        fch4_biofuel_post = nc.createVariable('fch4_biofuel_post', 'f4', ('time', 'lat', 'lon'))
        fch4_oilgas_prior = nc.createVariable('fch4_oilgas_prior', 'f4', ('time', 'lat', 'lon'))
        fch4_oilgas_post = nc.createVariable('fch4_oilgas_post', 'f4', ('time', 'lat', 'lon'))
        fch4_coal_prior = nc.createVariable('fch4_coal_prior', 'f4', ('time', 'lat', 'lon'))
        fch4_coal_post = nc.createVariable('fch4_coal_post', 'f4', ('time', 'lat', 'lon'))
        fch4_agri_prior = nc.createVariable('fch4_agri_prior', 'f4', ('time', 'lat', 'lon'))
        fch4_agri_post = nc.createVariable('fch4_agri_post', 'f4', ('time', 'lat', 'lon'))
        fch4_waste_prior = nc.createVariable('fch4_waste_prior', 'f4', ('time', 'lat', 'lon'))
        fch4_waste_post = nc.createVariable('fch4_waste_post', 'f4', ('time', 'lat', 'lon'))
        fch4_termite_prior = nc.createVariable('fch4_termite_prior', 'f4', ('time', 'lat', 'lon'))
        fch4_termite_post = nc.createVariable('fch4_termite_post', 'f4', ('time', 'lat', 'lon'))
        fch4_geol_prior = nc.createVariable('fch4_geol_prior', 'f4', ('time', 'lat', 'lon'))
        fch4_geol_post = nc.createVariable('fch4_geol_post', 'f4', ('time', 'lat', 'lon'))
        fch4_oce_prior = nc.createVariable('fch4_oce_prior', 'f4', ('time', 'lat', 'lon'))
        fch4_oce_post = nc.createVariable('fch4_oce_post', 'f4', ('time', 'lat', 'lon'))
        fch4_soils_prior = nc.createVariable('fch4_soils_prior', 'f4', ('time', 'lat', 'lon'))
        fch4_soils_post = nc.createVariable('fch4_soils_post', 'f4', ('time', 'lat', 'lon'))

        # --- set data
        lat[:] = self.d1_lats
        lon[:] = self.d1_lons
        cellarea[:] = self.garia_d2

        fch4_tot_prior[:, :, :] = apr[:, :, :]
        fch4_tot_post[:, :, :] = pst[:, :, :]

        fch4_wet_prior[:, :, :] = ctg_apr[0]
        fch4_wet_post[:, :, :] = ctg_pst[0]
        fch4_bb_prior[:, :, :] = ctg_apr[1]
        fch4_bb_post[:, :, :] = ctg_pst[1]
        fch4_biofuel_prior[:, :, :] = ctg_apr[2]
        fch4_biofuel_post[:, :, :] = ctg_pst[2]
        fch4_oilgas_prior[:, :, :] = ctg_apr[3]
        fch4_oilgas_post[:, :, :] = ctg_pst[3]
        fch4_coal_prior[:, :, :] = ctg_apr[4]
        fch4_coal_post[:, :, :] = ctg_pst[4]
        fch4_agri_prior[:, :, :] = ctg_apr[5]
        fch4_agri_post[:, :, :] = ctg_pst[5]
        fch4_waste_prior[:, :, :] = ctg_apr[6]
        fch4_waste_post[:, :, :] = ctg_pst[6]
        fch4_geol_prior[:, :, :] = ctg_apr[7]
        fch4_geol_post[:, :, :] = ctg_pst[7]
        fch4_termite_prior[:, :, :] = ctg_apr[8]
        fch4_termite_post[:, :, :] = ctg_pst[8]
        fch4_oce_prior[:, :, :] = ctg_apr[9]
        fch4_oce_post[:, :, :] = ctg_pst[9]
        fch4_soils_prior[:, :, :] = fl_sol #ctg_apr[10]
        fch4_soils_post[:, :, :] = fl_sol #ctg_apr[10]

        # --- date-time
        dd = []
        for ym in np.arange(self.years[0], self.years[1] + 1, 1):
            for mm in np.arange(1, 13, 1):
                dd.append(datetime(int(ym), int(mm), int(self.mdays[mm - 1])))
        dt = date2num(dd, units=time.units, calendar=time.calendar)
        time[:] = dt[:]

        # --- units
        cellarea.units = 'm2'
        fch4_tot_prior.units = 'g-CH4/m2/month'
        fch4_tot_post.units = 'g-CH4/m2/month'
        fch4_wet_prior.units = 'g-CH4/m2/month'
        fch4_wet_post.units = 'g-CH4/m2/month'
        fch4_bb_prior.units = 'g-CH4/m2/month'
        fch4_bb_post.units = 'g-CH4/m2/month'
        fch4_biofuel_prior.units = 'g-CH4/m2/month'
        fch4_biofuel_post.units = 'g-CH4/m2/month'
        fch4_oilgas_prior.units = 'g-CH4/m2/month'
        fch4_oilgas_post.units = 'g-CH4/m2/month'
        fch4_coal_prior.units = 'g-CH4/m2/month'
        fch4_coal_post.units = 'g-CH4/m2/month'
        fch4_agri_prior.units = 'g-CH4/m2/month'
        fch4_agri_post.units = 'g-CH4/m2/month'
        fch4_waste_prior.units = 'g-CH4/m2/month'
        fch4_waste_post.units = 'g-CH4/m2/month'
        fch4_geol_prior.units = 'g-CH4/m2/month'
        fch4_geol_post.units = 'g-CH4/m2/month'
        fch4_termite_prior.units = 'g-CH4/m2/month'
        fch4_termite_post.units = 'g-CH4/m2/month'
        fch4_oce_prior.units = 'g-CH4/m2/month'
        fch4_oce_post.units = 'g-CH4/m2/month'
        fch4_soils_prior.units = 'g-CH4/m2/month'
        fch4_soils_post.units = 'g-CH4/m2/month'

        print('*** SUCCESS writing for ', self.inv)
        nc.close()
