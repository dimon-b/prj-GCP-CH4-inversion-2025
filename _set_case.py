# -*- coding: utf-8 -*-
"""
Finalized:  26/xx/xx
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import os
import calendar
import math
from pylab import *


# === set class
class SetCase():
    def __init__(self):
        # --- period
        self.yr_s = 1999
        self.yr_e = 2025
        self.years = [self.yr_s, self.yr_e]
        self.months = calendar.month_abbr[1:]
        self.nmonth = len(self.months)
        self.nyear = len(np.arange(self.years[0], self.years[1], 1))
        self.ntime = self.nyear * self.nmonth

        # --- nc output period
        self.years_nc = [2000, 2024]
        self.nyear_nc = len(np.arange(self.years_nc[0], self.years_nc[1] + 1, 1))

        # --- others
        self.fs = 13
        self.abc = 'abcdefghijklmnopqrstuvwxyz'
        self.mdays = [16, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16]
        self.ndays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        # --- const
        self.R = 6371000

        # --- initiation
        self.def_path()
        self.def_grid()
        # todo
        # self.def_model()

    # - path
    def def_path(self):
        # - inp fluxes
        self.flx_inp_dir = 'D:/dbase/fluxes/gcp2025_flux_inp/'
        self.flx_inp_dir = '/S/data01/G5070/y0715/prj_GCP_v3/fluxes/gcp2025_flux_inp/'
        self.flx_inp_nc = 'GCP_Prior_CH4_fluxes.nc'

        # - prior fluxes
        self.flx_apr_dir = 'D:/dbase/fluxes/gcp2025_flux_prior/'
        self.flx_apr_nc = "gcp25_inv1.nc"  # total only to convert to gt3
        self.flx_apr_nc_full = "gcp25_inv1_full.nc"  # full components

        # - inv dirs
        self.inv_dir = '../inv_dir/'

        # -
        self.inc_OH_inp_dir = 'D:/dbase/fluxes/gcp2025_flux_inp/LOSS_FIELDS/OH_INCA/Scaled/'
        self.inc_OH_dir = 'D:/dbase/fluxes/gcp2025_OH/'

        # self.inc_OH_inp_dir = '/S/data01/G5070/y0715/prj_GCP_v3/fluxes/gcp2025_flux_inp/LOSS_FIELDS/OH_INCA/Scaled'

        # - obs
        self.obs_dir = 'D:/OneDrive - 国立大学法人千葉大学/prj_apack/obs/'
        # self.obspack_dir = ( self.obs_dir + 'o_orig/ObsPack/')
        self.wdcgg_dir = self.obs_dir + '/o_orig/WDCGG/'
        self.obsout_dir = '../inv_dir/obs/'
        self.sites_f = self.inv_dir + 'obsrvCH4_60sites.txt'

        # - model
        self.mod_dir = 'D:/dbase/ACTM/CH4_t42l67_CYC_JRA3Q_M_250920/'
        self.invq = ['Q03', 'Q04', 'Q05', 'Q06', 'Q07']
        self.invcases = ['inv1']
        self.unpcases = ['p30', 'p50', 'p99']
        self.unxcases = ['ctl', 'ux2', 'ux4']

        # - inv
        self.icase = 's060'
        self.inv_wrk_dir = '../results2025/'
        self.inv_run_dir = self.inv_wrk_dir
        self.inv_lss_dir = self.inv_run_dir + 'losscorr/'
        self.inv_apr_dir = self.inv_run_dir + 'priors/'
        self.inv_pst_dir = self.inv_run_dir + 'flux2d/'
        self.inv_ncd_dir = self.inv_run_dir + 'nc_out/'
        self.plt_dir = '../plots/'
        self.inv_apr_grd = self.inv_apr_dir + 'fch4_gcp2025_prior_inv1.grd'

        # - concentration
        self.conc_cases = [('pri_cyc', 'D:/dbase/ACTM_OH/run_ES4/CH4_t42l67_CYC_prior_M_240705/'),
                           ('pri_iav', 'D:/dbase/ACTM_OH/run_ES4/CH4_t42l67_IAV_prior_M_240707/'),
                           ('pst_cyc', 'D:/dbase/ACTM_OH/run_Oak/GWG_ch4gcp_OH-clm/'),
                           ('pst_iav', 'D:/dbase/ACTM_OH/run_Oak/GWG_ch4gcp_OH-iav/')
                           ]

        self.conc_cases = [('pri_gcp', 'D:/dbase/ACTM/CH4_t42l67_CYC_JRA3Q_M_260115/'),
                           ('pst_gcp', 'D:/dbase/ACTM/CH4_t42l67_CYC_JRA3Q_M_260131_post/'),
                           ]

        # self.set_dirs = [['D:/dbase/ACTM_OH/run_ES4/CH4_t42l67_CYC_prior_M_240705/',
        #                  'D:/dbase/ACTM_OH/run_ES4/CH4_t42l67_IAV_prior_M_240707/'],
        #                  ['D:/dbase/ACTM_OH/run_Oak/GWG_ch4gcp_OH-clm/',
        #                  'D:/dbase/ACTM_OH/run_Oak/GWG_ch4gcp_OH-iav/']]

    #
    def def_grid(self):
        self.d1_lons = np.arange(-180 + 0.5, 180, 1)
        self.d1_nlon = len(self.d1_lons)
        self.d1_lats = np.arange(-90 + 0.5, 90, 1)
        self.d1_nlat = len(self.d1_lats)

        dmdeg = 180. / self.d1_nlat
        r0 = 6371.0e3
        angle = np.array([(-89.5 + j * dmdeg) * math.pi / 180.0 for j in range(self.d1_nlat)])
        # self.garia_d1 = math.cos(angle)*(dmdeg*r0*math.pi/180.0)**2
        self.garia_d1 = np.array([math.cos(x) * (dmdeg * r0 * math.pi / 180.0) ** 2 for x in angle])

        # if not os.path.exists(self.plt_dir):
        #     os.makedirs(self.plt_dir)
        # self.all_sts = self.o_post_dir + '_full_all.csv'
        # self.fin_sts = self.o_post_dir + '_finn_all.csv'
        # self.msk_dir = '../mask/'
        # self.iregs = 54
        # self.tracer = 'ch4'
