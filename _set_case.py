# -*- coding: utf-8 -*-
"""
Created/Corrected:    25/07/16 16:54
Project:    ACTMpack project
@author:    Dmitry Belikov
"""
import os
from pylab import *


# === set class
class SetCase():
    def __init__(self):
        # # - period
        self.yr_s = 1999
        self.yr_e = 2024
        # self.yr_s = 1997
        # self.yr_e = 2022
        #
        # # - obs quality
        # self.site_ratio_lim = 0.0
        # self.air_ratio_lim = 0.0

        # --- others
        self.fs = 13
        self.abc = 'abcdefghijklmnopqrstuvwxyz'

        # --- const
        self.R = 6371000

        # --- initiation
        self.def_path()


    # - path
    def def_path(self):
        # - fluxes
        self.inp_dir = 'D:/dbase/fluxes/gcp2025_flux_inp/'
        self.out_dir = 'D:/dbase/fluxes/gcp2025_flux_out/'
        self.inp_flx = 'GCP_Prior_CH4_fluxes.nc'

        # - inv dirs
        self.inv_dir = '../inv_dir/'

        # - obs
        self.obs_dir = 'D:/OneDrive - 国立大学法人千葉大学/prj_apack/obs/'

        # self.obspack_dir = ( self.obs_dir +
        #                     'o_orig/ObsPack/obspack_ch4_1_GLOBALVIEWplus_v7.0_2024-10-29/data/txt/')
        self.wdcgg_dir = (self.obs_dir + '/o_orig/WDCGG/')
        self.obsout_dir = '../inv_dir/obs/'

        self.sites_f = self.inv_dir + 'obsrvCH4_60sites.txt'

        # - model
        self.mod_dir = 'D:/dbase/ACTM/CH4_t42l67_CYC_JRA3Q_M_250920/'
        self.invq = ['Q03', 'Q04', 'Q05', 'Q06', 'Q07']


        # self.plt_dir = '../plots/'
        # if not os.path.exists(self.plt_dir):
        #     os.makedirs(self.plt_dir)
        # self.all_sts = self.o_post_dir + '_full_all.csv'
        # self.fin_sts = self.o_post_dir + '_finn_all.csv'
        # self.msk_dir = '../mask/'
        # self.iregs = 54
        # self.tracer = 'ch4'
