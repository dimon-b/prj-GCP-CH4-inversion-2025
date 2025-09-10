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
        # self.yr_s = 1995
        # self.yr_e = 2024
        #
        # # - obs quality
        # self.site_ratio_lim = 0.0
        # self.air_ratio_lim = 0.0

        # --- initiation
        self.def_path()


        # - others
        self.fs = 13
        self.abc = 'abcdefghijklmnopqrstuvwxyz'


    # - path
    def def_path(self):
        # - inp dirs
        self.inp_dir = 'D:/dbase/fluxes/inp_flux_gcp2025/'
        self.inp_flx = 'GCP_Prior_CH4_fluxes.nc'
        # self.o_post_dir = self.wrk_dir + '/o_post/'
        # if not os.path.exists(self.o_post_dir):
        #     os.makedirs(self.o_post_dir)
        # self.plt_dir = '../plots/'
        # if not os.path.exists(self.plt_dir):
        #     os.makedirs(self.plt_dir)
        # self.all_sts = self.o_post_dir + '_full_all.csv'
        # self.fin_sts = self.o_post_dir + '_finn_all.csv'
        # self.msk_dir = '../mask/'
        # self.iregs = 54
        # self.tracer = 'ch4'


