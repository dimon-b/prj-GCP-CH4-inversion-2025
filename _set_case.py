# -*- coding: utf-8 -*-
"""
Finalized:  26/xx/xx
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import os
import numpy as np
import calendar
import math


# === set class
class SetCase():
    def __init__(self):
        # --- Inversion / F0-F7 period (inclusive: 27 years, 324 months)
        self.inv_year_s = 1999
        self.inv_year_e = 2025
        self.yr_s = self.inv_year_s       # legacy scripts and Fortran config
        self.yr_e = self.inv_year_e
        self.years = [self.yr_s, self.yr_e]
        self.months = calendar.month_abbr[1:]
        self.nmonth = len(self.months)
        self.nyear = self.yr_e - self.yr_s + 1
        self.ntime = self.nyear * self.nmonth

        # --- Output NetCDF period (NOT the F5/F7 inversion NetCDF)
        self.nc_out_year_s = 2000
        self.nc_out_year_e = 2024
        self.years_nc = [self.nc_out_year_s, self.nc_out_year_e]  # legacy alias
        self.nyear_nc = self.nc_out_year_e - self.nc_out_year_s + 1
        self.ntime_nc = self.nyear_nc * self.nmonth

        if self.ntime != 324 or self.ntime_nc != 300:
            raise ValueError('Expected inversion 324 months and output NetCDF 300 months')

        # --- others
        self.fs = 13
        self.abc = 'abcdefghijklmnopqrstuvwxyz'
        self.mdays = [16, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16]
        self.ndays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        # --- const
        self.R = 6371000

        # --- model
        self.invq = ['Q03', 'Q04', 'Q05', 'Q06', 'Q07']
        self.invcases = ['inv1']
        self.unpcases = ['p30', 'p50', 'p99']
        self.unxcases = ['ctl', 'ux2', 'ux4']
        # self.icase = 's060'
        self.icase = 's057'

        # --- One case for the complete Python workflow
        self.hcase = os.environ.get('CH4_CASE', 'CYC').strip().upper()
        if self.hcase not in ['CYC', 'INCA']:
            raise ValueError('CH4_CASE must be CYC or INCA: ' + self.hcase)

        # --- Fortran sources and executables inside c_gcpv3/
        self.fort_dir = './fortran/'

        # --- initiation
        self.def_path()
        self.def_grid()
        # todo
        # self.def_model()

    # --- def_path
    def def_path(self):
        # ===========================================================
        # - System-dependent directories
        if os.name == 'nt':
            self.flx_inp_dir = 'D:/OneDrive - 国立大学法人千葉大学/prj_GCP_v25/fluxes/gcp2025_flux_inp/'
            self.flx_apr_dir = 'D:/OneDrive - 国立大学法人千葉大学/prj_GCP_v25/fluxes/gcp2025_flux_prior/'
            self.inc_OH_dir = 'D:/OneDrive - 国立大学法人千葉大学/dbase/fluxes/gcp2025_OH/'
            self.obs_dir = 'D:/OneDrive - 国立大学法人千葉大学/prj_apack/obs/'
            actm_dir = 'D:/dbase/ACTM/'
        else:
            self.flx_inp_dir = '/S/data01/G5070/y0715/prj_GCP_v25/fluxes/gcp2025_flux_inp/'
            self.flx_apr_dir = '/S/data01/G5070/y0715/prj_GCP_v25/fluxes/gcp2025_flux_prior/'
            self.inc_OH_dir = '/S/data01/G5070/y0715/prj_GCP_v25/fluxes/gcp2025_OH/'
            self.obs_dir = '/S/data01/G5070/y0715/prj_apack/obs/'
            actm_dir = '../trout_MIROC/'

        # ===========================================================
        # - validation
        self.val_model_cases = {'cyc v2':
                                    {'file': '../submition/v2/MIROC4-ACTM_CH4conc_GMB_SURF_OH_Transcom.nc',
                                     'label': 'cyc v2',
                                     'color': 'b',
                                     'marker': 'o',
                                     },
                                'cyc new':
                                    {'file': '../results2025/CYC/nc_out/MIROC4-ACTM_CH4conc_GMB_SURF_OH_Transcom.nc',
                                     'label': 'cyc new',
                                     'color': 'r',
                                     'marker': '*',
                                     }, }
        self.val_df_name = 'model_obs_stats_inv2025.csv'
        self.val_obs_dir = '../validation/o_post_val/'
        self.val_sts = os.path.join(self.val_obs_dir, '_vall_all.csv')
        self.val_obs_diff_threshold = 100.0

        # ===========================================================
        # - inp fluxes
        self.flx_inp_nc = 'GCP_Prior_CH4_fluxes.nc'

        # - prior fluxes
        self.flx_apr_nc = 'gcp25_inv1.nc'  # total only to convert to gt3
        self.flx_apr_nc_full = 'gcp25_inv1_full.nc'  # full components
        self.flx_apr_grd = self.flx_apr_dir + 'fch4_gcp2025_prior_inv1.grd'
        self.inv_apr_grd = self.flx_apr_grd

        # - inv dirs
        self.inv_dir = '../inv_dir/'

        # --- Stage 1: keep the existing station and GF selection unchanged
        #     Do not reorder the 245 original GF indices.
        self.gf_file = os.environ.get('CH4_GF_FILE', '')
        self.actm_prior_dir = os.environ.get('CH4_PRIOR_MODEL_DIR', '')
        self.actm_post_dir = os.environ.get('CH4_POST_MODEL_DIR', '')

        # --- F2 runtime case input (original model runs)
        f2_defaults = {
            'CYC': '../trout_MIROC/CH4_t42l67_CYC_apr_H_261103/',
            'INCA': '../trout_MIROC/CH4_t42l67_INCA_apr_H_261104/',
        }
        self.f2_model_dir = os.environ.get('CH4_F2_MODEL_DIR', f2_defaults[self.hcase])
        self.f2_output_dir = os.environ.get('CH4_F2_OUTPUT_DIR', self.inv_dir + 'actm/')
        self.f2_output_prefix = 'miroc_ch4sfc_' + self.hcase.lower() + '_'

        # --- F3: combine ACTM concentrations with station observations
        self.f3_model_dir = self.f2_output_dir
        self.f3_obs_dir = self.obsout_dir if hasattr(self, 'obsout_dir') else self.inv_dir + 'obs/'
        self.f3_output_dir = self.inv_dir + 'joined/'
        self.f3_input_prefix = self.f2_output_prefix
        self.f3_output_prefix = 'obsrv_actm_' + self.hcase.lower() + '_'
        # --- Compiled Fortran executables (same binary for CYC/INCA)
        #     build_fortran.sh writes to ./fortran/bin/; optional override stays available.
        self.fort_exe_dir = os.environ.get('CH4_FORTRAN_BIN_DIR', os.path.join(self.fort_dir, 'bin/'))
        def exe_path(name):
            return os.path.abspath(os.path.join(self.fort_exe_dir, name))
        self.exe_f0 = os.environ.get('CH4_EXE_F0', exe_path('f0_fch4_gcp25_nc2gt3_t42'))
        self.exe_f2 = os.environ.get('CH4_EXE_F2', exe_path('f2_actm_ch4_at_obs_sites'))
        self.exe_f3 = os.environ.get('CH4_EXE_F3', exe_path('f3_combch4_obs_actm'))
        self.exe_f4 = os.environ.get('CH4_EXE_F4', exe_path('f4_ch4bu_mon_reg_t42'))
        self.exe_f5 = os.environ.get('CH4_EXE_F5', exe_path('f5_write_prior_ch4'))
        self.exe_f6 = os.environ.get('CH4_EXE_F6', exe_path('f6_nakazawa_filter'))
        self.exe_f7 = os.environ.get('CH4_EXE_F7', exe_path('f7_write_tdi_obsrvCH4dif'))

        # --- F0 converts the common 1999-2025 GCP NetCDF prior to GT3.
        self.f0_input_dir = self.flx_apr_dir
        self.f0_output_dir = self.flx_apr_dir
        self.f0_input_nc = 'gcp25_inv1.nc'
        self.f0_input_case = 'inv1'

        # --- F4/F5 are common to both OH cases if CH4 prior flux is identical.
        #     Keep original shared filenames to preserve downstream compatibility.
        self.f4_flux_file = os.environ.get('CH4_F4_FLUX_FILE', self.flx_apr_dir + 'fch4_gcp25_inv1_t42.gt3')
        self.f4_mask_file = self.inv_dir + 'mask/mask84_t42.bin'
        self.f4_output_file = self.inv_dir + 'uncert/ch4bu_reg54.dat'
        self.f5_input_file = self.f4_output_file
        self.f5_output_prefix = self.inv_dir + 'sourceprior/sourcepriorch4_p'

        # --- F6 Nakazawa digital filter: retain legacy fixed 1999-2025 period.
        self.f6_input_file = self.f3_output_dir + 'obsrv_actm_' + self.hcase.lower() + '_ch4dy.grd'
        self.f6_output_file = self.inv_dir + 'ch4obs_actm_' + self.hcase.lower() + '_fits.bin'

        # --- F7: case-specific outputs avoid overwriting CYC with INCA.
        self.f7_input_file = self.f6_output_file
        self.f7_output_summary = self.inv_dir + 'obsrvCH4_2025_' + self.hcase.lower() + '.txt'
        self.f7_output_prior = self.inv_dir + 'prior_data_' + self.hcase.lower() + '.dif'
        self.f7_output_nc_prefix = self.inv_dir + 'obsrvCH4_' + self.hcase.lower() + '_9925'

        # - OH
        self.inc_OH_inp_dir = self.flx_inp_dir + 'LOSS_FIELDS/OH_INCA/Scaled/'

        # - obs
        self.obspack_dir = (self.obs_dir + 'o_orig/ObsPack/obspack_ch4_1_GLOBALVIEWplus_v8.1_2026-07-13/data/txt/')
        self.wdcgg_dir = self.obs_dir + 'o_orig/WDCGG/'
        self.obsout_dir = '../inv_dir/obs/'
        self.sites_f = self.inv_dir + 'obsrvCH4_60sites.txt'

        # - inv
        self.inv_wrk_dir = '../results2025/' + self.hcase + '/'
        self.inv_run_dir = self.inv_wrk_dir

        # --- Bayesian inversion Python code and control template
        #     The original tdi_ch4 executable is not modified or included in this package.
        #     Copy/link it to ./inversion/tdi_ch4 or set CH4_EXE_INV.
        self.inv_code_dir = os.environ.get('CH4_INV_CODE_DIR', './inversion/')
        self.inv_control_template = os.environ.get('CH4_INV_CONTROL', os.path.join(self.inv_code_dir, 'control.example'))
        self.exe_inv = os.environ.get('CH4_EXE_INV', os.path.join(self.inv_code_dir, 'tdi_ch4'))
        self.inv_lsc_dir = self.inv_run_dir + 'losscorr/'
        self.inv_pst_dir = self.inv_run_dir + 'flux2d/'
        self.flx_pst_dir = self.inv_pst_dir
        self.flx_ncd_dir = self.inv_wrk_dir + 'nc_out/'
        self.inv_ncd_dir = self.flx_ncd_dir
        self.plt_dir = '../plots/'

        # --- Paths for code that expects these attributes
        #     Explicitly supplied model directories only.
        self.inv_mod_dir = self.actm_prior_dir

        # --- Do not activate case-dependent validation comparisons implicitly.
        self.post_model_conc = (
            '../results2025/CYC/nc_out/MIROC4-ACTM_CH4conc_GMB_SURF_OH_Transcom.nc'
            if self.hcase == 'CYC' else
            '../results2025/INCA/nc_out/MIROC4-ACTM_CH4conc_GMB_SURF_OH_INCA.nc'
        )

        # # - concentration
        # self.conc_cases = [
        #     ('apr_cyc', actm_dir + 'CH4_t42l67_CYC_M_260115/'),
        #     ('pst_cyc', actm_dir + 'CH4_t42l67_CYC_M_post_260302/'),
        #     ('apr_inc', actm_dir + 'CH4_t42l67_INCA_M_260307/'),
        #     ('pst_inc', actm_dir + 'CH4_t42l67_INCA_M_post_260311/'),
        #     ('pri_gcp', actm_dir + 'CH4_t42l67_CYC_JRA3Q_M_260115/'),
        #     ('pst_gcp', actm_dir + 'CH4_t42l67_CYC_JRA3Q_M_260131_post/')
        # ]

        # self.conc_cases = [
        #     ('pri_cyc', 'D:/dbase/ACTM_OH/run_ES4/CH4_t42l67_CYC_prior_M_240705/'),
        #     ('pri_iav', 'D:/dbase/ACTM_OH/run_ES4/CH4_t42l67_IAV_prior_M_240707/'),
        #     ('pst_cyc', 'D:/dbase/ACTM_OH/run_Oak/GWG_ch4gcp_OH-clm/'),
        #     ('pst_iav', 'D:/dbase/ACTM_OH/run_Oak/GWG_ch4gcp_OH-iav/')
        # ]

        # self.set_dirs = [
        #     ['D:/dbase/ACTM_OH/run_ES4/CH4_t42l67_CYC_prior_M_240705/',
        #      'D:/dbase/ACTM_OH/run_ES4/CH4_t42l67_IAV_prior_M_240707/'],
        #     ['D:/dbase/ACTM_OH/run_Oak/GWG_ch4gcp_OH-clm/',
        #      'D:/dbase/ACTM_OH/run_Oak/GWG_ch4gcp_OH-iav/']
        # ]

    # --- def_grid
    def def_grid(self):
        self.d1_lons = np.arange(-180 + 0.5, 180, 1)
        self.d1_nlon = len(self.d1_lons)
        self.d1_lats = np.arange(-90 + 0.5, 90, 1)
        self.d1_nlat = len(self.d1_lats)

        dmdeg = 180. / self.d1_nlat
        r0 = 6371.0e3
        angle = np.array([(-89.5 + j * dmdeg) * math.pi / 180.0 for j in range(self.d1_nlat)])
        # self.garia_d1 = math.cos(angle) * (dmdeg * r0 * math.pi / 180.0) ** 2
        self.garia_d1 = np.array([math.cos(x) * (dmdeg * r0 * math.pi / 180.0) ** 2 for x in angle])

        # if not os.path.exists(self.plt_dir):
        #     os.makedirs(self.plt_dir)
        # self.all_sts = self.o_post_dir + '_full_all.csv'
        # self.fin_sts = self.o_post_dir + '_finn_all.csv'
        # self.msk_dir = '../mask/'
        # self.iregs = 54
        # self.tracer = 'ch4'
