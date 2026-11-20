# -*- coding: utf-8 -*-
"""GCP/GMB 2025: shared configuration for stages D0-D7.

Fortran operations are named by D-stage and scientific action.
No stage_impl and no dependency on the caller's current directory.
"""

import calendar
import math
import os
from pathlib import Path

import numpy as np


# ============================================================
# === Shared configuration
class SetCase:

    def __init__(self):
        # --- Select CYC or INCA; _main.py and stage tests may set either variable
        self.hcase = os.environ.get('GCP_CASE', os.environ.get('CH4_CASE', 'CYC')).strip().upper()
        if self.hcase not in ('CYC', 'INCA'):
            raise ValueError('GCP_CASE / CH4_CASE must be CYC or INCA: ' + self.hcase)

        # --- Root paths are independent of the launch directory
        self.code_root = Path(__file__).resolve().parent
        self.project_root = self.code_root.parent
        self.fort_dir = str(self.code_root / 'fortran') + os.sep
        self.fort_exe_dir = str(Path(os.environ.get('CH4_FORTRAN_BIN_DIR',
                                                  str(self.code_root / 'fortran' / 'bin'))).expanduser().resolve()) + os.sep

        # --- Inversion period (D2-D7): 1999-2025 inclusive
        self.inv_year_s = 1999
        self.inv_year_e = 2025
        self.yr_s = self.inv_year_s
        self.yr_e = self.inv_year_e
        self.years = [self.yr_s, self.yr_e]
        self.months = list(calendar.month_abbr[1:])
        self.nmonth = 12
        self.nyear = self.yr_e - self.yr_s + 1
        self.ntime = self.nyear * self.nmonth

        # --- Original GCP source period / usual NetCDF diagnostic period
        self.nc_out_year_s = 2000
        self.nc_out_year_e = 2024
        self.years_nc = [self.nc_out_year_s, self.nc_out_year_e]
        self.nyear_nc = self.nc_out_year_e - self.nc_out_year_s + 1
        self.ntime_nc = self.nyear_nc * self.nmonth

        # --- D1 prior-flux period, independent of inversion period
        self.d1_year_s = 1998
        self.d1_year_e = 2026
        self.d1_years = [self.d1_year_s, self.d1_year_e]
        self.d1_nyear = self.d1_year_e - self.d1_year_s + 1
        self.d1_ntime = self.d1_nyear * self.nmonth

        if self.ntime != 324 or self.ntime_nc != 300 or self.d1_ntime != 348:
            raise ValueError('Unexpected period: inversion=324, output=300, D1=348 months')

        # --- Existing project conventions
        self.fs = 13
        self.abc = 'abcdefghijklmnopqrstuvwxyz'
        self.mdays = [16, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16]
        self.ndays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31]
        self.R = 6371000
        self.invq = ['Q03', 'Q04', 'Q05', 'Q06', 'Q07']
        self.invcases = ['inv1']
        self.unpcases = ['p30', 'p50', 'p99']
        self.unxcases = ['ctl', 'ux2', 'ux4']
        self.icase = 's057'

        self.def_path()
        self.def_grid()

    # --------------------------------------------------------
    # --- Normalize project-relative paths and maintain trailing '/' on dirs
    def project_path(self, path):
        path = Path(path).expanduser()
        if not path.is_absolute():
            path = self.project_root / path
        return str(path.resolve())

    def code_path(self, path):
        path = Path(path).expanduser()
        if not path.is_absolute():
            path = self.code_root / path
        return str(path.resolve())

    @staticmethod
    def directory(path):
        return str(path) + ('' if str(path).endswith(('/', '\\')) else os.sep)

    def executable_path(self, step, filename):
        override = os.environ.get('CH4_EXE_' + step, '').strip()
        if override:
            return str(Path(override).expanduser().resolve())
        return str((Path(self.fort_exe_dir) / filename).resolve())

    # ============================================================
    # === All project paths, ordered D0 -> D7
    def def_path(self):
        self.inv_dir = self.directory(self.project_path('inv_dir'))
        self.inv_wrk_dir = self.directory(self.project_path('results2025/' + self.hcase))
        self.inv_run_dir = self.inv_wrk_dir
        self.plt_dir = self.directory(self.project_path('plots'))

        # --- External input archives; preserve original server locations
        if os.name == 'nt':
            user_root = 'D:/OneDrive - 国立大学法人千葉大学/'
            self.flx_inp_dir = self.directory(user_root + 'prj_GCP_v25/fluxes/gcp2025_flux_inp')
            self.flx_apr_dir = self.directory(user_root + 'prj_GCP_v25/fluxes/gcp2025_flux_prior')
            self.inc_OH_dir = self.directory(user_root + 'dbase/fluxes/gcp2025_OH')
            self.obs_dir = self.directory(user_root + 'prj_apack/obs')
        else:
            self.flx_inp_dir = self.directory(self.project_root / 'fluxes/gcp2025_flux_inp')
            self.flx_apr_dir = self.directory(self.project_root / 'fluxes/gcp2025_flux_prior')
            self.inc_OH_dir = self.directory(self.project_root / 'fluxes/gcp2025_OH')
            self.obs_dir = '/S/data01/G5070/y0715/prj_apack/obs/'

        # ========================================================
        # === D0: INCA OH
        self.inc_OH_inp_dir = self.flx_inp_dir + 'LOSS_FIELDS/OH_INCA/Scaled/'

        # ========================================================
        # === D1: CH4 prior fluxes
        self.flx_inp_nc = 'GCP_Prior_CH4_fluxes.nc'
        self.flx_apr_nc = 'gcp25_inv1.nc'
        self.flx_apr_nc_full = 'gcp25_inv1_full.nc'
        self.flx_apr_grd = self.flx_apr_dir + 'fch4_gcp2025_prior_inv1.grd'
        self.inv_apr_grd = self.flx_apr_grd

        # --- Regrid 1-degree CH4 prior to T42
        self.d1_regrid_input_dir = self.flx_apr_dir
        self.d1_regrid_output_dir = self.flx_apr_dir
        self.d1_regrid_input_nc = self.flx_apr_nc
        self.d1_regrid_input_case = 'inv1'
        self.exe_d1_regrid = self.executable_path('D1_REGRID', 'df1_regrid_ch4_t42')
        self.d1_gt3_file = os.environ.get('CH4_D1_GT3_FILE', self.flx_apr_dir + 'fch4_gcp25_inv1_t42.gt3')

        # --- Aggregate 54 regional CH4 fluxes
        self.d1_region_mask_file = self.inv_dir + 'mask/mask84_t42.bin'
        self.d1_region_table_file = self.inv_dir + 'uncert/ch4bu_reg54.dat'
        self.exe_d1_regions = self.executable_path('D1_REGIONS', 'df1_aggregate_ch4_regions')

        # --- Create p30/p50/p99 prior uncertainty NetCDF
        self.d1_uncertainty_input_file = self.d1_region_table_file
        self.d1_uncertainty_output_prefix = self.inv_dir + 'sourceprior/sourcepriorch4_p'
        self.exe_d1_uncertainty = self.executable_path('D1_UNCERTAINTY', 'df1_write_prior_uncertainty')

        # ========================================================
        # === D2: Observations
        self.obspack_dir = self.obs_dir + 'o_orig/ObsPack/obspack_ch4_1_GLOBALVIEWplus_v8.1_2026-07-13/data/txt/'
        self.wdcgg_dir = self.obs_dir + 'o_orig/WDCGG/'
        self.obsout_dir = self.inv_dir + 'obs/'
        self.sites_f = self.inv_dir + 'obsrvCH4_60sites.txt'
        self.gf_file = os.environ.get('CH4_GF_FILE', '')

        # ========================================================
        # === D3: ACTM concentrations and observations
        d3_model_defaults = {
            'CYC': self.directory(self.project_path('trout_MIROC/CH4_t42l67_CYC_apr_H_261119')),
            'INCA': self.directory(self.project_path('trout_MIROC/CH4_t42l67_INCA_apr_H_261119')),
        }

        # --- Sample ACTM at stations
        self.d3_sample_model_dir = os.environ.get('CH4_D3_SAMPLE_MODEL_DIR', d3_model_defaults[self.hcase])
        self.d3_sample_output_dir = os.environ.get('CH4_D3_SAMPLE_OUTPUT_DIR', self.inv_dir + 'actm/')
        self.d3_sample_output_prefix = 'miroc_ch4sfc_' + self.hcase.lower() + '_'
        self.exe_d3_sample = self.executable_path('D3_SAMPLE', 'df3_sample_actm_sites')

        # --- Join station observations with ACTM
        self.d3_combine_model_dir = self.d3_sample_output_dir
        self.d3_combine_obs_dir = self.obsout_dir
        self.d3_combine_output_dir = self.inv_dir + 'joined/'
        self.d3_combine_input_prefix = self.d3_sample_output_prefix
        self.d3_combine_output_prefix = 'obsrv_actm_' + self.hcase.lower() + '_'
        self.exe_d3_combine = self.executable_path('D3_COMBINE', 'df3_combine_obs_actm')

        # --- Filter OBS / ACTM daily series
        self.d3_filter_input_file = self.d3_combine_output_dir + self.d3_combine_output_prefix + 'ch4dy.grd'
        self.d3_filter_output_file = self.inv_dir + 'ch4obs_actm_' + self.hcase.lower() + '_fits.bin'
        self.exe_d3_filter = self.executable_path('D3_FILTER', 'df3_filter_ch4_nakazawa')

        # --- Generate TDI observations and uncertainties
        self.d3_observations_input_file = self.d3_filter_output_file
        self.d3_observations_summary = self.inv_dir + 'obsrvCH4_2025_' + self.hcase.lower() + '.txt'
        self.d3_observations_prior = self.inv_dir + 'prior_data_' + self.hcase.lower() + '.dif'
        self.d3_observations_nc_prefix = self.inv_dir + 'obsrvCH4_' + self.hcase.lower() + '_9925'
        self.exe_d3_observations = self.executable_path('D3_OBSERVATIONS', 'df3_write_tdi_observations')

        # ========================================================
        # === D4: Bayesian TDI inversion, existing Fortran is not moved
        self.d4_tdi_dir = self.directory(self.inv_dir + 'tdi_inv')
        self.d4_files4run_dir = self.directory(self.d4_tdi_dir + 'files4run')
        self.d4_gf_dir = self.directory(self.inv_dir + 'gmats')
        self.inv_code_dir = self.d4_tdi_dir
        self.inv_control_template = os.environ.get('CH4_INV_CONTROL', self.d4_tdi_dir + 'control.example')
        self.exe_inv = os.environ.get('CH4_EXE_INV', self.d4_tdi_dir + 'tdi_ch4')
        # --- Start with one experiment; expand lists for the full 3 x 3 design
        self.d4_priors = ['p30']
        self.d4_uncs = ['ctl']
        self.d4_models = ['inv1']
        self.d4_diagnostics_dir = self.inv_wrk_dir + 'diagnostics/d4/'

        # ========================================================
        # === D5: Post-processing / posterior model fields
        self.inv_lsc_dir = self.inv_run_dir + 'losscorr/'
        self.inv_pst_dir = self.inv_run_dir + 'flux2d/'
        self.flx_pst_dir = self.inv_pst_dir
        self.flx_ncd_dir = self.inv_wrk_dir + 'nc_out/'
        self.inv_ncd_dir = self.flx_ncd_dir
        self.actm_prior_dir = os.environ.get('CH4_PRIOR_MODEL_DIR', self.d3_sample_model_dir)
        self.actm_post_dir = os.environ.get('CH4_POST_MODEL_DIR', '')
        self.inv_mod_dir = self.actm_prior_dir

        # ========================================================
        # === D6: Submission
        # --- Current d6_submit.py uses submission code outside SetCase.
        # --- Keep the configured directory, do not invent filenames.
        self.submit_dir = self.directory(self.project_path('submition'))

        # ========================================================
        # === D7: Validation
        self.val_model_cases = {
            'cyc v2': {
                'file': self.project_path('submition/v2/MIROC4-ACTM_CH4conc_GMB_SURF_OH_Transcom.nc'),
                'label': 'cyc v2', 'color': 'b', 'marker': 'o',
            },
            'cyc new': {
                'file': self.project_path('results2025/CYC/nc_out/MIROC4-ACTM_CH4conc_GMB_SURF_OH_Transcom.nc'),
                'label': 'cyc new', 'color': 'r', 'marker': '*',
            },
        }
        self.val_df_name = 'model_obs_stats_inv2025.csv'
        self.val_obs_dir = self.directory(self.project_path('validation/o_post_val'))
        self.val_sts = self.val_obs_dir + '_vall_all.csv'
        self.val_obs_diff_threshold = 100.0
        self.post_model_conc = (
            self.project_path('results2025/CYC/nc_out/MIROC4-ACTM_CH4conc_GMB_SURF_OH_Transcom.nc')
            if self.hcase == 'CYC' else
            self.project_path('results2025/INCA/nc_out/MIROC4-ACTM_CH4conc_GMB_SURF_OH_INCA.nc')
        )

    # ============================================================
    # === D1: one-degree lat/lon grid
    def def_grid(self):
        self.d1_lons = np.arange(-179.5, 180.0, 1.0)
        self.d1_nlon = len(self.d1_lons)
        self.d1_lats = np.arange(-89.5, 90.0, 1.0)
        self.d1_nlat = len(self.d1_lats)
        dmdeg = 180.0 / self.d1_nlat
        angle = np.deg2rad(self.d1_lats)
        self.garia_d1 = np.cos(angle) * (dmdeg * self.R * math.pi / 180.0) ** 2
