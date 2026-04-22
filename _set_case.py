# -*- coding: utf-8 -*-
"""
Finalized:  26/xx/xx
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import calendar
import math
from pylab import *
import cartopy.crs as ccrs
import cartopy.feature as cf


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
        # self.years_nc = [1999, 2025]
        self.nyear_nc = len(np.arange(self.years_nc[0], self.years_nc[1] + 1, 1))

        # --- others
        self.fs = 13
        self.abc = 'abcdefghijklmnopqrstuvwxyz'
        self.mdays = [16, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16]
        self.ndays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        # --- const
        self.R = 6371000

        # --- initiation
        # self.def_path()
        # self.def_grid()
        # todo
        # self.def_model()

    # --- path
    # def def_path(self):
        # - inp fluxes
        # self.flx_inp_dir = 'D:/OneDrive - 国立大学法人千葉大学/dbase/fluxes/gcp2025_flux_inp/'
        # self.flx_inp_dir = '/S/data01/G5070/y0715/prj_GCP_v25/fluxes/gcp2025_flux_inp/'
        self.flx_inp_dir = 'D:/OneDrive - 国立大学法人千葉大学/prj_GCP_v25/fluxes/gcp2025_flux_inp/'
        self.flx_inp_nc = 'GCP_Prior_CH4_fluxes.nc'

        # - prior fluxes
        self.flx_apr_dir = 'D:/OneDrive - 国立大学法人千葉大学/prj_GCP_v25/fluxes/gcp2025_flux_prior/'
        self.flx_apr_nc = "gcp25_inv1.nc"  # total only to convert to gt3
        self.flx_apr_nc_full = "gcp25_inv1_full.nc"  # full components
        self.flx_apr_grd = self.flx_apr_dir + 'fch4_gcp2025_prior_inv1.grd'

        # - post fluxes
        #self.flx_pst_dir = 'D:/dbase/fluxes/gcp2025_flux_post/'

        # - inv dirs
        self.inv_dir = '../inv_dir/'

        # -
        self.inc_OH_inp_dir = '/S/data01/G5070/y0715/prj_GCP_v25/fluxes/gcp2025_flux_inp/LOSS_FIELDS/OH_INCA/Scaled'
        self.inc_OH_dir = 'D:/OneDrive - 国立大学法人千葉大学/dbase/fluxes/gcp2025_OH/'

        # - obs
        self.obs_dir = 'D:/OneDrive - 国立大学法人千葉大学/prj_apack/obs/'
        self.wdcgg_dir = self.obs_dir + '/o_orig/WDCGG/'
        self.obsout_dir = '../inv_dir/obs/'
        self.sites_f = self.inv_dir + 'obsrvCH4_60sites.txt'

        # - model
        self.invq = ['Q03', 'Q04', 'Q05', 'Q06', 'Q07']
        self.invcases = ['inv1']
        self.unpcases = ['p30', 'p50', 'p99']
        self.unxcases = ['ctl', 'ux2', 'ux4']

        # - inv
        self.icase = 's060'
        self.hcase = 'CYC'
        # self.hcase = 'INCA'
        self.inv_wrk_dir = '../results2025/' + self.hcase + '/'
        # self.inv_mod_dir = '../trout_MIROC/CH4_t42l67_INCA_M_260212/'
        self.inv_mod_dir = '../trout_MIROC/CH4_t42l67_CYC_M_260115/'
        self.inv_run_dir = self.inv_wrk_dir
        self.inv_lsc_dir = self.inv_run_dir + 'losscorr/'
        self.inv_pst_dir = self.inv_run_dir + 'flux2d/'
        self.flx_pst_dir = self.inv_pst_dir
        self.flx_ncd_dir = self.inv_wrk_dir + '/nc_out/'

        self.plt_dir = '../plots/'


        # - concentration
        self.conc_cases = [('apr_cyc', 'D:/dbase/ACTM/CH4_t42l67_CYC_M_260115/'),
                           ('pst_cyc', 'D:/dbase/ACTM/CH4_t42l67_CYC_M_post_260302/'),
                           ('apr_inc', 'D:/dbase/ACTM/CH4_t42l67_INCA_M_260307/'),
                           ('pst_inc', 'D:/dbase/ACTM/CH4_t42l67_INCA_M_post_260311/'),
                           ]

    # ---
    # def def_grid(self):
        self.d1_lons = np.arange(-180 + 0.5, 180, 1)
        self.d1_nlon = len(self.d1_lons)
        self.d1_lats = np.arange(-90 + 0.5, 90, 1)
        self.d1_nlat = len(self.d1_lats)

        dmdeg = 180. / self.d1_nlat
        r0 = 6371.0e3
        angle = np.array([(-89.5 + j * dmdeg) * math.pi / 180.0 for j in range(self.d1_nlat)])
        self.garia_d1 = np.array([math.cos(x) * (dmdeg * r0 * math.pi / 180.0) ** 2 for x in angle])

# - add_map_feat
def add_map_feat(ax, mpl, sites):
    ax.set_extent(mpl[:4], crs=ccrs.PlateCarree())
    ax.xaxis.set_label_text('')
    ax.yaxis.set_label_text('')
    ax.set_xticks(np.arange(round(mpl[0], 0), mpl[1] + 0.001, mpl[4]), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(round(mpl[2], 0), mpl[3] + 0.001, mpl[5]), crs=ccrs.PlateCarree())


    resol = '10m'
    ax.grid(color='silver', linestyle=':', linewidth=0.5)
    # ax.add_feature(cf.NaturalEarthFeature(category='cultural', name='admin_0_boundary_lines_land', scale=resol,
    #                                       edgecolor='grey', facecolor='none', linestyle='--', alpha=0.99, linewidth=0.6,
    #                                       zorder=1))
    ax.add_feature(cf.NaturalEarthFeature('physical', 'lakes', scale=resol, edgecolor='navy',
                                          facecolor='#4292c6',  # '#9ecae1',  # cf.COLORS['water'],
                                          linewidth=0.75, zorder=2))
    ax.add_feature(cf.NaturalEarthFeature('physical', 'rivers_lake_centerlines',
                                          scale=resol, edgecolor='navy', linewidth=0.95, facecolor='none', zorder=1))

    # - Plot sites
    for name, (lat, lon) in sites:
        xlat, ylon = -2.2, 0.2
        if name == 'Karakalpakia':
            xlat, ylon = 0.2, 0.2
        if name == 'Nukus':
            xlat, ylon = -2.0, 0.2
        if name == 'Chimbay':
            xlat, ylon = -2.5, 0.2
        if name == 'Nukus':
            xlat, ylon = -1.9, -0.2
        if name == 'Khiva':
            xlat, ylon = -1.8, 0.2
        if name == 'Samarkand':
            xlat, ylon = -2.2, -0.5
        if name == 'Qarshi':
            xlat, ylon = -1.8, -0.5
        if name == 'Karakul':
            xlat, ylon = -2.2, -0.5
        if name == 'Navoi':
            xlat, ylon = -1.7, 0.2
        if name == 'Dzhangeldy':
            xlat, ylon = -3.5, -0.2
        if name == 'Bukhara':
            xlat, ylon = -2.5, -0.0
        if name == 'Tamdy':
            xlat, ylon = 0.2, 0.2
        if name == 'Fergana':
            xlat, ylon = -1.0, -0.5
        if name == 'Jizzakh':
            xlat, ylon = 0.2, -0.5
        if name == 'Andijan':
            xlat, ylon = 0.2, -0.2
        if name == 'Almalyk':
            xlat, ylon = 0.2, 0.2
        if name == 'Namangan':
            xlat, ylon = -1.2, 0.2
        if name == 'Kokand':
            xlat, ylon = -1.5, -0.5

        ax.plot(lon, lat, marker='o', color='r', markersize=2, transform=ccrs.PlateCarree())
        ax.plot(lon, lat, marker='o', color='w', markersize=1, transform=ccrs.PlateCarree())
        ax.text(lon + xlat, lat + ylon, name, fontsize=9, color='k', transform=ccrs.PlateCarree())