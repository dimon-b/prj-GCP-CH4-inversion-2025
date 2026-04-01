# -*- coding: utf-8 -*-
"""
Finalized:  26/xx/xx
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import os
import datetime
import pandas as pd
import matplotlib as mpl

import c1_conflux_nc
import i1_inca_OH
import c3_obs_files
import s0_server
import t0_submit
import v1_valid_surf

def main():
    # --- settings
    pd.set_option('display.expand_frame_repr', False)
    mpl.rcParams['figure.dpi'] = 300

    # --- Convert input fluxes nc -> nc
    run = 0
    if run:
        print('\t *-*-* Start c1_conflux_nc *-*-* ');        c1_conflux_nc.CnvFluxNc()

    # --- Convert input fluxes nc -> gt3
    #     Run FORTRAN f0_fch4_gcp25_nc2gt3_t42.f90

    # --- Forward MIROC4-ACTM with prior

    # --- Copy/Create obs file
    run = 0
    if run:
        print('\t *-*-* Start c3_obs_files *-*-* ');        c3_obs_files.ObsFilesTxt()

    # --- INCA OH
    run = 0
    if run:
        print('\t *-*-* Start i1_inca_OH *-*-* ');        i1_inca_OH.CnvOHNc()

    # --- Run FORTRAN f2_actm_ch4_at_obs_sites.f90
    # --- Run FORTRAN f3_combch4_obs_actm.f90
    # --- Run FORTRAN f4_ch4bu_mon_reg_t42.f90
    # --- Run FORTRAN f5_write_prior_ch4.f90
    # --- Run FORTRAN digi_filt/main_actmstn.f90
    # --- Run FORTRAN f7_write_tdi_obsrvCH4dif.f90

    # --- Run inversion
    #     codeint54

    # ---
    run = 0
    if run:
        print('\t *-*-*-*-* Start Server *-*-*-*-* ');     s0_server.ServerPart()

    # --- Forward MIROC4-ACTM with posterior

    # ---
    run = 1
    if run:
        print('\t *-*-*-*-* Start Submission *-*-*-*-* '); t0_submit.Write2Submit()

    # --- validation conc using ground-based obs
    run = 0
    if run:
        print('\t *-*-*-*-* Start Validation *-*-*-*-* '); v1_valid_surf.ConcValSurf()

    print('\nScript done!')


def est_time(led):
    if led == 0:
        os.system('cls' if os.name == 'nt' else 'clear')
        led = datetime.datetime.now()
        print('\nStart time       : ', led, '\n')
        return led

    else:
        print('\nEnd time        : ', datetime.datetime.now())
        print('Running time    : ', datetime.datetime.now() - led)
        print('Main done')
        return led


if __name__ == "__main__":
    led = est_time(0)
    main()
    led = est_time(led)
