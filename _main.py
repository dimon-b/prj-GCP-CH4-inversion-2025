# -*- coding: utf-8 -*-
"""
Created/Corrected:    25/09/09 09:09
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import os
import datetime
import pandas as pd
import matplotlib as mpl

import c1_conflux_nc
import c3_obs_files
import s0_server
import t0_submit

def main():
    # --- settings
    pd.set_option('display.expand_frame_repr', False)
    mpl.rcParams['figure.dpi'] = 300

    # --- Convert input fluxes nc -> nc
    run = 0
    if run:
        print('\t\t *-*-* Start c1_conflux_nc *-*-* ');        c1_conflux_nc.CnvFluxNc()

    # --- Convert input fluxes nc -> gt3
    #     Run FORTRAN f0_fch4_gcp25_nc2gt3_t42.f90

    # --- Forward MIROC4-ACTM

    # --- Copy/Create obs file
    run = 0
    if run:
        print('\t\t *-*-* Start c3_obs_files *-*-* ');        c3_obs_files.ObsFilesTxt()

    # --- Run FORTRAN f2_actm_ch4_at_obs_sites.f90
    # --- Run FORTRAN f3_combch4_obs_actm.f90
    # --- Run FORTRAN f4_ch4bu_mon_reg_t42.f90
    # --- Run FORTRAN f5_write_prior_ch4.f90
    # --- Run FORTRAN digi_filt/main_actmstn.f90
    # --- Run FORTRAN f7_write_tdi_obsrvCH4dif.f90

    # --- Run invertion
    #     codeint54

    run = 0
    if run:
        print('\n\t *-*-*-*-* Start Server *-*-*-*-* ');     s0_server.ServerPart()

    run = 1
    if run:
        print('\n\t *-*-*-*-* Start Submission *-*-*-*-* '); t0_submit.Write2Submit()


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
