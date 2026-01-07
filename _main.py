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



def main():
    # --- settings
    pd.set_option('display.expand_frame_repr', False)
    mpl.rcParams['figure.dpi'] = 300

    # --- Convert input fluxes nc -> nc
    run = 0
    if run:
        print('\t\t *-*-* Start c1_conflux_nc *-*-* ');        c1_conflux_nc.CnvFluxNc()

    # --- Copy/Create obs file
    run = 1
    if run:
        print('\t\t *-*-* Start c3_obs_files *-*-* ');        c3_obs_files.ObsFilesTxt()

    # todo
    # --- Convert input fluxes nc -> gt3
    #     print('\t\t *-*-* Start c2_conflux_gt3 *-*-* ');        c2_conflux_gt3.CnvFluxGt3()

    # --- Convert input fluxes nc -> gt3 using FORTRAN
    # --- Run model with a priori

    # # --- Create obsrvCH4_*.nc
    # run = 1
    # if run:
    #     print('\t\t *-*-* Start c3_obsmodel_nc *-*-* ');        c3_obsmodel_nc.ObsModNc()


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
