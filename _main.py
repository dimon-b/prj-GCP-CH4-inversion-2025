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

import c1_convert_flux



def main():
    # --- settings
    pd.set_option('display.expand_frame_repr', False)
    mpl.rcParams['figure.dpi'] = 300

    # --- Convert input fluxes
    run = 1
    if run:
        print('\t\t *-*-* Start c1_convert_flux *-*-* ');        c1_convert_flux.CnvFlux()

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
