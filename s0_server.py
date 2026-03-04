# -*- coding: utf-8 -*-
"""
Created:    26/01/11 05:20
Project:    GCP/GMB 2025 project [server part]
@author:    Dmitry Belikov
"""

import _set_case
import s1_Burden
import s2_LossCorr
import s3_write_grd
import s4_write_gt3
import s5_LossCorrPlt


class ServerPart(_set_case.SetCase):

    def __init__(self):
        super().__init__()

        run = 1
        if run:
            # --- FORTRAN-based Loss and Burden calculation, Loss Correction implementation
            #print('\n\t *-*-*-*-* Start s1_Burden *-*-*-*-* ');      s1_Burden.LossCorr()

            # --- FORTRAN-based Loss and Burden calculation, Loss Correction implementation
            print('\n\t *-*-*-*-* Start s2_LossCorr *-*-*-*-* ');      s2_LossCorr.LossCorr()

            # --- FORTRAN-based conversion of predicted.sources to 2D .grd files
            #print('\n\t *-*-*-*-* Start s3_write_grd *-*-*-*-* ');     s3_write_grd.Write_grd()

            # --- FORTRAN-based conversion of 2D .grd files to .gt3 format for MORIC4-ACTM post flux run
            #print('\n\t *-*-*-*-* Start s4_write_gt3 *-*-*-*-* ');     s4_write_gt3.Write_gt3()

            # --- check GT, Loss, Loss Correction
            print('\n\t *-*-*-*-* Start s5_LossCorrPlt *-*-*-*-* ');     s5_LossCorrPlt.LossCorrPlot()



