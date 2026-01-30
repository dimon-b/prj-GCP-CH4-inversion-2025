# -*- coding: utf-8 -*-
"""
Created:    26/01/11 05:20
Project:    GCP/GMB 2025 project [server part]
@author:    Dmitry Belikov
"""

import _set_case
import s1_LossCorr
import s2_write_grd
import s3_write_gt3
import s4_LossCorrP


class ServerPart(_set_case.SetCase):

    def __init__(self):
        super().__init__()

        run = 0
        if run:
            # --- FORTRAN-based Loss and Burden calculation, Loss Correction implementation
            print('\n\t *-*-*-*-* Start s1_LossCorr *-*-*-*-* ');      s1_LossCorr.LossCorr()

            # --- FORTRAN-based conversion of predicted.sources to 2D .grd files
            print('\n\t *-*-*-*-* Start s2_write_grd *-*-*-*-* ');     s2_write_grd.Write_grd()

            # --- FORTRAN-based conversion of 2D .grd files to .gt3 format for MORIC4-ACTM post flux run
            print('\n\t *-*-*-*-* Start s3_write_gt3 *-*-*-*-* ');     s3_write_gt3.Write_gt3()

            # --- check GT, Loss, Loss Correction
            print('\n\t *-*-*-*-* Start s4_LossCorrP *-*-*-*-* ');     s4_LossCorrP.LossCorrPlot()



