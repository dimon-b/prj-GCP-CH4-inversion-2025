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
import s5_write_nc_tot
import s6_write_nc_cat
import s1_write_nc_loss

class ServerPart(_set_case.SetCase):

    def __init__(self):
        super().__init__()

        run = 0
        if run:
            # --- get Loss and Burden; implement Loss Correction; run s1_ch4_burloss
            print('\n\t *-*-*-*-* Start s1_LossCorr *-*-*-*-* ');      s1_LossCorr.LossCorr()

            # --- convert predicted.sources to 2D grd files by FORTRAN; run s2_write_grd, s3_write_gt3
            print('\n\t *-*-*-*-* Start s2_write_grd *-*-*-*-* ');     s2_write_grd.Write_grd()
            print('\n\t *-*-*-*-* Start s3_write_gt3 *-*-*-*-* ');     s3_write_gt3.Write_gt3()

            # --- check GT, Loss, Loss Correction
            print('\n\t *-*-*-*-* Start s4_LossCorrP *-*-*-*-* ');     s4_LossCorrP.LossCorrPlot()

        # --- write total and category nc
        print('\n\t *-*-*-*-* Start s5_write_nc_tot *-*-*-*-* ');   s5_write_nc_tot.WriteNcTot()
        # print('\n\t *-*-*-*-* Start s1_write_nc_cat *-*-*-*-* ');   s1_write_nc_cat.WriteNcCat()

        # --- write loss for a priori and a posteriori; run s4_ch4_loss.f90
        #print('\n\t *-*-*-*-* Start write nc loss *-*-*-*-* '); s1_write_nc_loss.WriteNcLoss()


