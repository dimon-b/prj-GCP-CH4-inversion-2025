# -*- coding: utf-8 -*-
"""
Created:    26/01/11 05:20
Project:    GCP/GMB 2025 project [server part]
@author:    Dmitry Belikov
"""

import _set_case
import t1_write_nc_flux
import t3_write_nc_sink
import t4_write_3D_conc
import t5_write_comparison


class Write2Submit(_set_case.SetCase):

    def __init__(self):
        super().__init__()

        # --- write total and category fluxes
        print('\n\t *-*-*-*-* Start t1_write_nc_tot *-*-*-*-* ');        t1_write_nc_flux.WriteNcFlux()

        # --- write loss for a priori and a posteriori
        # print('\n\t *-*-*-*-* Start t3_write_nc_sink *-*-*-*-* ');        t3_write_nc_sink.WriteNcLoss()

        # --- write 3D conc
        # print('\n\t *-*-*-*-* Start t4_write_3D_conc *-*-*-*-* ');   t4_write_3D_conc.Write3dConc()

        # --- write comparison
        # print('\n\t *-*-*-*-* Start t5_write_comparison *-*-*-*-* ');   t5_write_comparison.WriteNcComp()
