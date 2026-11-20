# -*- coding: utf-8 -*-
"""Calculate atmospheric CH4 burden and chemical loss from forward ACTM output."""

from pathlib import Path

import _set_case
from p6_fortran_tools import full_path, run_fortran


# =============================================================
# === Compute CH4 burden and chemical loss from ACTM hourly GT3
class Burden(_set_case.SetCase):

    def __init__(self, setup=None, prior='p30', uncertainty='ctl', model_dir=None, dry_run=False, overwrite=False):
        if setup is None:
            super().__init__()
        else:
            self.__dict__.update(vars(setup))
        self.run_fort_burdlos(prior=prior, uncertainty=uncertainty, model_dir=model_dir, dry_run=dry_run, overwrite=overwrite)

    # --- run_fort_burdlos
    def run_fort_burdlos(self, prior='p30', uncertainty='ctl', model_dir=None, dry_run=False, overwrite=False):
        model_dir = model_dir or getattr(self, 'inv_mod_dir', None)
        if not model_dir:
            raise ValueError('Provide ACTM GT3 model directory via --model-dir or self.inv_mod_dir')

        output = full_path(self.inv_lsc_dir) / f'frt_burd_{self.icase}_{prior}_{uncertainty}_inv1.txt'
        return run_fortran(
            'pf1_ch4_burloss.f90', self, prior=prior, uncertainty=uncertainty, model_dir=model_dir,
            outputs=[output], output_dirs=[self.inv_lsc_dir],
            dry_run=dry_run, overwrite=overwrite,
        )


# --- Old class name for existing code
LossCorr = Burden
