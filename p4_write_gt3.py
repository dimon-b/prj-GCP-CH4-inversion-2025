# -*- coding: utf-8 -*-
"""Conservatively remap posterior 1-degree GRD fluxes to T42 GT3."""

import _set_case
from p6_fortran_tools import full_path, run_fortran


# =============================================================
# === Convert global posterior GRD to T42 GT3
class Write_gt3(_set_case.SetCase):

    def __init__(self, setup=None, prior='p30', uncertainty='ctl', dry_run=False, overwrite=False, lc_mode='file'):
        if setup is None:
            super().__init__()
        else:
            self.__dict__.update(vars(setup))
        self.run_fort_gt3(prior=prior, uncertainty=uncertainty, dry_run=dry_run, overwrite=overwrite, lc_mode=lc_mode)

    # --- run_fort_gt3
    def run_fort_gt3(self, prior='p30', uncertainty='ctl', dry_run=False, overwrite=False, lc_mode='file'):
        model = 'inv1'
        input_dir = full_path(self.inv_pst_dir) / prior
        if lc_mode not in ('file', 'identity'):
            raise ValueError('LC mode must be file or identity')
        suffix = '_raw' if lc_mode == 'identity' else ''
        stem = f'{self.icase}_{uncertainty}_{model}{suffix}'
        input_grd = input_dir / f'{stem}.grd'
        output_gt3 = input_dir / f'{stem}.gt3'

        run_fortran(
            'pf4_ch4_grd2gt3.f90', self, prior=prior, uncertainty=uncertainty,
            inputs=[input_grd], outputs=[output_gt3], output_dirs=[input_dir],
            dry_run=dry_run, overwrite=overwrite, lc_mode=lc_mode,
        )
        return output_gt3
