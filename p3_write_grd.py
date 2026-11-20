# -*- coding: utf-8 -*-
"""Apply TDI posterior corrections and annual loss correction to prior GRD fluxes."""

from pathlib import Path

import _set_case
from p6_fortran_tools import full_path, run_fortran


# =============================================================
# === Convert predicted.sources to global posterior GRD
class Write_grd(_set_case.SetCase):

    def __init__(self, setup=None, prior='p30', uncertainty='ctl', dry_run=False, overwrite=False, lc_mode='file'):
        if setup is None:
            super().__init__()
        else:
            self.__dict__.update(vars(setup))
        self.run_fort_grd(prior=prior, uncertainty=uncertainty, dry_run=dry_run, overwrite=overwrite, lc_mode=lc_mode)

    # --- run_fort_grd
    def run_fort_grd(self, prior='p30', uncertainty='ctl', dry_run=False, overwrite=False, lc_mode='file'):
        model = 'inv1'
        input_dir = full_path(self.inv_wrk_dir)
        output_dir = full_path(self.inv_pst_dir)
        predicted = input_dir / f'{self.icase}_{prior}' / f'{self.icase}_{uncertainty}_{model}' / 'predicted.sources'
        prior_grd = full_path(self.flx_apr_dir) / 'fch4_gcp2025_prior_inv1.grd'
        loss_factor = full_path(self.inv_lsc_dir) / f'gcp_LC_{self.icase}_{prior}_{uncertainty}_{model}.txt'
        mask_dir = full_path(self.inv_dir) / 'mask'
        if lc_mode not in ('file', 'identity'):
            raise ValueError('LC mode must be file or identity')
        suffix = '_raw' if lc_mode == 'identity' else ''
        output_grd = output_dir / prior / f'{self.icase}_{uncertainty}_{model}{suffix}.grd'
        budget = output_dir / f'check_GT_{self.icase}{suffix}_{prior}_{uncertainty}.txt'
        regional = output_dir / f'{self.icase}_{prior}_{uncertainty}_{model}{suffix}_fluxlnreg.dat'

        run_fortran(
            'pf3_ch4_grd.f90', self, prior=prior, uncertainty=uncertainty,
            inputs=[predicted, prior_grd, mask_dir / 'mask80.bin', mask_dir / 'ch4_regn_pulse.grd'] + ([loss_factor] if lc_mode == 'file' else []),
            outputs=[output_grd, budget, regional],
            output_dirs=[output_dir, output_grd.parent],
            dry_run=dry_run, overwrite=overwrite, lc_mode=lc_mode,
        )
        if not dry_run:
            expected_bytes = (self.yr_e - self.yr_s + 1) * 12 * 360 * 180 * 4
            if output_grd.stat().st_size != expected_bytes:
                raise RuntimeError(f'Invalid posterior GRD size: {output_grd.stat().st_size}, expected {expected_bytes}')
        return output_grd
