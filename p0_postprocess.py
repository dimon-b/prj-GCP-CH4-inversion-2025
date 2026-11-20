# -*- coding: utf-8 -*-
"""Run CH4 posterior processing (burden, loss correction, GRD, GT3, plots)."""

import argparse

import _set_case
import p1_burden
import p2_loss_corr
import p3_write_grd
import p4_write_gt3
import p5_loss_corr_plot


# =============================================================
# === Post-processing driver
class ServerPart(_set_case.SetCase):

    def __init__(self, case=None, icase=None, prior='p30', uncertainty='ctl',
                 steps=None, model_dir=None, post_nc=None, lc_site='SYO',
                 dry_run=False, overwrite=False, missing_2025='strict', lc_mode='file'):
        super().__init__()

        if case:
            self.hcase = case.upper()
            if self.hcase not in ['CYC', 'INCA']:
                raise ValueError('Case must be CYC or INCA')
        if icase:
            self.icase = icase
        if case:
            self.def_path()

        # --- Legacy default ran loss correction and its plots only
        if steps is None:
            steps = ['loss', 'plot']
        if isinstance(steps, str):
            steps = [name.strip().lower() for name in steps.split(',') if name.strip()]
        if 'all' in steps:
            steps = ['burden', 'loss', 'grd', 'gt3', 'plot']
        if not steps or any(s not in ('burden', 'loss', 'grd', 'gt3', 'plot') for s in steps):
            raise ValueError('Valid steps: burden,loss,grd,gt3,plot or all')

        print('\n=== POST process ==========================================')
        print('    Case      :', self.hcase, self.icase)
        print('    Prior/unc :', prior, uncertainty)
        print('    Steps     :', ','.join(steps))
        print('    Dry-run   :', dry_run)
        print('    Overwrite :', overwrite)
        print('    LC mode   :', lc_mode)
        if lc_mode == 'identity' and ('loss' in steps or 'plot' in steps):
            raise ValueError('Identity LC applies only to grd/gt3; run burden/loss/plot separately')

        for step in steps:
            if step == 'burden':
                p1_burden.Burden(self, prior=prior, uncertainty=uncertainty, model_dir=model_dir, dry_run=dry_run, overwrite=overwrite)
            if step == 'loss':
                p2_loss_corr.LossCorr(self, prior=prior, uncertainty=uncertainty, site=lc_site, post_model_conc=post_nc,
                                       dry_run=dry_run, overwrite=overwrite, missing_2025=missing_2025)
            if step == 'grd':
                p3_write_grd.Write_grd(self, prior=prior, uncertainty=uncertainty,
                                       dry_run=dry_run, overwrite=overwrite, lc_mode=lc_mode)
            if step == 'gt3':
                p4_write_gt3.Write_gt3(self, prior=prior, uncertainty=uncertainty,
                                       dry_run=dry_run, overwrite=overwrite, lc_mode=lc_mode)
            if step == 'plot':
                p5_loss_corr_plot.LossCorrPlot(self, prior=prior, uncertainty=uncertainty,
                                               site=lc_site, dry_run=dry_run, overwrite=overwrite)


# --- Direct command-line entry point (no _main.py changes required)
def main():
    parser = argparse.ArgumentParser(description='GCP2025 CH4 posterior processing')
    parser.add_argument('--case', choices=['CYC', 'INCA'], required=True)
    parser.add_argument('--icase', default='s057')
    parser.add_argument('--prior', default='p30')
    parser.add_argument('--unc', default='ctl')
    parser.add_argument('--steps', default='loss,plot')
    parser.add_argument('--model-dir', default=None)
    parser.add_argument('--post-nc', default=None)
    parser.add_argument('--lc-site', default='SYO')
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--lc-mode', choices=['file', 'identity'], default='file')
    parser.add_argument('--missing-2025', choices=['strict', 'identity'], default='strict')
    args = parser.parse_args()
    ServerPart(case=args.case, icase=args.icase, prior=args.prior, uncertainty=args.unc,
               steps=args.steps, model_dir=args.model_dir, post_nc=args.post_nc,
               lc_site=args.lc_site, dry_run=args.dry_run, overwrite=args.overwrite,
               missing_2025=args.missing_2025, lc_mode=args.lc_mode)


if __name__ == '__main__':
    main()
