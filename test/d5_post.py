# -*- coding: utf-8 -*-
"""D5: assemble posterior CH4 fluxes from D4 TDI corrections.

Run from c_gcpv3:
    python d5_post.py --case CYC --prior p30 --unc ctl
    python d5_post.py --case INCA --prior p30 --unc ctl
    python _main.py --case CYC --steps d5_post  # after registration

D4 and D1 are READ ONLY.  A single Fortran code writes 1-degree GRD,
conservative 128x64 GT3, and machine-readable monthly / regional budgets.
Loss correction is NOT calculated or applied in this version.
Diagnostics are ONLY under results2025/<CASE>/diagnostics/d5/.
"""
import argparse
import calendar
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from _set_case import SetCase


# ============================================================
# === D5: POST fluxes
class PostFlux(SetCase):

    def __init__(self, case=None):
        super().__init__()
        if case is not None and self.hcase != str(case).upper():
            raise ValueError('Case mismatch: ' + self.hcase + ' vs ' + str(case))
        self.case = self.hcase
        self.root = Path(self.inv_wrk_dir).resolve()
        self.diag_root = self.root / 'diagnostics' / 'd5'
        self.source_code = Path(__file__).resolve().parent / 'fortran' / 'df5_post_flux.f90'
        self.executable = Path(os.environ.get('CH4_EXE_POST',
                                str(Path(__file__).resolve().parent / 'fortran' / 'bin' / 'df5_post_flux'))).resolve()

    # --------------------------------------------------------
    # --- Resolve D1/D4 inputs from current SetCase; no stage_impl
    def get_paths(self, prior='p30', unc='ctl', model='inv1'):
        if prior not in self.unpcases or unc not in self.unxcases or model not in self.invcases:
            raise ValueError('Invalid experiment: ' + str((prior, unc, model)))
        run = self.root / (self.icase + '_' + prior) / (self.icase + '_' + unc + '_' + model)
        name = self.icase + '_' + unc + '_' + model
        suffix = '_noLC'
        output = self.root / 'flux2d' / prior / (name + suffix)
        diag = self.diag_root / (prior + '_' + unc + '_' + model + suffix)
        return {
            'run': run,
            'source': run / 'predicted.sources',
            'prior': Path(self.flx_apr_grd).resolve(),
            'mask': Path(self.inv_dir).resolve() / 'mask' / 'mask80.bin',
            'basis': Path(self.inv_dir).resolve() / 'mask' / 'ch4_regn_pulse.grd',
            'grd': Path(str(output) + '.grd'),
            'gt3': Path(str(output) + '.gt3'),
            'diag': diag,
            'budget': diag / 'monthly_budget.csv',
            'regional': diag / 'regional_budget.csv',
        }

    @staticmethod
    def check_file(path):
        path = Path(path)
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError('Missing or empty D5 input: ' + str(path))
        if path.is_symlink():
            raise ValueError('D5 does not accept symlinked input: ' + str(path))

    # --------------------------------------------------------
    # --- 1-degree GRD has no date metadata: verify record count explicitly
    def prior_start_record(self, filename):
        rec_bytes = 360 * 180 * 4
        size = Path(filename).stat().st_size
        if size % rec_bytes:
            raise ValueError('Prior GRD byte size is not a multiple of one 360x180 record: ' + str(filename))
        records = size // rec_bytes
        if (self.yr_s, self.yr_e, self.ntime) != (1999, 2025, 324):
            raise ValueError('Unexpected configured inversion period')
        if records == 324:
            return 1
        if records == 348:
            # D1 years 1998-2026, skip 1998 (12 records)
            return 13
        raise ValueError('Unsupported prior GRD record count: ' + str(records) + ' (expected 324 or 348)')

    def check_inputs(self, paths):
        for name in ('source', 'prior', 'mask', 'basis'):
            self.check_file(paths[name])
        rec_bytes = 360 * 180 * 4
        for name in ('mask', 'basis'):
            size = paths[name].stat().st_size
            if size != rec_bytes:
                raise ValueError(name + ' expected one 360x180 real32 record, got ' + str(size))
        return self.prior_start_record(paths['prior'])

    # --------------------------------------------------------
    # --- Compile only this Fortran source; no changes to D4 or legacy files
    def build_fortran(self):
        if not self.source_code.is_file():
            raise FileNotFoundError('Missing D5 Fortran source: ' + str(self.source_code))
        override = os.environ.get('CH4_EXE_POST', '').strip()
        if override:
            self.check_file(self.executable)
            return
        if self.executable.is_file() and self.executable.stat().st_mtime >= self.source_code.stat().st_mtime:
            return
        compiler = os.environ.get('FC', '').strip()
        if not compiler:
            compiler = 'gfortran' if shutil.which('gfortran') else 'ifort'
        if not shutil.which(compiler):
            raise RuntimeError('Fortran compiler unavailable: ' + compiler + '. Set FC or CH4_EXE_POST')
        flags = ['-O2', '-ffree-line-length-none'] if 'gfortran' in Path(compiler).name else ['-O2', '-free']
        self.executable.parent.mkdir(parents=True, exist_ok=True)
        cmd = [compiler] + flags + ['-o', str(self.executable), str(self.source_code)]
        print('D5 build:', ' '.join(cmd))
        subprocess.run(cmd, check=True, cwd=str(self.executable.parent))

    @staticmethod
    def namelist_quote(value):
        text = str(value)
        if '\n' in text or '\r' in text or "'" in text or len(text) > 1024:
            raise ValueError('Fortran namelist path invalid: ' + text)
        return "'" + text + "'"

    def make_namelist(self, paths, prior_start_record, staged):
        vals = {
            'source_file': paths['source'], 'prior_file': paths['prior'],
            'mask_file': paths['mask'], 'basis_file': paths['basis'],
            'out_grd': staged['grd'], 'out_gt3': staged['gt3'],
            'budget_file': staged['budget'], 'region_file': staged['regional'],
        }
        lines = ['&post_config']
        for key, val in vals.items():
            lines.append(' ' + key + ' = ' + self.namelist_quote(val))
        lines += [' initial_year = ' + str(self.yr_s),
                  ' nmonths = ' + str(self.ntime),
                  ' prior_start_record = ' + str(prior_start_record),
                  ' write_gt3 = .true.', '/']
        staged['config'].write_text('\n'.join(lines) + '\n', encoding='ascii')

    # --------------------------------------------------------
    # --- Save one diagnostic figure (no D3/F5 plots duplicated)
    def diagnostics(self, paths):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        df = pd.read_csv(paths['budget'])
        region = pd.read_csv(paths['regional'])
        if len(df) != self.ntime or len(region) != self.ntime*54:
            raise RuntimeError('Incomplete D5 budget output')
        for col in ('prior_Tg_month', 'correction_Tg_month',
                    'post_Tg_month', 'gt3_Tg_month'):
            if not np.isfinite(df[col]).all():
                raise RuntimeError('Invalid mass-budget value: ' + col)
        if not np.allclose(df['post_Tg_month'], df['prior_Tg_month'] + df['correction_Tg_month'],
                           rtol=3e-7, atol=1e-7):
            raise RuntimeError('D5 posterior mass-budget inconsistency')
        if not np.allclose(df['post_Tg_month'], df['gt3_Tg_month'], rtol=3e-7, atol=1e-7):
            raise RuntimeError('D5 GRD/GT3 mass-budget mismatch')
        if not np.allclose(region['tdi_Tg_per_365d_year'],region['gridded_Tg_per_365d_year'],
                           rtol=3e-7, atol=2e-8):
            raise RuntimeError('D5 regional correction budget mismatch')

        annual = df.groupby('year',as_index=False)[['prior_Tg_month','correction_Tg_month',
                       'post_Tg_month','gt3_Tg_month']].sum()
        annual.to_csv(paths['diag'] / 'annual_budget.csv',index=False)
        core = region[(region.year>=2000)&(region.year<=2024)]
        core_stats=core.groupby('region',as_index=False)['tdi_Tg_per_365d_year'].agg(['mean', 'std'])
        core_stats.to_csv(paths['diag'] / 'regional_2000_2024.csv')
        top=core_stats.assign(absm=lambda x:x['mean'].abs()).nlargest(12,'absm')

        fig,axes=plt.subplots(2,2,figsize=(14,9),constrained_layout=True)
        axes[0,0].plot(annual.year,annual.prior_Tg_month,label='Prior',lw=1.5)
        axes[0,0].plot(annual.year,annual.post_Tg_month,label='Posterior',lw=1.5)
        axes[0,0].set_ylabel('CH4 Tg / calendar year')
        axes[0,0].set_title('Global mass budget')
        axes[0,0].legend()
        axes[0,1].plot(annual.year,annual.correction_Tg_month,label='Regional correction')
        axes[0,1].axhline(0,color='gray',lw=.7)
        axes[0,1].set_ylabel('CH4 Tg / calendar year')
        axes[0,1].set_title('Regional source correction (no LossCorr)')
        axes[0,1].legend()
        timeline=pd.to_datetime(dict(year=df.year,month=df.month,day=1))
        axes[1,0].plot(timeline,df.post_Tg_month,label='Posterior')
        axes[1,0].plot(timeline,df.prior_Tg_month,alpha=.6,label='Prior')
        axes[1,0].set_ylabel('CH4 Tg / month')
        axes[1,0].set_title('Monthly flux (actual month lengths)')
        axes[1,0].legend()
        axes[1,1].barh(top.region.astype(str),top['mean'])
        axes[1,1].axvline(0,color='gray',lw=.7)
        axes[1,1].set_xlabel('Mean TDI regional correction (Tg / 365-day year)')
        axes[1,1].set_title('Largest signed regional means, 2000–2024')
        for ax in axes.flat:
            ax.grid(alpha=.2)
        fig.suptitle('D5 POST: '+self.case+' / no LossCorr / '+self.icase)
        fig.savefig(paths['diag'] / 'post_overview.png',dpi=160)
        plt.close(fig)

        maxerr=float(np.max(np.abs(df.gt3_rel_error)))
        neg=int(df.negative_cells.max())
        lines=['D5 POST / '+self.case+' / '+self.icase,
               'Loss correction: NOT APPLIED (prior + D4 regional corrections only)',
               'Period: '+str(self.yr_s)+'–'+str(self.yr_e),
               'Prior GRD: '+str(paths['prior']),
               'Predicted regional sources: '+str(paths['source']),
               'Global prior (2000–2024 mean annual): '+str(round(annual[(annual.year>=2000)&(annual.year<=2024)].prior_Tg_month.mean(),4))+' Tg/year',
               'Global posterior (2000–2024 mean annual): '+str(round(annual[(annual.year>=2000)&(annual.year<=2024)].post_Tg_month.mean(),4))+' Tg/year',
               'Max relative GRD→GT3 conservation error: '+str(maxerr),
               'Max monthly negative grid-cell count: '+str(neg),
               'Regional and global budgets: PASS',
               '1999 inversion boundary excluded from core means',
               '2025 months retained for ACTM continuity; missing observations are not filled.',
               'GRD and GT3 = kg CH4 m^-2 s^-1; budgets = Tg CH4/month or Tg/year.',
               'TDI native coefficients = Tg CH4 per 365-day year; not Tg/month.']
        (paths['diag'] / 'diagnostic_report.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
        print('\n'.join(lines))
        return df

    # --------------------------------------------------------
    # --- Never silently overwrite existing outputs
    def run(self, prior='p30', unc='ctl', model='inv1', dry_run=False, overwrite=False):
        paths = self.get_paths(prior, unc, model)
        start = self.check_inputs(paths)
        print('D5 post:',self.case,self.icase,prior,unc,model)
        print('  D4 corrections:',paths['source'])
        print('  Prior:',paths['prior'],'| start record:',start)
        print('  Output:',paths['grd'],'and',paths['gt3'])
        print('  Loss correction: DISABLED (not calculated or applied)')
        if dry_run:
            print('D5 DRY-RUN: no files written')
            return
        if (paths['grd'].exists() or paths['gt3'].exists()) and not overwrite:
            raise FileExistsError('Posterior output already exists; use --overwrite explicitly')
        self.build_fortran()
        paths['grd'].parent.mkdir(parents=True, exist_ok=True)
        paths['diag'].mkdir(parents=True, exist_ok=True)
        work = paths['diag'] / '_work'
        work.mkdir(exist_ok=True)
        staged = {
            'grd':work/'posterior.grd', 'gt3':work/'posterior.gt3',
            'budget':work/'monthly_budget.csv', 'regional':work/'regional_budget.csv',
            'config':work/'post_config.nml'
        }
        for key in ('grd','gt3','budget','regional'):
            if staged[key].exists():
                raise FileExistsError('Previous interrupted D5 output: ' + str(staged[key]))
        self.make_namelist(paths,start,staged)
        logfile=paths['diag']/'fortran.log'
        with logfile.open('w',encoding='utf-8') as log:
            subprocess.run([str(self.executable),str(staged['config'])],cwd=str(work),
                           stdout=log,stderr=subprocess.STDOUT,check=True)
        # --- Expected physical output size, no GRD record markers
        if staged['grd'].stat().st_size != self.ntime*360*180*4:
            raise RuntimeError('Unexpected posterior GRD byte size')
        if staged['gt3'].stat().st_size != self.ntime*(1024+8+128*64*4+8):
            raise RuntimeError('Unexpected GT3 byte size')
        # --- Move diagnostic CSV first and verify before publication of model flux
        staged['budget'].replace(paths['budget'])
        staged['regional'].replace(paths['regional'])
        self.diagnostics(paths)
        staged['grd'].replace(paths['grd'])
        staged['gt3'].replace(paths['gt3'])
        print('D5 POST: PASS',paths['grd'],paths['gt3'])


# ============================================================
# === CLI / _main.py integration

def run_stage(case='INCA', dry_run=False):
    old=os.environ.get('GCP_CASE')
    os.environ['GCP_CASE']=str(case).upper()
    try:
        PostFlux(case=case).run(dry_run=dry_run)
    finally:
        if old is None:
            os.environ.pop('GCP_CASE',None)
        else:
            os.environ['GCP_CASE']=old


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',required=True,choices=['CYC','INCA'])
    parser.add_argument('--prior',default='p30',choices=['p30','p50','p99'])
    parser.add_argument('--unc',default='ctl',choices=['ctl','ux2','ux4'])
    parser.add_argument('--model',default='inv1')
    parser.add_argument('--dry-run',action='store_true')
    parser.add_argument('--overwrite',action='store_true')
    args=parser.parse_args()
    os.environ['GCP_CASE']=args.case
    PostFlux(case=args.case).run(prior=args.prior,unc=args.unc,model=args.model,
                                dry_run=args.dry_run,overwrite=args.overwrite)
