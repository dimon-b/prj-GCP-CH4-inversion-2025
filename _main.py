# -*- coding: utf-8 -*-
"""GCP/GMB 2025: 8 large workflow stages, one case-specific entry point."""
import argparse
import importlib
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STAGES = {
    'd0_OH': ('d0_inca_oh', 'INCA OH NetCDF + GT3'),
    'd1_prior': ('d1_prior_fluxes', 'Prior fluxes'),
    'd2_obs': ('d2_observations', 'Observations'),
    'd3_obs+actm': ('d3_obs_actm', 'OBS + ACTM'),
    'd4_inversion': ('d4_inversion', 'TDI inversion'),
    'd4_diagnostics': ('d4_inversion', 'TDI checks and diagnostics'),
    'd5_post': ('d5_post', 'Post-processing'),
    'd6_submition': ('d6_submit', 'Submission'),
    'd7_validation': ('d7_validation', 'Validation'),
}


def get_steps(text):
    names = [s.strip() for s in text.split(',') if s.strip()]
    if not names:
        raise ValueError('Empty --steps. Use: ' + ','.join(STAGES))
    bad = [s for s in names if s not in STAGES]
    if bad:
        raise ValueError('Unknown --steps ' + str(bad) + '; valid: ' + ','.join(STAGES))
    if len(set(names)) != len(names):
        raise ValueError('Repeated stage in --steps')
    return names


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=['INCA', 'CYC'], required=True)
    parser.add_argument('--steps', required=True, help='Comma-separated stage names')
    parser.add_argument('--dry-run', action='store_true', help='Show stage actions without running them')
    args = parser.parse_args(argv)
    names = get_steps(args.steps)
    if args.case == 'CYC' and 'd0_OH' in names:
        parser.error('d0_OH prepares INCA-OH only')
    old = os.environ.get('GCP_CASE')
    os.environ['GCP_CASE'] = args.case
    start = time.monotonic()
    try:
        print('Case:', args.case, '| Steps:', ', '.join(names), '| Dry-run:', args.dry_run)
        for name in names:
            module_name, title = STAGES[name]
            print('\n' + '=' * 62)
            print(name, title)
            print('=' * 62)
            module = importlib.import_module(module_name)
            kwargs = dict(case=args.case, dry_run=args.dry_run)
            if name == 'd4_diagnostics':
                module.run_diagnostics_stage(**kwargs)
            else:
                module.run_stage(**kwargs)
        print('\nALL SELECTED STAGES: PASS')
        print('Elapsed seconds:', round(time.monotonic()-start, 1))
    finally:
        if old is None:
            os.environ.pop('GCP_CASE', None)
        else:
            os.environ['GCP_CASE'] = old


if __name__ == '__main__':
    main()
