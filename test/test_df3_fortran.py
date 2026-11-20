# -*- coding: utf-8 -*-
"""Read-only D3 production results check: python test/test_df3_fortran.py --case INCA

Not a numerical regression of Nakazawa; requires original routines and real data.
"""
import argparse
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from d3_obs_actm import ObsACTM


def check_file(name, path):
    p = Path(path).expanduser().resolve()
    if not p.is_file() or p.stat().st_size == 0:
        raise FileNotFoundError(name + ': ' + str(p))
    print('PASS:', name, p.name, p.stat().st_size, 'bytes')


def main(case, build=False):
    os.environ['GCP_CASE'] = case
    os.environ['CH4_CASE'] = case
    obj = ObsACTM(case)
    if build:
        for op in ('sample', 'combine', 'filter', 'observations'):
            print('D3 build:', obj.executable(op))
        print('D3 Fortran build PASS')
        return
    root = Path(obj.d3_sample_output_dir)
    for year in range(obj.yr_s, obj.yr_e+1):
        check_file('Sample ' + str(year), root / (obj.d3_sample_output_prefix + str(year) + '.grd'))
    for suffix in ('ch4yr.grd', 'ch4mm.grd', 'ch4dy.grd'):
        check_file('Combine', Path(obj.d3_combine_output_dir) / (obj.d3_combine_output_prefix + suffix))
    check_file('Filter', obj.d3_filter_output_file)
    check_file('Observations summary', obj.d3_observations_summary)
    check_file('Observations prior', obj.d3_observations_prior)
    for unc in obj.unxcases:
        if unc != 'ctl':
            check_file('Observations summary ' + unc, obj.d3_observations_summary + '.' + unc)
            check_file('Observations prior ' + unc, obj.d3_observations_prior + '.' + unc)
    for unc in obj.unxcases:
        for inv in obj.invcases:
            check_file('Observations NetCDF', str(obj.d3_observations_nc_prefix) + unc + inv + '.nc')
    print('D3 output inventory PASS (file presence/size only)')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=['CYC', 'INCA'], default='INCA')
    parser.add_argument('--build', action='store_true', help='Compile all D3 Fortran sources')
    args = parser.parse_args()
    main(args.case, build=args.build)
