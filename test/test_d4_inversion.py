# -*- coding: utf-8 -*-
"""Integration check for EXISTING TDI run outputs; never launches Fortran."""
import argparse
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from d4_inversion import Inversion


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=('INCA', 'CYC'), required=True)
    args = parser.parse_args()
    os.environ['GCP_CASE'] = args.case
    obj = Inversion(case=args.case)
    runs = list(obj.find_runs())
    if not runs:
        raise FileNotFoundError('No existing TDI output directories in ' + str(obj.result_root))
    bad = 0
    for path in runs:
        row = obj.diagnose_one(path)
        if row is None:
            continue
        print(row['status'], path)
        if row['status'] != 'PASS':
            print('  issues:', row['problems'])
            bad += 1
    if bad:
        raise RuntimeError('TDI integration check failed: ' + str(bad) + ' runs')
    print('D4 TDI integration: PASS')


if __name__ == '__main__':
    main()
