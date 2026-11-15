#!/usr/bin/env python3
# ============================================================
# === Check the active CH4 station list and the zero-constraint file

from pathlib import Path
import argparse
import re
import numpy as np


# ============================================================
# === Read station list; preserve the order of active IDs

def get_sites(filename):
    lines = Path(filename).read_text(encoding='utf-8').splitlines()
    nsites = int(lines[1].split()[0])
    sites = []
    for line in lines[2:]:
        values = line.split()
        if not values:
            continue
        sites.append((int(values[0]), int(values[1]), values[2]))
    if len(sites) != nsites:
        raise ValueError('Wrong number of stations in datflist')
    return sites


# ============================================================
# === Count valid data as in read_globalview.f90

def get_valid_months(filename, sites):
    from scipy.io import netcdf_file
    with netcdf_file(str(filename), 'r', mmap=False) as f:
        raw_names = f.variables['site_name'][:].copy()
        vals = f.variables['tracer_data'][:].copy()
    names = [b''.join(row).decode('ascii', errors='replace').strip() for row in raw_names]
    count = 0
    for sid, flag, site_name in sites:
        if flag != 9999:
            continue
        matches = [i for i, name in enumerate(names) if name == site_name]
        if not matches:
            raise ValueError('Station not found in F7: ' + site_name)
        ix = matches[-1]
        obs = vals[:, ix, 0]
        sig = vals[:, ix, 1]
        valid = np.isfinite(obs) & np.isfinite(sig) & (obs != -999.) & (sig != -999.) & (sig > 0.)
        count += int(valid.sum())
    return count


# ============================================================
# === Main

def main():
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument('--sites', default=str(root / 'inv_dir/obsrvCH4_60sites.txt'))
    parser.add_argument('--control', default=str(Path(__file__).parent / 'control.example'))
    parser.add_argument('--f7', default=None)
    args = parser.parse_args()

    sites = get_sites(args.sites)
    active = [s for s in sites if s[1] == 9999]
    inactive = {s[0] for s in sites if s[1] == 0}
    expected_inactive = {31, 90, 91, 93}
    if not expected_inactive.issubset(inactive):
        raise AssertionError('Duplicate IDs 31,90,91,93 are still active')
    if len(active) != 57:
        raise AssertionError('Expected 57 active stations, got ' + str(len(active)))
    active_names = [s[2] for s in active]
    if len(active_names) != len(set(active_names)):
        raise AssertionError('Duplicate active names remain')

    control = Path(args.control).read_text(encoding='utf-8')
    match = re.search(r"(?im)^\s*constf\s*=\s*['\"]([^'\"]+)['\"]", control)
    if match is None:
        raise ValueError('No constf in control template')
    const_name = match.group(1)
    const_file = Path(const_name)
    if not const_file.is_absolute():
        const_file = Path(args.control).resolve().parent / const_file
    const_lines = const_file.read_text(encoding='utf-8').splitlines()
    if len(const_lines) < 2 or int(const_lines[1].split()[0]) != 0:
        raise AssertionError('Observation constraints are not zero')

    print('Active stations :', len(active))
    print('Disabled IDs    :', sorted(expected_inactive))
    print('Constraint file :', const_file)
    print('nconsts         :', 0)
    if args.f7:
        nvalid = get_valid_months(args.f7, sites)
        print('Valid F7 months :', nvalid)
    print('CLEAN NETWORK CHECK: PASS')


if __name__ == '__main__':
    main()
