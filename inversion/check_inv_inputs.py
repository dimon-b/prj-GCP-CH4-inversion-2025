# -*- coding: utf-8 -*-
"""Check legacy F7/F5 NetCDF format before tdi_ch4 (NetCDF classic)."""
import argparse
from pathlib import Path

import numpy as np
from scipy.io import netcdf_file


# ============================================================
# === Read selected inversion station names (flag 9999)

def get_selected_sites(site_file):
    path = Path(site_file)
    selected = []
    for line in path.read_text(encoding='utf-8').splitlines()[2:]:
        fields = line.split()
        if len(fields) < 3:
            continue
        try:
            flag = int(fields[1])
        except ValueError:
            continue
        if flag == 9999:
            selected.append(fields[2])
    if not selected:
        raise ValueError('No active station rows flagged 9999 in ' + str(path))
    return selected


def get_name(row):
    return b''.join(row).decode('ascii', errors='replace').strip(' \x00')


# ============================================================
# === Inspect F7 monthly OBS-ACTM and uncertainties

def inspect_observations(fname, site_file=None, strict=False):
    fname = Path(fname)
    if not fname.is_file():
        raise FileNotFoundError(fname)

    with netcdf_file(str(fname), 'r', mmap=False) as nc:
        dims = nc.dimensions
        for key, value in [('station_count', 245), ('tlen', 324), ('data_type', 2), ('namelen', 45)]:
            if dims.get(key) != value:
                raise ValueError(str(fname) + ': wrong ' + key + ', expected ' + str(value) + ', got ' + str(dims.get(key)))
        required = ('tracer_data', 'site_name', 'station_count', 'tlen')
        for key in required:
            if key not in nc.variables:
                raise ValueError(str(fname) + ': missing variable ' + key)

        tracer_var = nc.variables['tracer_data']
        if tracer_var.dimensions != ('tlen', 'station_count', 'data_type'):
            raise ValueError(str(fname) + ': wrong tracer_data dimension order: ' + str(tracer_var.dimensions))
        arr = np.asarray(tracer_var[:], dtype=np.float64)
        time = np.asarray(nc.variables['tlen'][:], dtype=np.float64)
        names = [get_name(row) for row in nc.variables['site_name'][:]]
        stations = np.asarray(nc.variables['station_count'][:], dtype=int)

    if arr.shape != (324, 245, 2):
        raise ValueError(str(fname) + ': wrong tracer shape ' + str(arr.shape))
    if len(names) != 245 or not all(names):
        raise ValueError(str(fname) + ': missing/empty site names')
    # --- Legacy NetCDF contains repeated names; unique station IDs define distinct rows.
    #     read_globalview.f resolves the last matching station name; match that behavior.
    if len(set(stations)) != 245:
        raise ValueError(str(fname) + ': duplicate station numbers')
    if not np.all(np.diff(time) > 0):
        raise ValueError(str(fname) + ': time is not increasing')
    if not np.allclose(time[0], 1999.0 + 0.5 / 12, atol=1e-3):
        raise ValueError(str(fname) + ': wrong first time')
    if not np.allclose(time[-1], 2025.0 + 11.5 / 12, atol=1e-3):
        raise ValueError(str(fname) + ': wrong last time')

    differences = arr[:, :, 0]
    uncertainty = arr[:, :, 1]
    valid = np.isfinite(differences) & np.isfinite(uncertainty) & (differences != -999.) & (uncertainty > 0.) & (uncertainty != -999.)
    ninvalid_all = int(np.size(valid) - valid.sum())

    nselected = 0
    ninvalid_selected = 0
    selected_errors = []
    if site_file is not None:
        selected = get_selected_sites(site_file)
        index = {name: i for i, name in enumerate(names)}  # last matching index as in Fortran
        for name in selected:
            if name not in index:
                raise ValueError('Active site absent from F7 file: ' + name)
            k = index[name]
            bad = int((~valid[:, k]).sum())
            ninvalid_selected += bad
            if bad:
                selected_errors.append((name, bad))
        nselected = len(selected)
        if strict and ninvalid_selected:
            raise ValueError('F7 contains ' + str(ninvalid_selected) + ' missing/invalid selected-station months; '
                             'legacy read_globalview.f reads -999 without masking. First: ' + str(selected_errors[:5]))

    info = {'n_stations': 245, 'n_months': 324, 'n_active_stations': nselected,
            'n_invalid_all': ninvalid_all, 'n_invalid_active': ninvalid_selected,
            'active_missing_examples': selected_errors[:5], 'first_time': float(time[0]),
            'last_time': float(time[-1])}
    return info


# ============================================================
# === Inspect F5 source-prior NetCDF

def inspect_sources(fname):
    fname = Path(fname)
    if not fname.is_file():
        raise FileNotFoundError(fname)
    with netcdf_file(str(fname), 'r', mmap=False) as nc:
        dims = nc.dimensions
        for key, value in [('land_region', 54), ('month', 324), ('year', 27), ('presub_count', 4),
                           ('presubp1', 5), ('datatype_count', 2), ('namelen', 10)]:
            if dims.get(key) != value:
                raise ValueError(str(fname) + ': wrong ' + key + ', got ' + str(dims.get(key)))
        for key in ('sources_land', 'sources_presub', 'sources_offset', 'time', 'presub_name'):
            if key not in nc.variables:
                raise ValueError(str(fname) + ': missing variable ' + key)
        if nc.variables['sources_land'].dimensions != ('datatype_count', 'month', 'land_region'):
            raise ValueError(str(fname) + ': wrong sources_land dimension order')
        ts = np.asarray(nc.variables['time'][:], dtype=np.float64)
        land = np.asarray(nc.variables['sources_land'][:], dtype=np.float64)
    if land.shape != (2, 324, 54):
        raise ValueError(str(fname) + ': wrong sources_land shape')
    if not np.isfinite(land).all():
        raise ValueError(str(fname) + ': nonfinite source fluxes')
    if np.any(land[1] < 0):
        raise ValueError(str(fname) + ': negative source uncertainties')
    if not np.all(np.diff(ts) > 0):
        raise ValueError(str(fname) + ': non-monotonic prior time')
    if not np.isclose(ts[0], 1999 + 0.5 / 12, atol=1e-3) or not np.isclose(ts[-1], 2025 + 11.5 / 12, atol=1e-3):
        raise ValueError(str(fname) + ': wrong year interval')
    return {'n_regions': 54, 'n_months': 324, 'first_time': float(ts[0]), 'last_time': float(ts[-1])}


# ============================================================
# === CLI

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check F5/F7 classic NetCDF before inversion')
    parser.add_argument('--obs', required=True)
    parser.add_argument('--source', required=True)
    parser.add_argument('--sites', required=True)
    parser.add_argument('--strict', action='store_true')
    args = parser.parse_args()
    o = inspect_observations(args.obs, args.sites, strict=args.strict)
    s = inspect_sources(args.source)
    print('Observation F7:', o)
    print('Source F5     :', s)
    if o['n_invalid_active']:
        print('WARNING: selected station/months include missing data (-999). Review read_globalview.f before production.')
    print('Input format checks: PASS')
