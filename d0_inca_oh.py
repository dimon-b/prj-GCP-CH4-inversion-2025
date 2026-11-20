# -*- coding: utf-8 -*-
"""
INCA OH -> MIROC4-ACTM T42L67 (nominal, whole-period production candidate).

Uses original annual daily files. Daily OH is horizontally interpolated with
periodic longitude, then vertically interpolated in log(physical pressure),
using actual pmid(time, level, lat, lon) and ps(time, lat, lon) if available.
An explicit --vertical-mode nominal option permits a REFERENCE-PRESSURE,
NOT PHYSICAL-PRESSURE, sensitivity test from presnivs alone. No silent fallback.
Monthly means are computed AFTER each day's remapping.

Original tch4_INCA_OH_t42.nc is never modified.
The nominal mode remains APPROXIMATE: real pmid/Ps are not provided in OH-only files.
"""

import argparse
import calendar
import math
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


# =============================================================
# === Names / units
# =============================================================

def get_paths():
    project = Path(__file__).resolve().parent.parent
    inp = project / 'fluxes/gcp2025_flux_inp/LOSS_FIELDS/OH_INCA/Scaled'
    ref = project / 'fluxes/gcp2025_OH/tch4_actm_gcp21oh.nc'
    out = project / 'fluxes/gcp2025_OH/tch4_INCA_OH_t42_v2.nc'
    return inp, ref, out


def get_name(ds, candidates, required=True):
    names = {name.lower(): name for name in ds.variables}
    for name in candidates:
        if name.lower() in names:
            return names[name.lower()]
    if required:
        raise ValueError(f'Missing variable {candidates}. Available: {list(ds.variables)}')
    return None


def get_dim(ds, candidates):
    names = {name.lower(): name for name in ds.dims}
    for name in candidates:
        if name.lower() in names:
            return names[name.lower()]
    raise ValueError(f'Missing dimension {candidates}. Available: {list(ds.dims)}')


def get_oh_factor(units, assumed=None):
    # --- Parse the concentration unit AFTER the molecular-count prefix.
    # 'molec. m-3' must NOT be misclassified as 'cm-3' because 'molec'
    # itself ends with 'c'. This caused a potential 1e6 conversion error.
    raw = str(units or assumed or '').lower().strip()
    raw = raw.replace('molecules', 'molec').replace('molecule', 'molec')
    if not raw.startswith('molec'):
        raise ValueError(f'OH units not confirmed: {units!r}')
    dimensional = raw[len('molec'):].strip(' ./\t').replace(' ', '')
    dimensional = dimensional.replace('**', '').replace('^', '')
    if dimensional in ('cm-3', 'cm3-1', 'cm3'):
        return 1.0
    if dimensional in ('m-3', 'm3-1', 'm3'):
        return 1e-6
    raise ValueError(f'Unknown OH concentration units: {units!r}, parsed dimension {dimensional!r}')


def get_pressure_factor(units, assumed=None):
    units = (units or assumed or '').lower().strip().replace(' ', '')
    if units in ('pa', 'pascal', 'pascals'):
        return 1.0
    if units in ('hpa', 'mbar', 'millibar', 'millibars'):
        return 100.0
    raise ValueError(f'Pressure units not confirmed: {units!r}. Use --assume-pressure-unit only if justified.')


def get_source_metadata(ds, args):
    tdim = get_dim(ds, ['time_counter', 'time'])
    zdim = get_dim(ds, ['presnivs', 'level', 'lev', 'z'])
    ydim = get_dim(ds, ['lat', 'latitude', 'y'])
    xdim = get_dim(ds, ['lon', 'longitude', 'x'])
    oh_name = get_name(ds, ['OH'])
    mode = getattr(args, 'vertical_mode', 'physical')
    oh_factor = get_oh_factor(ds[oh_name].attrs.get('units'), getattr(args, 'assume_oh_unit', None))

    expected_oh = {tdim, zdim, ydim, xdim}
    if set(ds[oh_name].dims) != expected_oh:
        raise ValueError(f'OH dims {ds[oh_name].dims}; expected {expected_oh}')

    cfg = dict(t=tdim, z=zdim, y=ydim, x=xdim, oh=oh_name,
               oh_factor=oh_factor, mode=mode)

    if mode == 'physical':
        p_name = get_name(ds, ['pmid', 'p_mid', 'pressure_mid', 'pressure'], required=False)
        ps_name = getattr(args, 'ps_var', None) or get_name(
            ds, ['ps', 'psol', 'sp', 'psurf', 'surface_pressure', 'surfacepressure'],
            required=False)
        if p_name is None or ps_name is None or ps_name not in ds:
            raise ValueError(
                'Physical-pressure mapping requires pmid(time, level, lat, lon) '
                'AND ps(time, lat, lon). They are absent in this OH-only file. '
                'LMDZ presnivs is PSEUDO-pressure, not actual column pressure. '
                'Find companion meteorological fields, or explicitly select '
                '--vertical-mode nominal for an approximate sensitivity experiment.')
        if set(ds[p_name].dims) != expected_oh:
            raise ValueError(f'pmid dims {ds[p_name].dims}; expected {expected_oh}')
        if set(ds[ps_name].dims) != {tdim, ydim, xdim}:
            raise ValueError(f'ps dims {ds[ps_name].dims}; expected time/lat/lon')
        cfg.update(p=p_name, ps=ps_name,
                   p_factor=get_pressure_factor(ds[p_name].attrs.get('units'),
                                                getattr(args, 'assume_pressure_unit', None)),
                   ps_factor=get_pressure_factor(ds[ps_name].attrs.get('units'),
                                                 getattr(args, 'assume_pressure_unit', None)))

    elif mode == 'nominal':
        # presnivs is pseudo-pressure at reference surface pressure in LMDZ.
        # It is NOT actual pmid(t, z, y, x). This is intentionally approximate.
        if zdim not in ds.coords or ds[zdim].dims != (zdim,):
            raise ValueError('Nominal mode requires one-dimensional presnivs coordinate')
        factor = get_pressure_factor(ds[zdim].attrs.get('units'),
                                     getattr(args, 'assume_pressure_unit', None))
        pref = float(getattr(args, 'nominal_reference_pressure_pa', 101325.0))
        if not 80000.0 <= pref <= 110000.0:
            raise ValueError('Invalid nominal reference pressure; expected 80000..110000 Pa')
        pressure = np.asarray(ds[zdim].values, dtype=np.float64) * factor
        if not np.isfinite(pressure).all() or (pressure <= 0).any():
            raise ValueError('Invalid nominal presnivs pressure coordinate')
        if len(np.unique(pressure)) != len(pressure):
            raise ValueError('Duplicate source pressure levels')
        if pressure.max() > 120000.0:
            raise ValueError('Source pseudo-pressure is outside atmospheric range')
        if pressure.max() < 50000.0:
            raise ValueError('Source pseudo-pressure does not reach the lower troposphere')
        cfg.update(presnivs_pa=pressure, nominal_ref_pa=pref)

    else:
        raise ValueError(f'Unknown vertical mode: {mode}')
    return cfg


def inspect_daily(path, args):
    print(f'INPUT: {path}')
    with xr.open_dataset(path, decode_times=True) as ds:
        print('dims:', dict(ds.sizes))
        for name in ds.data_vars:
            var = ds[name]
            print(f'{name:24s} dims={var.dims} units={var.attrs.get("units", "UNSPECIFIED")}')
        print('coords:', list(ds.coords))
        zdim = get_dim(ds, ['presnivs', 'level', 'lev', 'z'])
        if zdim in ds.coords:
            z = ds[zdim]
            values = z.to_numpy().astype(float)
            print(f'{zdim} attributes: {dict(z.attrs)}')
            print(f'{zdim} range: {np.nanmin(values):.9g} .. {np.nanmax(values):.9g}')
            print(f'{zdim} first levels:', np.array2string(values[:8], precision=7))
            print(f'{zdim} last levels:', np.array2string(values[-8:], precision=7))
            print('NOTE: LMDZ presnivs is reference pseudo-pressure, not column-resolved pmid.')
        try:
            cfg = get_source_metadata(ds, args)
        except (ValueError, KeyError) as e:
            print('INPUT CHECK: BLOCKED:', e)
            print('NO output file created.')
            return False
        print('INPUT CHECK: PASS, vertical mode =', cfg['mode'])
        print('OH factor to molecules cm-3:', cfg['oh_factor'])
        if cfg['mode'] == 'physical':
            print('physical pmid:', cfg['p'], 'factor to Pa:', cfg['p_factor'])
            print('physical ps:', cfg['ps'], 'factor to Pa:', cfg['ps_factor'])
            block = ds[[cfg['oh'], cfg['p'], cfg['ps']]].isel({cfg['t']: 0}).load()
            for name in [cfg['oh'], cfg['p'], cfg['ps']]:
                a = block[name].to_numpy()
                finite = a[np.isfinite(a)]
                if finite.size:
                    print(f'{name}: finite={finite.size}/{a.size} min={finite.min():.6g} max={finite.max():.6g}')
        else:
            print('WARNING: nominal pseudo-pressure mapping is an APPROXIMATION.')
            print('reference pressure, Pa:', cfg['nominal_ref_pa'])
            print('DO NOT use this as physically corrected INCA-OH until compared with')
            print('column-resolved pressure fields and forward CH4 sensitivities.')
        return True


# =============================================================
# === Regrid one source timestep
# =============================================================

def horizontal_regrid(block, cfg, lats, lons):
    names = {cfg['z']: 'z', cfg['y']: 'y', cfg['x']: 'x'}
    block = block.rename(names)
    x = np.mod(block['x'].values.astype(float), 360.0)
    block = block.assign_coords(x=x).sortby('x').sortby('y')

    # --- De-duplicate an existing 360-degree endpoint
    u, first = np.unique(np.round(block['x'].values, 7), return_index=True)
    block = block.isel(x=np.sort(first))
    x = block['x'].values
    if len(x) < 2:
        raise ValueError('Need at least two distinct source longitudes')

    # --- Periodic longitude seam; NO extrapolated 0/360 gap
    left = block.isel(x=[-1]).assign_coords(x=[float(x[-1]) - 360.0])
    right = block.isel(x=[0]).assign_coords(x=[float(x[0]) + 360.0])
    ext = xr.concat([left, block, right], dim='x')
    result = ext.interp(y=lats, x=np.mod(lons, 360.0), method='linear')
    return result


def interp_vertical(oh, pmid_pa, ps_pa, sigma, max_outside_pct=10.0):
    # --- Arrays (source_level, lat, lon) and (lat, lon)
    if oh.ndim != 3 or pmid_pa.shape != oh.shape or ps_pa.shape != oh.shape[1:]:
        raise ValueError(f'Invalid source shapes: OH={oh.shape}, pmid={pmid_pa.shape}, ps={ps_pa.shape}')

    nz, ny, nx = oh.shape
    ncol = ny * nx
    p = pmid_pa.reshape(nz, ncol).astype(float)
    v = oh.reshape(nz, ncol).astype(float)
    ps = ps_pa.reshape(ncol).astype(float)

    if not np.isfinite(ps).all() or (ps < 30000).any() or (ps > 115000).any():
        raise ValueError('Invalid surface pressure. Expected physical Ps in Pa, roughly 30–115 kPa.')

    valid = np.isfinite(p) & (p > 0) & np.isfinite(v) & (v >= 0)
    nvalid = valid.sum(axis=0)
    if np.any(nvalid < 2):
        raise ValueError(f'Insufficient valid OH/pmid levels in {(nvalid < 2).sum()} columns')

    # --- Sort independently within each column (not by nominal presnivs)
    idx = np.argsort(np.where(valid, p, np.inf), axis=0)
    p = np.take_along_axis(p, idx, axis=0)
    v = np.take_along_axis(v, idx, axis=0)

    target = np.asarray(sigma, dtype=float)[:, None] * ps[None, :]
    if not (np.isfinite(target) & (target > 0)).all():
        raise ValueError('Invalid ACTM target pressures')

    col = np.arange(ncol)[None, :]
    p_min = p[0, :][None, :]
    p_max = p[nvalid - 1, np.arange(ncol)][None, :]
    outside = (target < p_min) | (target > p_max)
    outside_pct = 100.0 * outside.sum() / outside.size
    if outside_pct > max_outside_pct:
        raise ValueError(f'{outside_pct:.2f}% target levels outside source pressure range; '
                         f'max allowed {max_outside_pct:.2f}%. Inspect vertical coverage.')

    # --- Vectorized lower_bound over pressure for all target levels/columns
    lo = np.zeros(target.shape, dtype=np.int32)
    hi = np.broadcast_to((nvalid - 1)[None, :], target.shape).copy().astype(np.int32)
    for _ in range(int(math.ceil(math.log2(nz))) + 2):
        mid = (lo + hi) // 2
        go_right = p[mid, col] < target
        lo = np.where(go_right, np.minimum(mid + 1, nvalid[None, :] - 1), lo)
        hi = np.where(go_right, hi, mid)

    upper = lo
    lower = np.maximum(upper - 1, 0)
    p0 = p[lower, col]
    p1 = p[upper, col]
    v0 = v[lower, col]
    v1 = v[upper, col]

    # --- Constant edge continuation beyond coverage, always reported above
    denom = np.log(p1) - np.log(p0)
    frac = np.divide(np.log(target) - np.log(p0), denom,
                     out=np.zeros_like(target), where=(denom > 1e-12))
    frac = np.clip(frac, 0.0, 1.0)
    result = v0 + (v1 - v0) * frac
    if not np.isfinite(result).all() or (result < 0).any():
        raise ValueError('Vertical interpolation produced NaN, infinite or negative OH')
    return result.reshape((len(sigma), ny, nx)).astype(np.float32), outside_pct


def convert_one_day(ds, cfg, idx, lats, lons, sigma, max_outside_pct):
    fields = [cfg['oh']]
    if cfg['mode'] == 'physical':
        fields += [cfg['p'], cfg['ps']]
    block = ds[fields].isel({cfg['t']: idx}, drop=True).load()
    block = horizontal_regrid(block, cfg, lats, lons)
    oh = block[cfg['oh']].transpose('z', 'y', 'x').to_numpy() * cfg['oh_factor']
    if cfg['mode'] == 'physical':
        pmid = block[cfg['p']].transpose('z', 'y', 'x').to_numpy() * cfg['p_factor']
        ps = block[cfg['ps']].transpose('y', 'x').to_numpy() * cfg['ps_factor']
    else:
        ny, nx = len(lats), len(lons)
        pmid = np.broadcast_to(cfg['presnivs_pa'][:, None, None], oh.shape)
        ps = np.full((ny, nx), cfg['nominal_ref_pa'], dtype=np.float64)
    return interp_vertical(oh, pmid, ps, sigma, max_outside_pct=max_outside_pct)


def convert_one_month_nominal(ds, cfg, idxs, lats, lons, sigma, max_outside_pct):
    # --- With fixed presnivs and fixed reference pressure, remapping weights
    # --- are time-independent: interpolate(monthly mean) == mean(interpolate(daily)).
    # --- This avoids thousands of costly daily interpolation operations.
    mean_source = ds[[cfg['oh']]].isel({cfg['t']: idxs}).mean(cfg['t'], skipna=False)
    block = horizontal_regrid(mean_source, cfg, lats, lons)
    oh = block[cfg['oh']].transpose('z', 'y', 'x').to_numpy() * cfg['oh_factor']
    pmid = np.broadcast_to(cfg['presnivs_pa'][:, None, None], oh.shape)
    ps = np.full((len(lats), len(lons)), cfg['nominal_ref_pa'], dtype=np.float64)
    return interp_vertical(oh, pmid, ps, sigma, max_outside_pct=max_outside_pct)


# =============================================================
# === Output
# =============================================================

class NCWriter:
    def __init__(self, path, lat, lon, sigma, vertical_mode, nominal_ref_pa=None):
        try:
            from netCDF4 import Dataset
        except ImportError:
            raise RuntimeError('netCDF4 Python module is required for the production output')

        self.ds = Dataset(path, 'w', format='NETCDF4')
        self.ds.createDimension('time', None)
        self.ds.createDimension('level', len(sigma))
        self.ds.createDimension('y', len(lat))
        self.ds.createDimension('x', len(lon))

        self.t = self.ds.createVariable('time', 'f8', ('time',))
        self.t.units = 'days since 1998-01-01 00:00:00'
        self.t.calendar = 'proleptic_gregorian'
        self.z = self.ds.createVariable('level', 'f4', ('level',))
        self.y = self.ds.createVariable('y', 'f4', ('y',))
        self.x = self.ds.createVariable('x', 'f4', ('x',))
        self.z[:] = sigma
        self.y[:] = lat
        self.x[:] = lon
        self.z.units = '1'
        self.z.long_name = 'ACTM sigma level'
        self.y.units = 'degrees_north'
        self.x.units = 'degrees_east'
        self.oh = self.ds.createVariable('OH', 'f4', ('time', 'level', 'y', 'x'),
                                         zlib=True, complevel=2, shuffle=True,
                                         chunksizes=(1, len(sigma), len(lat), len(lon)),
                                         fill_value=np.float32(-9999.0))
        self.oh.units = 'molecules cm-3'
        self.oh.long_name = 'INCA hydroxyl radical, daily-remapped monthly mean'
        self.ds.vertical_mapping = vertical_mode
        if vertical_mode == 'nominal':
            self.ds.nominal_reference_pressure_pa = float(nominal_ref_pa)
            self.ds.IMPORTANT = ('APPROXIMATE LMDZ pseudo-pressure mapping; '
                                 'NOT valid column-resolved physical pressure. '
                                 'For scientific production use obtain pmid and ps.')
        self.ds.history = ('INCA OH v2: daily-remapped monthly means in physical mode; '
                           'monthly-first equivalent in nominal mode; periodic longitude; '
                           'log-pressure vertical interpolation; no extra OH scaling. '
                           'Vertical mode: ' + vertical_mode)
        self.i = 0

    def write(self, year, month, oh):
        date = (pd.Timestamp(year=year, month=month, day=16) - pd.Timestamp('1998-01-01')).days
        self.t[self.i] = float(date)
        self.oh[self.i, :, :, :] = oh
        self.i += 1
        self.ds.sync()

    def close(self):
        self.ds.close()


def build_months(ds, cfg, year, lats, lons, sigma, args):
    tdim = cfg['t']
    time = ds[tdim]
    dates = pd.DataFrame({
        'year': time.dt.year.to_numpy().astype(int),
        'month': time.dt.month.to_numpy().astype(int),
        'day': time.dt.day.to_numpy().astype(int),
    })
    if not np.all(dates.year == year):
        raise ValueError(f'{year}: source dates include other years')

    # --- In production, all Gregorian daily samples must be present, once each.
    # --- A 20-day threshold alone can silently average incomplete months.
    strict_daily = not getattr(args, 'allow_incomplete_months', False)
    keys = dates[['month', 'day']]
    if keys.duplicated().any():
        duplicate = keys.loc[keys.duplicated(), :].head(5)
        raise ValueError(f'{year}: duplicate daily dates: {duplicate.to_dict("records")}')
    if strict_daily:
        expected_count = 366 if calendar.isleap(year) else 365
        if len(dates) != expected_count:
            raise ValueError(f'{year}: expected {expected_count} daily records, got {len(dates)}')

    monthly = []
    for month in range(1, 13):
        idxs = np.flatnonzero(dates.month.to_numpy() == month)
        ndays = dates.iloc[idxs].day.nunique()
        expected_days = calendar.monthrange(year, month)[1]
        if strict_daily:
            found_days = set(dates.iloc[idxs].day.astype(int))
            missing_days = sorted(set(range(1, expected_days + 1)) - found_days)
            if ndays != expected_days or missing_days:
                raise ValueError(f'{year}-{month:02d}: missing daily records {missing_days}; '
                                 f'expected {expected_days} days, got {ndays}')
        elif ndays < args.min_days:
            raise ValueError(f'{year}-{month:02d}: only {ndays} distinct days; need {args.min_days}')
        if cfg['mode'] == 'nominal':
            mean, mean_outside = convert_one_month_nominal(
                ds, cfg, idxs, lats, lons, sigma, args.max_outside_pct)
        else:
            acc = np.zeros((len(sigma), len(lats), len(lons)), dtype=np.float64)
            outside = 0.0
            for i in idxs:
                v, frac = convert_one_day(ds, cfg, int(i), lats, lons, sigma,
                                          args.max_outside_pct)
                acc += v
                outside += frac
            mean = (acc / len(idxs)).astype(np.float32)
            mean_outside = outside / len(idxs)
        if not np.isfinite(mean).all() or (mean < 0).any():
            raise ValueError(f'{year}-{month:02d}: invalid monthly OH')
        print(f'  {year}-{month:02d} samples={len(idxs):3d}, days={ndays:2d}, '
              f'outside={mean_outside:.2f}%, mean OH={mean.mean():.4g}')
        monthly.append(mean)
    return monthly


def run(args):
    ref = Path(args.reference).expanduser().resolve()
    src_dir = Path(args.input_dir).expanduser().resolve()
    out = Path(args.output).expanduser().resolve()

    legacy = ref.parent / 'tch4_INCA_OH_t42.nc'
    if out == legacy.resolve():
        raise ValueError('Refusing to overwrite the original INCA OH')
    if out.exists() or out.with_name(out.name + '.part').exists():
        raise FileExistsError(f'Output already exists: {out} (or its .part). Never overwrite a run.')
    if args.start_year > args.end_year:
        raise ValueError('Invalid year range')
    if args.min_days < 1 or args.min_days > 31:
        raise ValueError('min-days must be 1..31')
    if args.extend_edges and (args.start_year != 2000 or args.end_year != 2024):
        raise ValueError('--extend-edges requires 2000–2024 source years')

    missing = [str(src_dir / f'INCA_DM_{year}_OH_scaled.nc')
               for year in range(args.start_year, args.end_year + 1)
               if not (src_dir / f'INCA_DM_{year}_OH_scaled.nc').is_file()]
    if missing:
        raise FileNotFoundError('Missing input years:\n' + '\n'.join(missing))

    with xr.open_dataset(ref, decode_times=False) as ds:
        lats = ds['y'].to_numpy().astype(float)
        lons = ds['x'].to_numpy().astype(float)
        sigma = ds['level'].to_numpy().astype(float)
    if not (np.isfinite(sigma).all() and (sigma > 0).all() and (sigma <= 1).all()):
        raise ValueError('ACTM levels must be sigma in (0,1]')

    # --- Preflight the first daily file before opening an output
    f0 = src_dir / f'INCA_DM_{args.start_year}_OH_scaled.nc'
    with xr.open_dataset(f0, decode_times=True) as ds:
        first_cfg = get_source_metadata(ds, args)
        reference_grid = {
            name: ds[first_cfg[name]].to_numpy().copy()
            for name in ('z', 'y', 'x')
        }
        reference_oh_factor = first_cfg['oh_factor']
    print('Vertical mode:', first_cfg['mode'])
    if first_cfg['mode'] == 'nominal':
        print('WARNING: reference pseudo-pressure approximation, NOT physical-pressure interpolation')
        print('Reference pressure [Pa]:', first_cfg['nominal_ref_pa'])

    out.parent.mkdir(parents=True, exist_ok=True)
    temp = out.with_name(out.name + '.part')
    writer = None
    try:
        writer = NCWriter(temp, lats, lons, sigma, first_cfg['mode'],
                          first_cfg.get('nominal_ref_pa'))
        n = 0
        for year in range(args.start_year, args.end_year + 1):
            fname = src_dir / f'INCA_DM_{year}_OH_scaled.nc'
            print('Processing:', fname)
            with xr.open_dataset(fname, decode_times=True) as ds:
                cfg = get_source_metadata(ds, args)
                if cfg['mode'] != first_cfg['mode'] or cfg['oh_factor'] != reference_oh_factor:
                    raise ValueError(f'{year}: inconsistent OH units or vertical mode between files')
                for name in ('z', 'y', 'x'):
                    current = ds[cfg[name]].to_numpy()
                    original = reference_grid[name]
                    if current.shape != original.shape or not np.allclose(current, original, rtol=0, atol=1e-5):
                        raise ValueError(f'{year}: source coordinate {name} changed between years')
                months = build_months(ds, cfg, year, lats, lons, sigma, args)

            if args.extend_edges and year == 2000:
                for edge_year in [1998, 1999]:
                    for m, v in enumerate(months, 1):
                        writer.write(edge_year, m, v)
                        n += 1
            for m, v in enumerate(months, 1):
                writer.write(year, m, v)
                n += 1
            if args.extend_edges and year == 2024:
                for edge_year in [2025, 2026]:
                    for m, v in enumerate(months, 1):
                        writer.write(edge_year, m, v)
                        n += 1

        writer.close()
        writer = None
        expected = 12 * (args.end_year - args.start_year + 1 + (4 if args.extend_edges else 0))
        if n != expected:
            raise ValueError(f'Month count mismatch: {n} vs {expected}')
        with xr.open_dataset(temp, engine='netcdf4', decode_times=True) as ds:
            if ds.sizes['time'] != expected or ds['OH'].shape[1:] != (len(sigma), len(lats), len(lons)):
                raise ValueError('Output validation failed')
            # --- Validate every written month, not just first/middle/last.
            times = pd.DatetimeIndex(ds['time'].to_numpy())
            if times.has_duplicates or not times.is_monotonic_increasing:
                raise ValueError('Output time coordinate is duplicated or unordered')
            if len(set((int(x.year), int(x.month)) for x in times)) != expected:
                raise ValueError('Output months are not unique')
            for ti in range(expected):
                a = ds['OH'].isel(time=ti).to_numpy()
                if not np.isfinite(a).all() or (a < 0).any():
                    raise ValueError(f'Invalid OH in output month {ti}')
            if args.extend_edges:
                for month in range(12):
                    first = ds['OH'].isel(time=24+month).to_numpy()  # 2000
                    last = ds['OH'].isel(time=expected-36+month).to_numpy()  # 2024
                    for year_offset in [0, 12]:  # 1998 and 1999
                        a = ds['OH'].isel(time=year_offset+month).to_numpy()
                        if not np.array_equal(a, first):
                            raise ValueError(f'1998/1999 extension differs at month {month+1}')
                    for year_offset in [expected-24, expected-12]:  # 2025 and 2026
                        a = ds['OH'].isel(time=year_offset+month).to_numpy()
                        if not np.array_equal(a, last):
                            raise ValueError(f'2025/2026 extension differs at month {month+1}')
        temp.replace(out)
        print('CONVERSION PASS:', out)
        print('Months:', n, 'Grid:', len(sigma), len(lats), len(lons))
    except BaseException:
        if writer is not None:
            writer.close()
        if temp.exists():
            temp.unlink()
        raise


def main():
    inp, ref, out = get_paths()
    parser = argparse.ArgumentParser(description='Strict, non-destructive INCA OH v2 conversion')
    parser.add_argument('--input-dir', default=str(inp))
    parser.add_argument('--reference', default=str(ref))
    parser.add_argument('--output', default=None, help='output NetCDF; default depends on vertical mode')
    parser.add_argument('--start-year', type=int, default=2000)
    parser.add_argument('--end-year', type=int, default=2024)
    parser.add_argument('--vertical-mode', choices=['physical', 'nominal'], default='physical',
                        help='physical: requires actual pmid and ps (default); '
                             'nominal: explicit approximate LMDZ pseudo-pressure mapping')
    parser.add_argument('--nominal-reference-pressure-pa', type=float, default=101325.0,
                        help='reference pressure for nominal mode only (Pa)')
    parser.add_argument('--ps-var', help='actual surface pressure variable name, if nonstandard')
    parser.add_argument('--assume-oh-unit', choices=['molecules/m3', 'molecules/cm3'])
    parser.add_argument('--assume-pressure-unit', choices=['Pa', 'hPa'])
    parser.add_argument('--min-days', type=int, default=20,
                        help='minimum distinct days if --allow-incomplete-months is set')
    parser.add_argument('--allow-incomplete-months', action='store_true',
                        help='EXPERIMENTAL: permit missing days (disabled by default)')
    parser.add_argument('--max-outside-pct', type=float, default=10.0)
    parser.add_argument('--extend-edges', action='store_true', help='1998=1999=2000; 2025=2026=2024')
    parser.add_argument('--inspect', action='store_true', help='read one daily file; no output')
    parser.add_argument('--inspect-year', type=int, default=2000)
    parser.add_argument('--run', action='store_true', help='explicitly create separate v2 NetCDF')
    args = parser.parse_args()
    if args.output is None:
        if args.vertical_mode == 'nominal':
            args.output = str(out.with_name('tch4_INCA_OH_t42_nominal.nc'))
        else:
            args.output = str(out)

    if args.run:
        run(args)
    else:
        path = Path(args.input_dir) / f'INCA_DM_{args.inspect_year}_OH_scaled.nc'
        ok = inspect_daily(path, args)
        if ok:
            print('No files changed. Use --run --extend-edges only after reviewing the conversion mode.')
        else:
            print('No files changed. Find companion pressure fields for physical conversion,')
            print('or inspect the approximate mode with --vertical-mode nominal.')


def run_nominal_full_period():
    # --- Use only for an EXPLICIT run from _main.py. No silent fallback to physical mode.
    from types import SimpleNamespace
    inp, ref, out = get_paths()
    dest = out.with_name('tch4_INCA_OH_t42_nominal_1998-2026.nc')
    args = SimpleNamespace(
        input_dir=str(inp), reference=str(ref), output=str(dest),
        start_year=2000, end_year=2024,
        vertical_mode='nominal', nominal_reference_pressure_pa=101325.0,
        ps_var=None, assume_oh_unit=None, assume_pressure_unit=None,
        min_days=20, max_outside_pct=10.0,
        allow_incomplete_months=False, extend_edges=True)
    run(args)





# =============================================================
# === GT3: preserved from the verified independent converter
import shutil
import struct
import subprocess
def gt3_paths():
    project = Path(__file__).resolve().parent.parent
    inp = project / 'fluxes/gcp2025_OH/tch4_INCA_OH_t42_nominal_1998-2026.nc'
    out = project / 'fluxes/gcp2025_OH/tch4_INCA_OH_t42_nominal_1998-2026.gt3'
    src = Path(__file__).resolve().parent / 'fortran/df0_inca_oh_nc2gt3_t42.f90'
    return inp, out, src


def inspect_nc(path):
    if not path.is_file():
        raise FileNotFoundError(f'Input OH NetCDF not found: {path}')

    with xr.open_dataset(path, decode_times=True) as ds:
        if 'OH' not in ds:
            raise ValueError('Input NetCDF has no OH variable')
        oh = ds['OH']
        if oh.dims != ('time', 'level', 'y', 'x'):
            raise ValueError(f'Unexpected OH dimensions: {oh.dims}')
        if oh.shape[1:] != (67, 64, 128):
            raise ValueError(f'Unexpected OH grid: {oh.shape}')
        units = str(oh.attrs.get('units', '')).strip()
        if units not in ('molecules cm-3', 'molec/cm3', 'molecules/cm3'):
            raise ValueError(f'Unexpected OH units: {units!r}')
        dates = pd.DatetimeIndex(ds['time'].values)
        if not len(dates) or len(dates) % 12:
            raise ValueError('Expected a complete number of years')
        expected = pd.period_range(dates[0].to_period('M'), periods=len(dates), freq='M')
        actual = pd.PeriodIndex(dates, freq='M')
        if not expected.equals(actual):
            raise ValueError('Missing, repeated or out-of-order OH months')
        if dates[0].month != 1 or dates[-1].month != 12:
            raise ValueError('OH time must begin in January and end in December')
        longitude = ds['x'].values
        if not np.allclose(longitude, np.arange(128)*360/128, atol=0.001, rtol=0):
            raise ValueError('Longitude is not the expected T42 grid')
        latitude = ds['y'].values
        if not (np.diff(latitude)>0).all() and not (np.diff(latitude)<0).all():
            raise ValueError('Latitude grid is not monotonic')
        level = ds['level'].values
        if not (np.isfinite(level).all() and ((level>0)&(level<=1)).all()):
            raise ValueError('Invalid ACTM sigma levels')
        print('INPUT OK:', path)
        print('Months:', len(dates), 'period:', actual[0], '..', actual[-1])
        print('Grid:', 'x'.join(str(v) for v in oh.shape[1:]), 'units:', units)
        return int(dates[0].year), len(dates)


# =============================================================
# === GT3 structural verification
# =============================================================

def verify_gt3(path, start_year, nmonths):
    head_bytes = 64*16
    data_bytes = 128*64*67*4
    expected_size = nmonths*(head_bytes+data_bytes+16)
    actual_size = path.stat().st_size
    if actual_size != expected_size:
        raise ValueError(f'GT3 size mismatch: {actual_size} vs {expected_size}')

    with path.open('rb') as f:
        for index in range(nmonths):
            head_len = struct.unpack('>I', f.read(4))[0]
            if head_len != head_bytes:
                raise ValueError(f'Invalid GT3 header record #{index}: {head_len}')
            hdr = f.read(head_bytes)
            end_head = struct.unpack('>I', f.read(4))[0]
            if end_head != head_bytes:
                raise ValueError(f'Invalid GT3 header trailer #{index}')
            words = [hdr[i:i+16].decode('ascii') for i in range(0,head_bytes,16)]
            year = start_year+index//12
            month = index%12+1
            if words[2].strip() != 'NDOH' or words[15].strip() != 'molec/cm3':
                raise ValueError(f'Unexpected GT3 OH header #{index}')
            if words[26][:6] != f'{year:04d}{month:02d}':
                raise ValueError(f'Wrong GT3 date #{index}: {words[26]!r}')
            if int(words[63]) != 128*64*67:
                raise ValueError(f'Wrong GT3 data size header #{index}')
            data_len = struct.unpack('>I', f.read(4))[0]
            if data_len != data_bytes:
                raise ValueError(f'Wrong GT3 data length #{index}')
            f.seek(data_bytes, 1)
            end_data = struct.unpack('>I', f.read(4))[0]
            if end_data != data_bytes:
                raise ValueError(f'Wrong GT3 data trailer #{index}')
        if f.read(1):
            raise ValueError('Unexpected trailing GT3 bytes')
    print('GT3 CHECK: PASS', nmonths, 'monthly records;', actual_size, 'bytes')


def build_fortran(source, binary, compiler):
    nf = shutil.which('nf-config')
    if not nf:
        raise RuntimeError('nf-config not found; load Intel NetCDF Fortran environment first')
    cc = shutil.which(compiler)
    if not cc:
        raise RuntimeError(f'Fortran compiler not found: {compiler}')
    flags = subprocess.check_output([nf, '--fflags'], text=True).split()
    libs = subprocess.check_output([nf, '--flibs'], text=True).split()
    opts = ['-O2', '-traceback'] if compiler == 'ifort' else ['-O2', '-fbacktrace']
    binary.parent.mkdir(parents=True,exist_ok=True)
    cmd = [cc]+opts+flags+[str(source),'-o',str(binary)]+libs
    print('Build:', ' '.join(cmd))
    subprocess.run(cmd,check=True,cwd=str(binary.parent))
    return binary




# =============================================================
# === D0: single large-stage entry point
def run_stage(case='INCA', dry_run=False):
    case = str(case).upper()
    if case != 'INCA':
        raise ValueError('D0 prepares INCA-OH only; use --case INCA')

    inp, ref, old_out = get_paths()
    nc = old_out.with_name('tch4_INCA_OH_t42_nominal_1998-2026.nc')
    gt3 = nc.with_suffix('.gt3')
    source = Path(__file__).resolve().parent / 'fortran/df0_inca_oh_nc2gt3_t42.f90'
    print('D0 INCA-OH: daily -> monthly T42L67 NetCDF -> GT3')
    print('NetCDF:', nc)
    print('GT3   :', gt3)
    print('Vertical mapping: nominal pseudo-pressure, not column pmid/Ps')

    if dry_run:
        print('DRY-RUN: no conversion and no output files changed')
        return

    if nc.is_file():
        start_year, nmonths = inspect_nc(nc)
        if start_year != 1998 or nmonths != 348:
            raise ValueError('Existing OH NetCDF is not the expected 1998-2026 period')
        print('NetCDF already exists and passed metadata checks: not overwriting')
    else:
        run_nominal_full_period()
        start_year, nmonths = inspect_nc(nc)
        if start_year != 1998 or nmonths != 348:
            raise ValueError('Generated OH NetCDF does not have 348 months, 1998-2026')

    if gt3.exists():
        verify_gt3(gt3, 1998, 348)
        print('GT3 already exists and passed structural checks: not overwriting')
        return

    if not source.is_file():
        raise FileNotFoundError(source)
    temp = gt3.with_suffix('.gt3.part')
    if temp.exists():
        raise FileExistsError('Incomplete GT3 already exists: ' + str(temp))
    binary = source.parent / 'bin/df0_inca_oh_nc2gt3_t42'
    import os as _os
    compiler = _os.environ.get('FC', 'ifort')
    build_fortran(source, binary, compiler)
    try:
        subprocess.run([str(binary), str(nc), str(temp), '1998'], check=True)
        verify_gt3(temp, 1998, 348)
        temp.rename(gt3)
    except BaseException:
        if temp.exists():
            temp.unlink()
        raise
    print('D0 COMPLETE:', gt3)


if __name__ == '__main__':
    main()
