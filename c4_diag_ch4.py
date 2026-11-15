# -*- coding: utf-8 -*-
"""
Project:    GCP/GMB 2025 project
Purpose:    F6/F7 CH4 scientific diagnostics, one plot per active station.
Notes:      Source formats follow f6_nakazawa_filter.f90 and f7_write_tdi_obsrvCH4dif.f90.
"""
import os
import argparse
import calendar
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import netcdf_file


# ============================================================
# === Constants from Fortran (do not silently infer a different format)
N_SITE = 245
YEAR_S = 1999
YEAR_E = 2025
N_DAY = sum(366 if calendar.isleap(y) else 365 for y in range(YEAR_S, YEAR_E + 1))
N_MONTH = (YEAR_E - YEAR_S + 1) * 12
FILL = -999.0
LAG_WEAK_MARGIN = 0.05  # diagnostic heuristic, not a significance test


# ============================================================
# === Dates and site metadata
def get_dates():
    days = pd.date_range(str(YEAR_S) + '-01-01', str(YEAR_E) + '-12-31', freq='D')
    months = pd.date_range(str(YEAR_S) + '-01-01', str(YEAR_E) + '-12-01', freq='MS')
    if len(days) != N_DAY or len(months) != N_MONTH:
        raise ValueError('F6/F7 calendar does not match 1999-2025')
    return days, months


def get_sites(filename):
    # --- Fortran F7 skips the first two lines of obsrvCH4_60sites.txt.
    rows = []
    for line in Path(filename).read_text(encoding='utf-8', errors='replace').splitlines()[2:]:
        fields = line.split()
        if len(fields) < 6:
            continue
        try:
            station_no = int(fields[0])
            flag = int(fields[1])
            lat = float(fields[3])
            lon = float(fields[4])
            alt = float(fields[5])
        except ValueError:
            continue
        rows.append((station_no, flag, fields[2], lat, lon, alt))
    if len(rows) != N_SITE:
        raise ValueError('Expected 245 station rows, found ' + str(len(rows)) + ': ' + str(filename))
    sites = pd.DataFrame(rows, columns=['station_id', 'flag', 'site', 'lat', 'lon', 'alt'])
    sites['index'] = np.arange(N_SITE)
    sites['active'] = sites['flag'] == 9999
    if sites['station_id'].duplicated().any():
        raise ValueError('Duplicate station IDs in ' + str(filename))
    # --- Repeated labels are present in the original 245-site file (e.g. MLO, POC).
    #     Station number and row position, not site_name alone, identify a station.
    if not sites['site'].all():
        raise ValueError('Empty station name in ' + str(filename))
    return sites


# ============================================================
# === Fortran direct-access float32 binaries, without record markers
#     Input F6:  record(day) = (ntrac=2, nsite=245) => 1960 bytes/day
#     Output F6: record(day) = (ntrac=2, kind=3, nsite=245) => 5880 bytes/day

def get_binary(filename, shape, endian='<'):
    filename = Path(filename)
    if not filename.is_file():
        raise FileNotFoundError(filename)
    dtype = np.dtype(endian + 'f4')
    expected = int(np.prod(shape, dtype=np.int64)) * dtype.itemsize
    actual = filename.stat().st_size
    if actual != expected:
        raise ValueError(str(filename) + ': wrong size: ' + str(actual) + ' bytes; expected ' + str(expected))
    return np.memmap(filename, dtype=dtype, mode='r', shape=shape, order='C')


def get_f6_files(input_file, output_file, endian='<'):
    inp = get_binary(input_file, (N_DAY, 2, N_SITE), endian)
    out = get_binary(output_file, (N_DAY, 2, 3, N_SITE), endian)
    return inp, out


# ============================================================
# === F7 NetCDF classic. The Python axis order is (month, station, type).
def get_name(row):
    if isinstance(row, bytes):
        return row.decode('ascii', errors='replace').strip(' \x00')
    return b''.join(row).decode('ascii', errors='replace').strip(' \x00')


def get_f7_file(filename):
    filename = Path(filename)
    if not filename.is_file():
        raise FileNotFoundError(filename)
    try:
        with netcdf_file(str(filename), mode='r', mmap=False) as nc:
            dims = dict(nc.dimensions)
            for name, length in [('station_count', N_SITE), ('tlen', N_MONTH), ('data_type', 2), ('namelen', 45)]:
                if dims.get(name) != length:
                    raise ValueError(str(filename) + ': wrong dimension ' + name + ': ' + str(dims.get(name)))
            for name in ['site_name', 'station_count', 'tlen', 'tracer_data']:
                if name not in nc.variables:
                    raise ValueError(str(filename) + ': missing variable ' + name)
            var = nc.variables['tracer_data']
            if tuple(var.dimensions) != ('tlen', 'station_count', 'data_type'):
                raise ValueError(str(filename) + ': wrong tracer_data dimension order: ' + str(var.dimensions))
            data = np.array(var[:], dtype='f8')
            numbers = np.array(nc.variables['station_count'][:], dtype='i4')
            names = [get_name(row) for row in nc.variables['site_name'][:]]
            time = np.array(nc.variables['tlen'][:], dtype='f8')
    except (TypeError, ValueError) as exc:
        if isinstance(exc, ValueError) and 'wrong ' in str(exc):
            raise
        raise ValueError(str(filename) + ': cannot read as NetCDF classic. F7 uses nf_create(NF_CLOBBER). ' + str(exc)) from exc
    if data.shape != (N_MONTH, N_SITE, 2):
        raise ValueError(str(filename) + ': unexpected tracer_data shape ' + str(data.shape))
    if not all(names):
        raise ValueError(str(filename) + ': empty NetCDF station name')
    # --- Duplicate site_name values are expected for some distinct station IDs.
    if len(set(numbers)) != N_SITE:
        raise ValueError(str(filename) + ': duplicated station numbers')
    expect = YEAR_S + (np.arange(N_MONTH) + 0.5) / 12.0
    if not np.allclose(time, expect, atol=0.0003, rtol=0):
        raise ValueError(str(filename) + ': time values do not match 1999-2025 monthly midpoints')
    return {'filename': str(filename), 'data': data, 'numbers': numbers, 'names': names, 'time': time}


def get_valid(value, low=1201.0, high=5001.0):
    return np.isfinite(value) & (value >= low) & (value <= high)


def get_f7_valid(data):
    return np.isfinite(data[:, :, 0]) & np.isfinite(data[:, :, 1]) & (data[:, :, 0] != FILL) & (data[:, :, 1] > 0.0)


# ============================================================
# === Monthly reconstruction using EXACT F7 pairing and RMS rules
#     F7 uses fitted OBS and ACTM on common valid days for the difference.
#     Monthly RMS(raw OBS - fitted OBS) is pooled over the same calendar
#     month across all 27 years, then uncertainty = RMS * scale + 5 ppb.
#     If raw monthly observation is unavailable, uncertainty doubles.

def get_f6_monthly(fits, days=None, months=None):
    if days is None or months is None:
        days, months = get_dates()
    if fits.shape != (len(days), 2, 3, N_SITE):
        raise ValueError('F6 fits shape must be (days, 2, 3, 245)')

    nmonth = len(months)
    fit_diff = np.full((nmonth, N_SITE), np.nan, dtype='f8')
    raw_diff = np.full((nmonth, N_SITE), np.nan, dtype='f8')
    obs_raw_mean = np.full((nmonth, N_SITE), np.nan, dtype='f8')
    obs_fit_mean = np.full((nmonth, N_SITE), np.nan, dtype='f8')
    mod_fit_mean = np.full((nmonth, N_SITE), np.nan, dtype='f8')
    paired_fit_days = np.zeros((nmonth, N_SITE), dtype='i4')
    paired_raw_days = np.zeros((nmonth, N_SITE), dtype='i4')
    sumsq = np.zeros((12, N_SITE), dtype='f8')
    nsum = np.zeros((12, N_SITE), dtype='i4')

    offset = 0
    for nm, date in enumerate(months):
        count = calendar.monthrange(date.year, date.month)[1]
        data = fits[offset:offset + count]
        if len(data) != count:
            raise ValueError('F6 daily file ends before ' + str(date))
        obs_raw = np.asarray(data[:, 0, 0, :], dtype='f8')
        mod_raw = np.asarray(data[:, 1, 0, :], dtype='f8')
        obs_fit = np.asarray(data[:, 0, 1, :], dtype='f8')
        mod_fit = np.asarray(data[:, 1, 1, :], dtype='f8')
        mask_fit = get_valid(obs_fit) & get_valid(mod_fit)
        mask_raw = get_valid(obs_raw) & get_valid(mod_raw)
        nf = mask_fit.sum(axis=0)
        nr = mask_raw.sum(axis=0)
        paired_fit_days[nm] = nf
        paired_raw_days[nm] = nr
        ok_fit = nf > 0
        ok_raw = nr > 0
        if np.any(ok_fit):
            mean_o = np.divide(np.where(mask_fit, obs_fit, 0).sum(axis=0), nf, out=np.full(N_SITE, np.nan), where=ok_fit)
            mean_m = np.divide(np.where(mask_fit, mod_fit, 0).sum(axis=0), nf, out=np.full(N_SITE, np.nan), where=ok_fit)
            obs_fit_mean[nm] = mean_o
            mod_fit_mean[nm] = mean_m
            fit_diff[nm] = mean_o - mean_m
        if np.any(ok_raw):
            mean_o = np.divide(np.where(mask_raw, obs_raw, 0).sum(axis=0), nr, out=np.full(N_SITE, np.nan), where=ok_raw)
            mean_m = np.divide(np.where(mask_raw, mod_raw, 0).sum(axis=0), nr, out=np.full(N_SITE, np.nan), where=ok_raw)
            obs_raw_mean[nm] = mean_o
            raw_diff[nm] = mean_o - mean_m
        # --- F7 calculates stdev using OBS raw+fit, regardless of ACTM.
        residual_mask = get_valid(obs_raw) & get_valid(obs_fit)
        residual = np.where(residual_mask, obs_raw - obs_fit, 0.0)
        sumsq[date.month - 1] += np.square(residual).sum(axis=0)
        nsum[date.month - 1] += residual_mask.sum(axis=0)
        offset += count

    if offset != len(days):
        raise ValueError('F6 calendar mismatch: ' + str(offset) + ' / ' + str(len(days)))
    rms = np.sqrt(np.divide(sumsq, nsum, out=np.full_like(sumsq, np.nan), where=nsum > 0))
    return {'fit_diff': fit_diff, 'raw_diff': raw_diff, 'obs_raw_mean': obs_raw_mean,
            'obs_fit_mean': obs_fit_mean, 'mod_fit_mean': mod_fit_mean, 'rms': rms,
            'paired_fit_days': paired_fit_days, 'paired_raw_days': paired_raw_days}


def get_expected_f7(monthly, scale, months=None):
    if months is None:
        _, months = get_dates()
    diff = monthly['fit_diff']
    idx = np.asarray(months.month) - 1
    sigma = monthly['rms'][idx] * float(scale) + 5.0
    valid_raw_month = get_valid(monthly['obs_raw_mean'])
    sigma = np.where(valid_raw_month, sigma, sigma * 2.0)
    valid = get_valid(monthly['obs_fit_mean']) & get_valid(monthly['mod_fit_mean']) & np.isfinite(sigma) & (sigma >= 0)
    expected = np.full((len(months), N_SITE, 2), FILL, dtype='f8')
    expected[:, :, 0][valid] = diff[valid]
    expected[:, :, 1][valid] = sigma[valid]
    return expected


# ============================================================
# === Seasonal phase using Nakazawa fitted minus trend, paired OBS/ACTM days
#     The detrended monthly climatologies avoid confusing the secular CH4
#     growth trend with a shift in the seasonal cycle.
def get_f6_phase(fits, sites, days=None, min_days_per_calendar_month=60):
    if days is None:
        days, _ = get_dates()
    m_idx = np.asarray(days.month)
    obs_clim = np.full((12, N_SITE), np.nan, dtype='f8')
    mod_clim = np.full((12, N_SITE), np.nan, dtype='f8')
    n_days = np.zeros((12, N_SITE), dtype='i4')

    for month in range(1, 13):
        chunk = fits[m_idx == month]
        obs_fit = np.asarray(chunk[:, 0, 1, :], dtype='f8')
        mod_fit = np.asarray(chunk[:, 1, 1, :], dtype='f8')
        obs_trend = np.asarray(chunk[:, 0, 2, :], dtype='f8')
        mod_trend = np.asarray(chunk[:, 1, 2, :], dtype='f8')
        valid = get_valid(obs_fit) & get_valid(mod_fit) & get_valid(obs_trend) & get_valid(mod_trend)
        n = valid.sum(axis=0)
        n_days[month - 1] = n
        good = n >= min_days_per_calendar_month
        obs_clim[month - 1] = np.divide(np.where(valid, obs_fit - obs_trend, 0).sum(axis=0), n,
                                         out=np.full(N_SITE, np.nan), where=good)
        mod_clim[month - 1] = np.divide(np.where(valid, mod_fit - mod_trend, 0).sum(axis=0), n,
                                         out=np.full(N_SITE, np.nan), where=good)

    records = []
    for i, row in sites.iterrows():
        a, b = obs_clim[:, i], mod_clim[:, i]
        profile = get_lag_correlations(a, b)
        ranked = sorted(((lag, corr) for lag, corr in profile.items() if np.isfinite(corr)),
                        key=lambda v: v[1], reverse=True)
        lag = ranked[0][0] if ranked else np.nan
        corr = ranked[0][1] if ranked else np.nan
        second = ranked[1][1] if len(ranked) > 1 else np.nan
        margin = corr - second if np.isfinite(second) else np.nan
        record = {'station_id': int(row.station_id), 'site': row.site, 'active': bool(row.active),
                  'lat': float(row.lat), 'lon': float(row.lon),
                  'n_calendar_months': int(np.sum(np.isfinite(a) & np.isfinite(b))),
                  'min_days_per_month': int(n_days[:, i].min()),
                  'seasonal_lag_months': lag, 'best_correlation': corr,
                  'second_correlation': second, 'lag_margin': margin,
                  'lag_weak': bool(np.isfinite(margin) and margin < LAG_WEAK_MARGIN),
                  'obs_seasonal_amplitude_ppb': (float(np.max(a) - np.min(a)) if np.isfinite(a).all() else np.nan),
                  'actm_seasonal_amplitude_ppb': (float(np.max(b) - np.min(b)) if np.isfinite(b).all() else np.nan)}
        for lag_number in range(-3, 4):
            record['corr_' + str(lag_number)] = profile[lag_number]
        records.append(record)
    return pd.DataFrame(records)


# ============================================================
# === Statistical helpers
def get_statistics(values):
    v = np.asarray(values, dtype='f8')
    v = v[np.isfinite(v)]
    if not len(v):
        return (np.nan, np.nan, np.nan, np.nan)
    return (float(np.mean(v)), float(np.std(v)), float(np.min(v)), float(np.max(v)))


def get_lag_correlations(obs, mod):
    # --- Positive lag: ACTM seasonal cycle occurs later than OBS.
    #     Cyclic 12-month shift; seven candidate lags (-3 .. +3).
    obs = np.asarray(obs, dtype='f8')
    mod = np.asarray(mod, dtype='f8')
    result = {lag: np.nan for lag in range(-3, 4)}
    if obs.shape != (12,) or mod.shape != (12,):
        raise ValueError('Seasonal climatology must contain 12 months')
    if not np.isfinite(obs).all() or not np.isfinite(mod).all():
        return result
    a = obs - np.mean(obs)
    b = mod - np.mean(mod)
    if np.std(a) < 0.1 or np.std(b) < 0.1:
        return result
    for lag in range(-3, 4):
        result[lag] = float(np.corrcoef(a, np.roll(b, -lag))[0, 1])
    return result


def get_best_lag(obs, mod):
    profile = get_lag_correlations(obs, mod)
    finite = [(lag, corr) for lag, corr in profile.items() if np.isfinite(corr)]
    return max(finite, key=lambda p: p[1]) if finite else (np.nan, np.nan)


def get_f6_stats(inp, fits, sites):
    records = []
    for i, row in sites.iterrows():
        a = np.asarray(inp[:, 0, i], dtype='f8')
        b = np.asarray(inp[:, 1, i], dtype='f8')
        a_out = np.asarray(fits[:, 0, 0, i], dtype='f8')
        b_out = np.asarray(fits[:, 1, 0, i], dtype='f8')
        af = np.asarray(fits[:, 0, 1, i], dtype='f8')
        bf = np.asarray(fits[:, 1, 1, i], dtype='f8')
        at = np.asarray(fits[:, 0, 2, i], dtype='f8')
        bt = np.asarray(fits[:, 1, 2, i], dtype='f8')
        mask = get_valid(af) & get_valid(bf)
        raw = get_valid(a) & get_valid(b)
        match_a = get_valid(a_out)
        match_b = get_valid(b_out)
        preserved_a = match_a & (np.abs(a - a_out) > 0.001)
        preserved_b = match_b & (np.abs(b - b_out) > 0.001)
        resid = get_valid(a_out) & get_valid(af)
        stats = get_statistics((af - bf)[mask])
        rmse = np.sqrt(np.mean((a_out[resid] - af[resid]) ** 2)) if resid.any() else np.nan
        records.append({'station_id': row.station_id, 'site': row.site, 'active': bool(row.active),
                        'lat': row.lat, 'lon': row.lon,
                        'n_input_obs': int(get_valid(a, 900, 5000).sum()),
                        'n_input_model': int(get_valid(b, 900, 5000).sum()),
                        'n_raw_pairs': int(raw.sum()), 'n_fit_pairs': int(mask.sum()),
                        'n_obs_trend': int(get_valid(at).sum()), 'n_model_trend': int(get_valid(bt).sum()),
                        'n_raw_preservation_errors': int(preserved_a.sum() + preserved_b.sum()),
                        'mean_fit_obs_minus_model': stats[0], 'std_fit_obs_minus_model': stats[1],
                        'obs_raw_fit_rmse': float(rmse)})
    return pd.DataFrame(records)


def get_f7_stats(nc, sites):
    numbers = nc['numbers']
    names = nc['names']
    data = nc['data']
    if list(numbers) != list(sites.station_id) or names != list(sites.site):
        raise ValueError(nc['filename'] + ': NetCDF station order/names differ from site list')
    valid = get_f7_valid(data)
    summary = []
    for i, row in sites.iterrows():
        v = valid[:, i]
        mean, std, low, high = get_statistics(data[:, i, 0][v])
        sig = get_statistics(data[:, i, 1][v])
        summary.append({'file': Path(nc['filename']).name, 'station_id': row.station_id,
                        'site': row.site, 'active': bool(row.active), 'lat': row.lat,
                        'n_valid': int(v.sum()), 'n_missing': int((~v).sum()),
                        'mean_obs_minus_model': mean, 'std_obs_minus_model': std,
                        'min_obs_minus_model': low, 'max_obs_minus_model': high,
                        'mean_uncertainty': sig[0], 'max_uncertainty': sig[3],
                        'n_uncertainty_gt20': int(np.sum(v & (data[:, i, 1] > 20))),
                        'n_uncertainty_le20': int(np.sum(v & (data[:, i, 1] <= 20)))})
    return pd.DataFrame(summary)


def get_f7_monthly(nc, sites, months):
    data = nc['data']
    valid = get_f7_valid(data)
    active = sites.active.to_numpy(dtype=bool)
    summary = []
    for i, date in enumerate(months):
        for label, group in [('active', active), ('all', np.ones(N_SITE, dtype=bool))]:
            selected = valid[i] & group
            d = data[i, :, 0][selected]
            u = data[i, :, 1][selected]
            summary.append({'file': Path(nc['filename']).name, 'month': date.strftime('%Y-%m'),
                            'group': label, 'n_stations': int(group.sum()), 'n_valid': int(selected.sum()),
                            'mean_obs_minus_model': get_statistics(d)[0], 'median_obs_minus_model': float(np.median(d)) if len(d) else np.nan,
                            'mean_uncertainty': get_statistics(u)[0]})
    return pd.DataFrame(summary)


def get_comparison(nc, expected, sites, months, tolerance=0.01):
    actual = nc['data']
    v_actual = get_f7_valid(actual)
    v_expected = get_f7_valid(expected)
    missing_mismatch = v_actual != v_expected
    numerical = v_actual & v_expected
    delta = np.abs(actual - expected)
    bad0 = numerical & (delta[:, :, 0] > tolerance)
    bad1 = numerical & (delta[:, :, 1] > tolerance)
    bad = missing_mismatch | bad0 | bad1
    nbad = int(bad.sum())
    where = np.argwhere(bad)
    examples = []
    for month_i, site_i in where[:200]:
        row = sites.iloc[site_i]
        examples.append({'file': Path(nc['filename']).name, 'month': months[month_i].strftime('%Y-%m'),
                         'site': row.site, 'active': bool(row.active),
                         'actual_difference': actual[month_i, site_i, 0], 'expected_difference': expected[month_i, site_i, 0],
                         'actual_uncertainty': actual[month_i, site_i, 1], 'expected_uncertainty': expected[month_i, site_i, 1]})
    info = {'file': Path(nc['filename']).name, 'n_missing_pattern_mismatch': int(missing_mismatch.sum()),
            'n_difference_mismatch': int(bad0.sum()), 'n_uncertainty_mismatch': int(bad1.sum()),
            'n_total_bad_months': nbad, 'tolerance_ppb': tolerance,
            'max_abs_diff_ppb': float(delta[:, :, 0][numerical].max()) if numerical.any() else np.nan,
            'max_abs_uncertainty_ppb': float(delta[:, :, 1][numerical].max()) if numerical.any() else np.nan}
    return info, examples


# ============================================================
# === Legacy inversion reads every month (including FILL) from firstdat
#     through lastdat and resolves duplicated names to their LAST match.
#     The output-NetCDF interval (2000-2024) is unrelated to this period.
def get_inversion_period(setup, control_file=None):
    path = control_file if control_file is not None else getattr(setup, 'inv_control_template', None)
    if path:
        path = Path(path).expanduser().resolve()
        if path.is_file():
            content = path.read_text(encoding='utf-8', errors='replace')
            values = {}
            for key in ('firstdat', 'lastdat'):
                match = re.search(r'(?im)^\s*' + key + r'\s*=\s*(\d+)\b', content)
                if not match:
                    raise ValueError(str(path) + ': missing ' + key)
                values[key] = int(match.group(1))
            first, last = values['firstdat'], values['lastdat']
            source = 'control: ' + str(path)
        elif control_file is not None:
            raise FileNotFoundError(path)
        else:
            first = int(getattr(setup, 'inv_year_s', getattr(setup, 'yr_s', YEAR_S)))
            last = int(getattr(setup, 'inv_year_e', getattr(setup, 'yr_e', YEAR_E)))
            source = 'SetCase fallback (control missing: ' + str(path) + ')'
    else:
        first = int(getattr(setup, 'inv_year_s', getattr(setup, 'yr_s', YEAR_S)))
        last = int(getattr(setup, 'inv_year_e', getattr(setup, 'yr_e', YEAR_E)))
        source = 'SetCase fallback (no control configured)'
    if not (YEAR_S <= first <= last <= YEAR_E):
        raise ValueError('Inversion data period outside F7 1999-2025: ' + str(first) + '-' + str(last))
    return first, last, source


def get_inversion_missing(nc, sites, months, firstdat, lastdat):
    # --- Reproduce the station-name lookup in read_globalview.f:
    #     the loop retains the LAST matching site_name, not the station ID.
    name_to_last_index = {name: i for i, name in enumerate(nc['names'])}
    start = (firstdat - YEAR_S) * 12
    end = (lastdat - YEAR_S + 1) * 12
    valid = get_f7_valid(nc['data'])
    result = {}
    examples = []
    n_invalid = 0
    n_fill = 0
    for i, row in sites.loc[sites.active].iterrows():
        if row.site not in name_to_last_index:
            raise ValueError('Inversion active name not found in F7: ' + str(row.site))
        j = name_to_last_index[row.site]
        ok = valid[start:end, j]
        values = nc['data'][start:end, j, :]
        missing = int((~ok).sum())
        fill = int(((values == FILL).any(axis=1)).sum())
        result[int(row.station_id)] = {'missing': missing, 'fill': fill,
                                       'netcdf_station_id_used': int(sites.iloc[j].station_id)}
        n_invalid += missing
        n_fill += fill
        if missing and len(examples) < 12:
            bad_month = months[start + np.flatnonzero(~ok)[0]].strftime('%Y-%m')
            examples.append((int(row.station_id), str(row.site), bad_month, missing))
    return {'by_id': result, 'n_invalid': n_invalid, 'n_fill': n_fill,
            'examples': examples, 'n_requested': int(sites.active.sum()),
            'n_months': end - start, 'period': str(firstdat) + '-' + str(lastdat)}


def get_duplicate_active_sites(sites, nc=None):
    groups = []
    for name, rows in sites.loc[sites.active].groupby('site', sort=True):
        if len(rows) < 2:
            continue
        ids = rows.station_id.astype(int).to_list()
        identical = None
        if nc is not None:
            indices = rows.index.to_list()
            base = nc['data'][:, indices[0], :]
            identical = all(np.array_equal(base, nc['data'][:, j, :], equal_nan=True) for j in indices[1:])
        groups.append({'site': name, 'ids': ids, 'rows_identical': identical,
                       'inversion_last_id': ids[-1]})
    return groups


# ============================================================
# === Detect numerically identical F7 time series even with DIFFERENT names.
#     Compare both OBS-ACTM and uncertainty for every month, including masks.
#     This is independent of the station-name duplication check above.
def get_identical_active_series(sites, nc, min_valid=12):
    if nc is None:
        return []
    data = nc['data']
    if data.shape != (N_MONTH, len(sites), 2):
        raise ValueError('Unexpected F7 series shape in duplicate check: ' + str(data.shape))
    active = sites.index[sites.active].to_list()
    valid = get_f7_valid(data)
    duplicates = []
    for pos, i in enumerate(active):
        if int(valid[:, i].sum()) < min_valid:
            continue
        for j in active[pos + 1:]:
            if int(valid[:, j].sum()) < min_valid:
                continue
            same_mask = np.array_equal(valid[:, i], valid[:, j])
            if not same_mask:
                continue
            same_difference = np.array_equal(data[:, i, 0], data[:, j, 0], equal_nan=True)
            if not same_difference:
                continue
            same_uncertainty = np.array_equal(data[:, i, 1], data[:, j, 1], equal_nan=True)
            duplicates.append({'id_a': int(sites.iloc[i].station_id),
                               'id_b': int(sites.iloc[j].station_id),
                               'name_a': str(sites.iloc[i].site),
                               'name_b': str(sites.iloc[j].site),
                               'n_valid': int(valid[:, i].sum()),
                               'same_uncertainty': bool(same_uncertainty),
                               'same_name': bool(sites.iloc[i].site == sites.iloc[j].site)})
    return duplicates


# ============================================================
# === Select stations by inclusive STATION IDs, not file row positions
#     --stations 1 42 means all active stations from 1 through 42.
def get_plot_indices(sites, selectors=None):
    if selectors is None:
        subset = sites.loc[sites.active]
    else:
        if len(selectors) != 2:
            raise ValueError('--stations requires exactly two station IDs: START END')
        start, end = map(int, selectors)
        lo = int(sites.station_id.min())
        hi = int(sites.station_id.max())
        if not (lo <= start <= end <= hi):
            raise ValueError('Station range must satisfy ' + str(lo) + ' <= START <= END <= ' + str(hi))
        subset = sites.loc[sites.active & sites.station_id.between(start, end)]
    return subset.index.to_list()


# ============================================================
# === A station's 12-month climatology after subtracting Nakazawa trend
#     Same paired-days approach and 60-day minimum as get_f6_phase().
def get_station_climatology(fits, days, site_i, min_days=60):
    obs = np.asarray(fits[:, 0, 1, site_i], dtype='f8')
    mod = np.asarray(fits[:, 1, 1, site_i], dtype='f8')
    otr = np.asarray(fits[:, 0, 2, site_i], dtype='f8')
    mtr = np.asarray(fits[:, 1, 2, site_i], dtype='f8')
    valid = get_valid(obs) & get_valid(mod) & get_valid(otr) & get_valid(mtr)
    mm = np.asarray(days.month)
    a = np.full(12, np.nan)
    b = np.full(12, np.nan)
    n = np.zeros(12, dtype='i4')
    for month in range(1, 13):
        ok = valid & (mm == month)
        n[month - 1] = int(ok.sum())
        if n[month - 1] >= min_days:
            a[month - 1] = np.mean(obs[ok] - otr[ok])
            b[month - 1] = np.mean(mod[ok] - mtr[ok])
    return a, b, n


# ============================================================
# === Prefer ctl/inv1 as the single comparison series in panel D.
#     ctl / ux2 / ux4 have the same concentration differences by design.
def get_primary_f7(files):
    if not files:
        return None
    names = sorted(files, key=lambda name: (0 if 'ctlinv1' in name.lower() else
                                            1 if 'ctl' in name.lower() else 2, name))
    return files[names[0]]


def get_uncertainty_files(files):
    found = {}
    for unc in ('ctl', 'ux2', 'ux4'):
        matches = sorted(name for name in files if unc in name.lower())
        if matches:
            found[unc] = files[matches[0]]
    return found


# ============================================================
# === Station plot: 3 x 2 panels, full inversion period
#     Raw points and fitted curves are deliberately visually distinct.
def plot_station(outdir, site_i, sites, days, months, inp, fits, monthly, files, hcase='CYC'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    row = sites.iloc[site_i]
    site = str(row.site)
    station_id = int(row.station_id)
    safe = ''.join(ch if ch.isalnum() or ch in '-_' else '_' for ch in site)[:75]
    figure, ax = plt.subplots(3, 2, figsize=(18, 12), constrained_layout=True)

    # --- A: Full 1999-2025 daily series, four original/fitted curves
    raw_obs = np.asarray(inp[:, 0, site_i], dtype='f8')
    raw_mod = np.asarray(inp[:, 1, site_i], dtype='f8')
    fit_obs = np.asarray(fits[:, 0, 1, site_i], dtype='f8')
    fit_mod = np.asarray(fits[:, 1, 1, site_i], dtype='f8')
    raw_obs[~get_valid(raw_obs)] = np.nan
    raw_mod[~get_valid(raw_mod)] = np.nan
    fit_obs[~get_valid(fit_obs)] = np.nan
    fit_mod[~get_valid(fit_mod)] = np.nan

    ax[0, 0].plot(days, fit_obs, lw=1.05, color='tab:red', label='OBS fitted', zorder=3)
    ax[0, 0].plot(days, fit_mod, lw=1.05, color='tab:blue', label='ACTM fitted', zorder=3)
    ax[0, 0].scatter(days, raw_obs, s=2, alpha=0.35, color='black', label='OBS raw', rasterized=True, zorder=4)
    ax[0, 0].scatter(days, raw_mod, s=2, alpha=0.22, color='gray', label='ACTM raw', rasterized=True, zorder=2)
    ax[0, 0].set_title('A. Raw observations/model and Nakazawa fits')
    ax[0, 0].set_ylabel('CH4 (ppb)')
    ax[0, 0].legend(ncol=2, fontsize=8)

    # --- B: Detrended monthly climatology and relative seasonal lag
    o_clim, m_clim, n_clim = get_station_climatology(fits, days, site_i)
    lag, corr = get_best_lag(o_clim, m_clim)
    xx = np.arange(1, 13)
    ax[0, 1].plot(xx, o_clim, '-o', ms=3, color='tab:red', label='OBS fitted - trend')
    ax[0, 1].plot(xx, m_clim, '-o', ms=3, color='tab:blue', label='ACTM fitted - trend')
    ax[0, 1].axhline(0, color='gray', lw=0.7)
    ax[0, 1].set_xticks(xx)
    ax[0, 1].set_xticklabels(list(calendar.month_abbr)[1:])
    ax[0, 1].set_title('B. Seasonal climatology (ACTM lag: ' +
                       (str(int(lag)) + ' mo' if np.isfinite(lag) else 'N/A') + ')')
    ax[0, 1].set_ylabel('CH4 seasonal anomaly (ppb)')
    ax[0, 1].legend(fontsize=8)

    # --- C: Residual of Nakazawa fit to raw OBS
    residual = raw_obs - fit_obs
    residual[~(np.isfinite(raw_obs) & np.isfinite(fit_obs))] = np.nan
    err = residual[np.isfinite(residual)]
    rmse = np.sqrt(np.mean(err ** 2)) if len(err) else np.nan
    ax[1, 0].scatter(days, residual, s=2, alpha=0.5, color='tab:purple', rasterized=True)
    ax[1, 0].axhline(0, color='gray', lw=0.7)
    ax[1, 0].set_title('C. OBS raw - OBS fitted | RMSE = ' +
                       (f'{rmse:.2f}' if np.isfinite(rmse) else 'N/A') + ' ppb')
    ax[1, 0].set_ylabel('Filter residual (ppb)')

    # --- D: Monthly OBS - ACTM passed to inversion; compare with reconstruction
    primary = get_primary_f7(files)
    rec = np.asarray(monthly['fit_diff'][:, site_i], dtype='f8')
    ax[1, 1].plot(months, rec, lw=1, color='tab:orange', label='From F6 fitted values')
    if primary is not None:
        obs_minus_mod = np.asarray(primary['data'][:, site_i, 0], dtype='f8').copy()
        obs_minus_mod[~get_f7_valid(primary['data'])[:, site_i]] = np.nan
        ax[1, 1].plot(months, obs_minus_mod, '.', ms=2.5, color='black', label='F7 NetCDF')
    ax[1, 1].axhline(0, color='gray', lw=0.7)
    ax[1, 1].set_title('D. Monthly OBS - ACTM (F6 reconstruction vs F7)')
    ax[1, 1].set_ylabel('OBS - ACTM (ppb)')
    ax[1, 1].legend(fontsize=8)

    # --- E: Observation uncertainty by F7 inversion option
    uncertainties = get_uncertainty_files(files)
    colors = {'ctl': 'tab:green', 'ux2': 'tab:orange', 'ux4': 'tab:red'}
    for unc, nc in uncertainties.items():
        sigma = np.asarray(nc['data'][:, site_i, 1], dtype='f8').copy()
        sigma[~get_f7_valid(nc['data'])[:, site_i]] = np.nan
        ax[2, 0].plot(months, sigma, lw=0.9, color=colors[unc], label=unc)
    ax[2, 0].axhline(20, ls='--', lw=0.9, color='gray', label='20 ppb')
    ax[2, 0].set_title('E. Monthly F7 uncertainty')
    ax[2, 0].set_ylabel('Uncertainty (ppb)')
    ax[2, 0].legend(fontsize=8)

    # --- F: Real daily support for monthly values (0..31 days)
    ax[2, 1].plot(months, monthly['paired_raw_days'][:, site_i],
                  lw=0.85, label='Raw paired OBS/ACTM', color='black')
    ax[2, 1].plot(months, monthly['paired_fit_days'][:, site_i],
                  lw=0.85, label='Fitted paired OBS/ACTM', color='tab:blue')
    ax[2, 1].set_ylim(-0.5, 32)
    ax[2, 1].set_title('F. Data coverage per calendar month')
    ax[2, 1].set_ylabel('Paired days / month')
    ax[2, 1].legend(fontsize=8)

    for row_axes in ax:
        for a in row_axes:
            a.grid(alpha=0.2)
            if a in (ax[0, 0], ax[1, 0], ax[1, 1], ax[2, 0], ax[2, 1]):
                a.xaxis.set_major_locator(mdates.YearLocator(5))
                a.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
                a.set_xlim(days[0], days[-1])
    figure.suptitle(hcase + ' | ' + site + ' | ID ' + str(station_id) +
                   ' | Lat ' + f'{row.lat:.2f}' + ', Lon ' + f'{row.lon:.2f}', fontsize=14)
    filename = outdir / (f'station_{station_id:03d}_' + safe + '.png')
    figure.savefig(filename, dpi=140)
    plt.close(figure)
    return filename


# ============================================================
# === One human-readable TXT table for active stations
#     All statistics are derived, never used to modify F6/F7 inputs.
def write_station_stats(outdir, sites, selected, f6_stats=None, phase=None, f7_files=None,
                        monthly=None, months=None, inversion_checks=None, duplicates=None,
                        identical_series=None, inv_period=None):
    f7_files = f7_files or {}
    inversion_checks = inversion_checks or {}
    duplicates = duplicates or []
    identical_series = identical_series or []
    primary = get_primary_f7(f7_files)
    primary_name = Path(primary['filename']).name if primary is not None else None
    inv_check = inversion_checks.get(primary_name) if primary_name else None
    columns = (' ID  Site                                    Lat      Lon  Nobs  Nfit  NF7  '
               'Nmiss NoRaw FirstF7 LastF7   InvBad BiasF7  RMSEF7 SigmaF7 AmpOBS AmpACTM Lag')
    lines = ['CH4 F6/F7 DIAGNOSTICS',
             'F6 period: 1999-2025; F7: monthly; IDs from obsrvCH4_60sites.txt',
             'BiasF7 = mean(OBS-ACTM); RMSEF7 = sqrt(mean((OBS-ACTM)^2)); units = ppb',
             'Lag > 0 means ACTM seasonal cycle occurs later than OBS',
             'Nobs = raw OBS days; Nfit = paired fitted days; NF7 = valid F7 months',
             'Nmiss = invalid F7 months (out of 324); NoRaw = valid F7 months with zero raw paired days',
             'InvBad = invalid requested months in inversion period (using read_globalview last-name-match)',
             'F7 primary file: ' + (primary_name or '(none)'),
             columns, '-' * len(columns)]
    for i in selected:
        site = sites.iloc[i]
        rec = f6_stats.iloc[i] if f6_stats is not None else None
        ph = phase.iloc[i] if phase is not None else None
        nobs = int(rec.n_input_obs) if rec is not None else 0
        nfit = int(rec.n_fit_pairs) if rec is not None else 0
        nvalid, nmiss, no_raw = 0, np.nan, np.nan
        first_f7 = last_f7 = '-'
        bias = rmse = sigma = np.nan
        if primary is not None:
            values = primary['data'][:, i, :]
            mask = get_f7_valid(primary['data'])[:, i]
            nvalid = int(mask.sum())
            nmiss = int((~mask).sum())
            if nvalid:
                diff = values[mask, 0]
                bias = float(np.mean(diff))
                rmse = float(np.sqrt(np.mean(diff * diff)))
                sigma = float(np.mean(values[mask, 1]))
                if months is not None:
                    idx = np.flatnonzero(mask)
                    first_f7 = months[idx[0]].strftime('%Y-%m')
                    last_f7 = months[idx[-1]].strftime('%Y-%m')
                if monthly is not None:
                    no_raw = int(np.sum(mask & (monthly['paired_raw_days'][:, i] == 0)))
        amp_obs = ph.obs_seasonal_amplitude_ppb if ph is not None else np.nan
        amp_mod = ph.actm_seasonal_amplitude_ppb if ph is not None else np.nan
        lag = ph.seasonal_lag_months if ph is not None else np.nan
        inv_bad = (inv_check['by_id'].get(int(site.station_id), {}).get('missing')
                   if inv_check is not None else None)
        def num(x, n=7):
            return f'{x:{n}.1f}' if np.isfinite(x) else ' ' * (n - 3) + 'N/A'
        lag_string = f'{int(lag):3d}' if np.isfinite(lag) else 'N/A'
        miss_string = f'{nmiss:5d}' if np.isfinite(nmiss) else '  N/A'
        raw_string = f'{no_raw:5d}' if np.isfinite(no_raw) else '  N/A'
        inv_string = f'{inv_bad:6d}' if inv_bad is not None else '   N/A'
        lines.append(f'{int(site.station_id):3d}  {str(site.site)[:39]:39s} '
                     f'{float(site.lat):7.2f} {float(site.lon):8.2f} '
                     f'{nobs:5d} {nfit:5d} {nvalid:4d} {miss_string} {raw_string} '
                     f'{first_f7:7s} {last_f7:7s} {inv_string} '
                     f'{num(bias)} {num(rmse)} {num(sigma)} {num(amp_obs)} {num(amp_mod)} {lag_string:>3s}')

    lines += ['', 'SEASONAL LAG CORRELATIONS (Nakazawa fitted minus trend)',
              'Cyclic 12-month climatology, shift ACTM relative to OBS; lag positive = ACTM later',
              'Delta = best r minus second-best r; WEAK if Delta < ' + str(LAG_WEAK_MARGIN) +
              ' (heuristic only; not statistical significance)',
              ' ID   r(-3)  r(-2)  r(-1)   r(0)  r(+1)  r(+2)  r(+3)  Best  Delta  QC']
    for i in selected:
        sid = int(sites.iloc[i].station_id)
        if phase is None:
            continue
        row = phase.iloc[i]
        vals = [row['corr_' + str(lag)] for lag in range(-3, 4)]
        lag = row.seasonal_lag_months
        delta = row.lag_margin
        lag_str = f'{int(lag):+3d}' if np.isfinite(lag) else 'N/A'
        delta_str = f'{delta:5.3f}' if np.isfinite(delta) else '  N/A'
        corr_str = ' '.join(f'{v:6.3f}' if np.isfinite(v) else '   N/A' for v in vals)
        status = ('WEAK' if row.lag_weak else 'OK') if np.isfinite(lag) else 'N/A'
        lines.append(f'{sid:3d} {corr_str} {lag_str:>5s} {delta_str:>6s}  {status}')

    lines += ['', 'DUPLICATE ACTIVE STATION NAMES (read_globalview uses LAST NetCDF name match)']
    if duplicates:
        for group in duplicates:
            equal = ('YES' if group['rows_identical'] else 'NO') if group['rows_identical'] is not None else 'N/A'
            lines.append('IDs=' + ','.join(map(str, group['ids'])) + ' | ' + group['site'] +
                         ' | F7 rows identical=' + equal + ' | matched ID=' + str(group['inversion_last_id']))
    else:
        lines.append('None among active stations')

    lines += ['', 'IDENTICAL ACTIVE F7 SERIES (exact monthly OBS-ACTM; station IDs not names)']
    if identical_series:
        for pair in identical_series:
            lines.append('IDs=' + str(pair['id_a']) + ',' + str(pair['id_b']) +
                         ' | Nvalid=' + str(pair['n_valid']) +
                         ' | sigma equal=' + ('YES' if pair['same_uncertainty'] else 'NO') +
                         ' | names equal=' + ('YES' if pair['same_name'] else 'NO'))
            lines.append('    ' + pair['name_a'] + ' | ' + pair['name_b'])
    else:
        lines.append('No fully identical active difference series with identical validity masks')

    lines += ['', 'INVERSION INPUT VALIDITY (all active requested stations, each NetCDF variant)']
    if inv_period:
        lines.append('Period: ' + inv_period[0] + '; source: ' + inv_period[1])
    if inversion_checks:
        for filename, info in inversion_checks.items():
            lines.append(filename + ': invalid=' + str(info['n_invalid']) +
                         ', literal -999=' + str(info['n_fill']) +
                         ', station-month requests=' + str(info['n_requested'] * info['n_months']))
            for sid, site, date, count in info['examples']:
                lines.append('    ID=' + str(sid) + ' ' + date + ' | invalid months=' + str(count) + ' | ' + site)
    else:
        lines.append('No F7 NetCDF checked')
    lines.append('WARNING: read_globalview.f does not mask -999; do not run inversion with invalid active months.')
    filename = outdir / 'ch4_station_stats.txt'
    filename.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return filename


# ============================================================
# === Run read-only diagnostics for one CYC/INCA case
#     Default: plot every active station; --stations START END limits IDs.
def run_diagnostics(hcase, setup, outdir=None, plots=True, stations=None, check_f7=True,
                    skip_f6=False, skip_f7=False, endian='<', control_file=None):
    hcase = hcase.upper()
    if hcase not in ['CYC', 'INCA']:
        raise ValueError('Expected CYC or INCA')
    if skip_f6 and skip_f7:
        raise ValueError('Cannot skip both F6 and F7')
    sites = get_sites(setup.sites_f)
    firstdat, lastdat, inv_source = get_inversion_period(setup, control_file)
    selected = get_plot_indices(sites, stations)
    days, months = get_dates()
    outdir = Path(outdir) if outdir else Path(setup.inv_wrk_dir) / 'diagnostics'
    outdir = outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print('\n\tCase          :', hcase)
    print('\tActive sites  :', int(sites.active.sum()))
    print('\tSelected sites:', len(selected), 'in', stations if stations else 'all active')
    print('\tOutput        :', outdir)
    if not selected:
        print('\tWARNING: no active station in selected range')

    inp = fits = monthly = f6_stats = phase = None
    if not skip_f6:
        print('\n\tF6 input  :', Path(setup.f6_input_file).resolve())
        print('\tF6 output :', Path(setup.f6_output_file).resolve())
        inp, fits = get_f6_files(setup.f6_input_file, setup.f6_output_file, endian)
        monthly = get_f6_monthly(fits, days, months)
        f6_stats = get_f6_stats(inp, fits, sites)
        phase = get_f6_phase(fits, sites, days)
        nbad = int(f6_stats.n_raw_preservation_errors.sum())
        print('\tF6 raw preservation mismatches:', nbad)
        nh_lag = phase.loc[phase.active & (phase.lat >= 50) & (phase.seasonal_lag_months == 1)]
        print('\tHigh-latitude NH active ACTM +1 mo lags:', len(nh_lag))

    f7_files = {}
    mismatch_examples = []
    mismatch_total = 0
    inv_checks = {}
    if not skip_f7:
        prefix = str(Path(setup.f7_output_nc_prefix).expanduser().resolve())
        matches = sorted(Path(prefix).parent.glob(Path(prefix).name + '*.nc'))
        if not matches:
            print('\tWARNING: no F7 NetCDF found for', prefix)
        for filename in matches:
            if hcase.lower() not in filename.name.lower():
                continue
            nc = get_f7_file(filename)
            stat = get_f7_stats(nc, sites)
            count_missing = int(stat.loc[stat.active, 'n_missing'].sum())
            count_high = int(stat.loc[stat.active, 'n_uncertainty_gt20'].sum())
            print('\tF7', filename.name, '| missing:', count_missing, '| sigma>20:', count_high)
            if check_f7 and monthly is not None:
                unc_tag = filename.name[len(Path(prefix).name):-3][:3]
                if unc_tag in ('ctl', 'ux2', 'ux4'):
                    expected = get_expected_f7(monthly, {'ctl': 1, 'ux2': 2, 'ux4': 4}[unc_tag], months)
                    result, examples = get_comparison(nc, expected, sites, months)
                    print('\tF7 vs F6 mismatch months:', result['n_total_bad_months'], '(' + unc_tag + ')')
                    mismatch_total += result['n_total_bad_months']
                    mismatch_examples.extend(examples)
                else:
                    print('\tWARNING: F7 crosscheck skipped for unknown uncertainty:', filename.name)
            inv_checks[filename.name] = get_inversion_missing(nc, sites, months, firstdat, lastdat)
            print('\tINV missing requests:', inv_checks[filename.name]['n_invalid'], '| literal -999:', inv_checks[filename.name]['n_fill'])
            f7_files[filename.name] = nc

    primary = get_primary_f7(f7_files)
    duplicates = get_duplicate_active_sites(sites, primary)
    identical_series = get_identical_active_series(sites, primary)
    if duplicates:
        print('\tWARNING: repeated ACTIVE station names:', len(duplicates))
        for group in duplicates:
            print('\t ', group['site'], '| IDs:', group['ids'], '| last match:', group['inversion_last_id'], '| F7 identical:', group['rows_identical'])
    if identical_series:
        print('\tWARNING: identical ACTIVE F7 difference series:', len(identical_series))
        for pair in identical_series:
            print('\t ', pair['id_a'], pair['id_b'], '| sigma same:', pair['same_uncertainty'])
    print('\tInversion data interval:', firstdat, '-', lastdat, '|', inv_source)
    if not inv_source.startswith('control:'):
        print('\tWARNING: inversion interval not verified from an actual control file')
    if any(info['n_invalid'] for info in inv_checks.values()):
        print('\tWARNING: invalid active station-months would be read by legacy read_globalview.f')
    if primary is not None and monthly is not None:
        valid = get_f7_valid(primary['data'])
        rawless = valid & (monthly['paired_raw_days'] == 0)
        n_rawless = int(rawless[:, sites.active.to_numpy()].sum())
        print('\tValid F7 active months without paired raw days:', n_rawless)
    if phase is not None:
        weak = phase.loc[phase.active & phase.lag_weak]
        print('\tActive stations with weak lag separation (<0.05):', len(weak))

    stats_file = write_station_stats(outdir, sites, selected, f6_stats, phase, f7_files,
                                     monthly=monthly, months=months, inversion_checks=inv_checks,
                                     duplicates=duplicates, identical_series=identical_series,
                                     inv_period=(str(firstdat) + '-' + str(lastdat), inv_source))
    print('\n\tStation statistics:', stats_file)
    if mismatch_total:
        print('\tWARNING: F6 -> F7 inconsistencies (first 10):')
        for example in mismatch_examples[:10]:
            print('\t ', example['file'], example['month'], example['site'],
                  'expected=', example['expected_difference'], 'actual=', example['actual_difference'])

    if plots and fits is not None:
        for i in selected:
            path = plot_station(outdir, i, sites, days, months, inp, fits, monthly, f7_files, hcase)
            print('\tPlot:', path.name)
    elif plots and fits is None:
        print('\tNo station plots: --skip-f6 was specified')
    print('\n\tDiagnostics saved:', outdir)
    return {'outdir': outdir, 'station_ids': [int(sites.iloc[i].station_id) for i in selected],
            'comparison_errors': mismatch_total, 'n_identical_active_series': len(identical_series),
            'statistics': stats_file}


# ============================================================
# === CLI
def main():
    parser = argparse.ArgumentParser(description='Read-only diagnostics for F6 binary and F7 NetCDF files')
    parser.add_argument('--case', choices=['CYC', 'INCA', 'ALL'], default='CYC')
    parser.add_argument('--stations', nargs=2, type=int, metavar=('START', 'END'), default=None, help='Inclusive station ID range; default all active stations')
    parser.add_argument('--output', default=None, help='Custom report directory, if one case')
    parser.add_argument('--no-plots', action='store_true')
    parser.add_argument('--no-crosscheck', action='store_true', help='Skip F7 reconstruction from F6')
    parser.add_argument('--skip-f6', action='store_true', help='Read only F7 NetCDF')
    parser.add_argument('--skip-f7', action='store_true', help='Read only F6 binary')
    parser.add_argument('--control', default=None, help='Actual inversion control file (firstdat/lastdat)')
    parser.add_argument('--endian', choices=['little', 'big'], default='little')
    args = parser.parse_args()
    if args.skip_f6 and args.skip_f7:
        parser.error('Cannot skip both F6 and F7')
    if args.case == 'ALL' and args.output:
        parser.error('--output with ALL would mix CYC and INCA reports; leave it unset')
    import _set_case
    for hcase in (['CYC', 'INCA'] if args.case == 'ALL' else [args.case]):
        os.environ['CH4_CASE'] = hcase
        setup = _set_case.SetCase()
        run_diagnostics(hcase, setup, outdir=args.output, plots=not args.no_plots,
                        stations=args.stations, check_f7=not args.no_crosscheck,
                        skip_f6=args.skip_f6, skip_f7=args.skip_f7,
                        endian='<' if args.endian == 'little' else '>',
                        control_file=args.control)


if __name__ == '__main__':
    main()
