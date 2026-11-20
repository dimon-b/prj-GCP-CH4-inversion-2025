# -*- coding: utf-8 -*-
"""D3: ACTM stations -> OBS/ACTM -> Nakazawa -> TDI observations.

All Python orchestration is implemented here. Fortran executables read
D3 operation-specific namelists; no stage_impl/helper imports.
"""
import os
import subprocess
import shutil
import tempfile
from datetime import datetime
from pathlib import Path
import calendar
import numpy as np
import pandas as pd
from scipy.io import netcdf_file

from _set_case import SetCase


# ============================================================
# === D3: observations and ACTM concentrations
class ObsACTM(SetCase):

    def __init__(self, case=None):
        super().__init__()
        self.case = self.hcase.upper()
        if case is not None and self.case != str(case).upper():
            raise ValueError('SetCase mismatch: ' + self.case + ' != ' + str(case))
        self.root = Path(__file__).resolve().parent
        self.config_dir = (Path(self.inv_wrk_dir) / 'config').resolve()
        self.bin_dir = Path(self.fort_exe_dir).resolve()

    # --------------------------------------------------------
    # --- D3 namelist names must match their Fortran group names
    NAMES = {
        'sample': 'd3_sample_actm_sites',
        'combine': 'd3_combine_obs_actm',
        'filter': 'd3_filter_ch4_nakazawa',
        'observations': 'd3_write_tdi_observations',
    }

    def write_config(self, operation, fields, config_path=None):
        name = self.NAMES[operation]
        for key, val in fields.items():
            if "'" in str(val) or '\n' in str(val):
                raise ValueError('Invalid namelist field: ' + key)
        content = '&' + name + '\n'
        content += " hcase = '" + self.case + "'\n"
        for key, val in fields.items():
            if isinstance(val, int):
                content += ' ' + key + ' = ' + str(val) + '\n'
            else:
                content += " " + key + " = '" + str(val) + "'\n"
        content += '/\n'
        self.config_dir.mkdir(parents=True, exist_ok=True)
        dest = Path(config_path) if config_path is not None else self.config_dir / (name + '.nml')
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(content, encoding='utf-8')
        print('D3', operation, 'namelist:', dest)
        return dest

    # --- Compile D3 Fortran when binary is missing or sources have changed
    def executable(self, operation):
        name = 'df3_' + self.NAMES[operation].removeprefix('d3_')
        configured = getattr(self, 'exe_d3_' + operation, '')
        exe = Path(configured).expanduser().resolve() if configured else self.bin_dir / name
        if exe.name != name:
            raise ValueError('D3 expects ' + name + ', got ' + exe.name)

        source = self.root / 'fortran' / (name + '.f90')
        if not source.is_file():
            raise FileNotFoundError('Missing D3 Fortran source: ' + str(source))

        sources = [source]
        if operation == 'filter':
            sources += [self.root / 'fortran' / 'lsfbf4_main.f',
                        self.root / 'fortran' / 'lsfbf4_sub.f']
        for item in sources:
            if not item.is_file():
                raise FileNotFoundError('Missing Fortran dependency: ' + str(item))

        rebuild = not exe.is_file()
        if not rebuild:
            rebuild = any(p.stat().st_mtime_ns > exe.stat().st_mtime_ns for p in sources)
        if rebuild:
            self.build_executable(operation, sources, exe)
        if not exe.is_file():
            raise FileNotFoundError('D3 compilation did not create executable: ' + str(exe))
        return exe

    def build_executable(self, operation, sources, exe):
        compiler = (os.environ.get('FC') or shutil.which('ifort') or
                    shutil.which('ifx') or shutil.which('gfortran'))
        if not compiler:
            raise RuntimeError('No Fortran compiler: install ifort, ifx, or gfortran, or set FC')
        compiler_name = Path(compiler).name.lower()
        cmd = [compiler, '-O2']
        if 'gfortran' in compiler_name:
            cmd += ['-std=legacy', '-ffree-line-length-none',
                    '-ffixed-line-length-none', '-fallow-argument-mismatch']
        if operation == 'observations':
            nf_config = shutil.which('nf-config')
            if not nf_config:
                raise RuntimeError('df3_write_tdi_observations requires nf-config (NetCDF-Fortran)')
            flags = subprocess.check_output([nf_config, '--fflags'], text=True).split()
            libs = subprocess.check_output([nf_config, '--flibs'], text=True).split()
        else:
            flags, libs = [], []

        exe.parent.mkdir(parents=True, exist_ok=True)
        cmd += flags + [str(p) for p in sources] + ['-o', str(exe)] + libs
        print('Build D3', operation, ':', ' '.join(cmd), flush=True)
        subprocess.run(cmd, cwd=str(exe.parent), check=True)

    def run_exe(self, operation, config, work_dir=None):
        exe = self.executable(operation)
        print('Run D3', operation, exe)
        if work_dir is not None:
            work_dir.mkdir(parents=True, exist_ok=True)
        subprocess.run([str(exe), str(config)], cwd=work_dir, check=True)

    @staticmethod
    def require_file(path):
        p = Path(path).expanduser().resolve()
        if not p.is_file() or p.stat().st_size == 0:
            raise FileNotFoundError('Missing/empty D3 file: ' + str(p))
        return p

    # --------------------------------------------------------
    # --- D3 sample: ACTM concentrations at station locations
    def get_actm_at_sites(self):
        model = Path(self.d3_sample_model_dir).expanduser().resolve()
        sites = self.require_file(self.sites_f)
        if not model.is_dir() or self.case not in model.name.upper():
            raise FileNotFoundError('Wrong/missing ACTM directory for ' + self.case + ': ' + str(model))
        out = Path(self.d3_sample_output_dir).expanduser().resolve()
        out.mkdir(parents=True, exist_ok=True)
        cfg = self.write_config('sample', {
            'obs_file': str(sites), 'i_dir': str(model) + os.sep,
            'o_dir': str(out) + os.sep, 'o_fil': self.d3_sample_output_prefix,
            's_year': int(self.yr_s), 'f_year': int(self.yr_e),
        })
        self.run_exe('sample', cfg)
        for year in range(self.yr_s, self.yr_e + 1):
            self.require_file(out / (self.d3_sample_output_prefix + str(year) + '.grd'))
        print('D3 sample PASS:', self.yr_e - self.yr_s + 1, 'station annual files')

    # --------------------------------------------------------
    # --- D3 combine: join station observations and model concentrations
    def combine_obs_actm(self):
        model = Path(self.d3_combine_model_dir).expanduser().resolve()
        obs = Path(self.d3_combine_obs_dir).expanduser().resolve()
        out = Path(self.d3_combine_output_dir).expanduser().resolve()
        sites = self.require_file(self.sites_f)
        if not obs.is_dir():
            raise FileNotFoundError(obs)
        for year in range(self.yr_s, self.yr_e + 1):
            self.require_file(model / (self.d3_combine_input_prefix + str(year) + '.grd'))
        out.mkdir(parents=True, exist_ok=True)
        cfg = self.write_config('combine', {
            'obs_file': str(sites), 'obs_dir': str(obs) + os.sep,
            'i_dir': str(model) + os.sep, 'o_dir': str(out) + os.sep,
            'i_fil': self.d3_combine_input_prefix,
            'o_yr': self.d3_combine_output_prefix + 'ch4yr.grd',
            'o_mm': self.d3_combine_output_prefix + 'ch4mm.grd',
            'o_dy': self.d3_combine_output_prefix + 'ch4dy.grd',
            's_year': int(self.yr_s), 'f_year': int(self.yr_e),
        })
        self.run_exe('combine', cfg)
        for suffix in ('ch4yr.grd', 'ch4mm.grd', 'ch4dy.grd'):
            self.require_file(out / (self.d3_combine_output_prefix + suffix))
        print('D3 combine PASS')

    # --------------------------------------------------------
    # --- Commit validated outputs, preserving the previous version
    def commit_outputs(self, staged, targets, operation):
        staged = [Path(p) for p in staged]
        targets = [Path(p) for p in targets]
        if len(staged) != len(targets) or not staged:
            raise ValueError('D3 output lists are inconsistent')
        if len(set(targets)) != len(targets):
            raise ValueError('D3 output paths are not unique')
        for path in staged:
            self.require_file(path)

        existing = [path for path in targets if path.is_file()]
        backup_dir = None
        backups = {}
        if existing:
            stamp = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
            backup_dir = Path(self.inv_wrk_dir).resolve() / 'backups' / (operation + '_' + stamp)
            backup_dir.mkdir(parents=True, exist_ok=False)
            for index, path in enumerate(targets):
                if path.is_file():
                    backup = backup_dir / (str(index).zfill(2) + '_' + path.name)
                    shutil.copy2(path, backup)
                    backups[path] = backup

        completed = []
        try:
            for source, target in zip(staged, targets):
                target.parent.mkdir(parents=True, exist_ok=True)
                source.replace(target)
                completed.append(target)
        except Exception:
            for target in completed:
                if target in backups:
                    shutil.copy2(backups[target], target)
                elif target.exists():
                    target.unlink()
            raise

        if backup_dir is not None:
            print('D3', operation, 'backup:', backup_dir)

    # --------------------------------------------------------
    # --- F6 output layout: 9862 direct-access records, 245 stations, 2 tracers,
    #     3 float32 fields per station/tracer (see df3_filter_ch4_nakazawa.f90)
    def check_filter_output(self, path):
        output = self.require_file(path)
        expected = 9862 * 245 * 2 * 3 * 4
        actual = output.stat().st_size
        if actual != expected:
            raise ValueError('D3 filter output size mismatch: ' + str(actual) +
                             ' != ' + str(expected) + ' bytes: ' + str(output))
        return output

    # --------------------------------------------------------
    # --- D3 filter: Nakazawa digital filter
    def filter_ch4(self):
        inp = self.require_file(self.d3_filter_input_file)
        out = Path(self.d3_filter_output_file).expanduser().resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        fields = {'i_file': str(inp), 'o_file': str(out)}
        self.write_config('filter', fields)

        # --- Write to a temporary location; old output stays intact if Fortran fails
        with tempfile.TemporaryDirectory(prefix='.d3_filter_', dir=out.parent) as tmp:
            staged = Path(tmp) / out.name
            run_fields = {'i_file': str(inp), 'o_file': str(staged)}
            cfg = self.write_config('filter', run_fields, Path(tmp) / 'd3_filter_ch4_nakazawa.nml')
            work = (Path(self.inv_wrk_dir) / 'diagnostics' / 'd3' / 'filter_work').resolve()
            self.run_exe('filter', cfg, work)
            self.check_filter_output(staged)
            self.commit_outputs([staged], [out], 'filter')

        print('D3 filter PASS')

    # --------------------------------------------------------
    # --- Expected outputs of the Fortran TDI writer
    def observation_products(self, summary, prior, prefix):
        products = [summary, prior]
        for unc in self.unxcases:
            if unc != 'ctl':
                products += [Path(str(summary) + '.' + unc), Path(str(prior) + '.' + unc)]
        for unc in self.unxcases:
            for inv in self.invcases:
                products.append(Path(str(prefix) + unc + inv + '.nc'))
        return products

    # --------------------------------------------------------
    # --- D3 observations: TDI inversion observations
    def write_tdi_observations(self):
        sites = self.require_file(self.sites_f)
        fitted = self.require_file(self.d3_observations_input_file)
        summary = Path(self.d3_observations_summary).expanduser().resolve()
        prior = Path(self.d3_observations_prior).expanduser().resolve()
        prefix = Path(self.d3_observations_nc_prefix).expanduser().resolve()
        targets = self.observation_products(summary, prior, prefix)

        # --- The temporary outputs must be on the same filesystem as final outputs
        if len({p.parent for p in targets}) != 1:
            raise ValueError('D3 observation outputs must share one directory')
        summary.parent.mkdir(parents=True, exist_ok=True)
        fields = {
            'obs_file': str(sites), 'o_obsrv': str(summary),
            'i_infit': str(fitted), 'o_prior': str(prior),
            'o_pr_nc': str(prefix),
            'cfg_s_year': int(self.yr_s), 'cfg_f_year': int(self.yr_e),
        }
        self.write_config('observations', fields)

        with tempfile.TemporaryDirectory(prefix='.d3_observations_', dir=summary.parent) as tmp:
            staged_summary = Path(tmp) / summary.name
            staged_prior = Path(tmp) / prior.name
            staged_prefix = Path(tmp) / prefix.name
            run_fields = dict(fields)
            run_fields['o_obsrv'] = str(staged_summary)
            run_fields['o_prior'] = str(staged_prior)
            run_fields['o_pr_nc'] = str(staged_prefix)
            cfg = self.write_config('observations', run_fields, Path(tmp) / 'd3_write_tdi_observations.nml')
            self.run_exe('observations', cfg)
            staged = self.observation_products(staged_summary, staged_prior, staged_prefix)
            self.commit_outputs(staged, targets, 'observations')

        print('D3 observations PASS:', len(targets), 'outputs')

    # --------------------------------------------------------
    # --- Diagnostics after completing D3, or on existing results
    def diagnostics(self, stations=None, plots=True):
        if (self.yr_s, self.yr_e) != (YEAR_S, YEAR_E):
            raise ValueError('D3 diagnostics expects 1999-2025')
        out = Path(self.d3_diagnostics_dir).expanduser().resolve() if hasattr(self, 'd3_diagnostics_dir') else Path(self.inv_wrk_dir).resolve() / 'diagnostics' / 'd3'
        return generate_d3_diagnostics(self.case, self, outdir=out, plots=plots, stations=stations)

    # --------------------------------------------------------
    # --- Complete D3
    def run(self):
        self.get_actm_at_sites()
        self.combine_obs_actm()
        self.filter_ch4()
        self.write_tdi_observations()
        self.diagnostics()
        print('\nD3 OBS + ACTM: PASS')


# ============================================================
# === Called from _main.py

def run_stage(case='INCA', dry_run=False):
    case = str(case).upper()
    if case not in ('INCA', 'CYC'):
        raise ValueError('Unknown case: ' + case)
    old_gcp, old_ch4 = os.environ.get('GCP_CASE'), os.environ.get('CH4_CASE')
    os.environ['GCP_CASE'] = case
    os.environ['CH4_CASE'] = case
    try:
        obj = ObsACTM(case=case)
        if dry_run:
            print('D3 dry-run:', case, obj.yr_s, obj.yr_e)
            print('sample -> combine -> filter -> observations -> diagnostics; no writes or Fortran calls')
            return
        obj.run()
    finally:
        for key, value in [('GCP_CASE', old_gcp), ('CH4_CASE', old_ch4)]:
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


# ============================================================
# === D3 scientific diagnostics (integrated from the original c4 diagnostics)
# ============================================================
# === Constants from Fortran (do not silently infer a different format)
N_SITE = 245
YEAR_S = 1999
YEAR_E = 2025
N_DAY = sum(366 if calendar.isleap(y) else 365 for y in range(YEAR_S, YEAR_E + 1))
N_MONTH = (YEAR_E - YEAR_S + 1) * 12
FILL = -999.0


# ============================================================
# === Dates and site metadata
def get_dates():
    days = pd.date_range(str(YEAR_S) + '-01-01', str(YEAR_E) + '-12-31', freq='D')
    months = pd.date_range(str(YEAR_S) + '-01-01', str(YEAR_E) + '-12-01', freq='MS')
    if len(days) != N_DAY or len(months) != N_MONTH:
        raise ValueError('D3 calendar does not match 1999-2025')
    return days, months


def get_sites(filename):
    # --- Fortran TDI skips the first two lines of obsrvCH4_60sites.txt.
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
#     Input filter:  record(day) = (ntrac=2, nsite=245) => 1960 bytes/day
#     Output filter: record(day) = (ntrac=2, kind=3, nsite=245) => 5880 bytes/day

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


def get_filter_files(input_file, output_file, endian='<'):
    inp = get_binary(input_file, (N_DAY, 2, N_SITE), endian)
    out = get_binary(output_file, (N_DAY, 2, 3, N_SITE), endian)
    return inp, out


# ============================================================
# === TDI NetCDF classic. The Python axis order is (month, station, type).
def get_name(row):
    if isinstance(row, bytes):
        return row.decode('ascii', errors='replace').strip(' \x00')
    return b''.join(row).decode('ascii', errors='replace').strip(' \x00')


def get_tdi_netcdf(filename):
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
        raise ValueError(str(filename) + ': cannot read as NetCDF classic. TDI uses nf_create(NF_CLOBBER). ' + str(exc)) from exc
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


def get_tdi_valid(data):
    return np.isfinite(data[:, :, 0]) & np.isfinite(data[:, :, 1]) & (data[:, :, 0] != FILL) & (data[:, :, 1] > 0.0)


# ============================================================
# === Monthly reconstruction using EXACT TDI pairing and RMS rules
#     TDI uses fitted OBS and ACTM on common valid days for the difference.
#     Monthly RMS(raw OBS - fitted OBS) is pooled over the same calendar
#     month across all 27 years, then uncertainty = RMS * scale + 5 ppb.
#     If raw monthly observation is unavailable, uncertainty doubles.

def get_filter_monthly(fits, days=None, months=None):
    if days is None or months is None:
        days, months = get_dates()
    if fits.shape != (len(days), 2, 3, N_SITE):
        raise ValueError('filter fits shape must be (days, 2, 3, 245)')

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
            raise ValueError('filter daily file ends before ' + str(date))
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
        # --- TDI calculates stdev using OBS raw+fit, regardless of ACTM.
        residual_mask = get_valid(obs_raw) & get_valid(obs_fit)
        residual = np.where(residual_mask, obs_raw - obs_fit, 0.0)
        sumsq[date.month - 1] += np.square(residual).sum(axis=0)
        nsum[date.month - 1] += residual_mask.sum(axis=0)
        offset += count

    if offset != len(days):
        raise ValueError('filter calendar mismatch: ' + str(offset) + ' / ' + str(len(days)))
    rms = np.sqrt(np.divide(sumsq, nsum, out=np.full_like(sumsq, np.nan), where=nsum > 0))
    return {'fit_diff': fit_diff, 'raw_diff': raw_diff, 'obs_raw_mean': obs_raw_mean,
            'obs_fit_mean': obs_fit_mean, 'mod_fit_mean': mod_fit_mean, 'rms': rms,
            'paired_fit_days': paired_fit_days, 'paired_raw_days': paired_raw_days}


def get_expected_tdi(monthly, scale, months=None):
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
def get_filter_phase(fits, sites, days=None, min_days_per_calendar_month=60):
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
        lag, corr = get_best_lag(a, b)
        records.append({'station_id': int(row.station_id), 'site': row.site, 'active': bool(row.active),
                        'lat': float(row.lat), 'lon': float(row.lon),
                        'n_calendar_months': int(np.sum(np.isfinite(a) & np.isfinite(b))),
                        'min_days_per_month': int(n_days[:, i].min()),
                        'seasonal_lag_months': lag, 'best_correlation': corr,
                        'obs_seasonal_amplitude_ppb': (float(np.max(a) - np.min(a)) if np.isfinite(a).all() else np.nan),
                        'actm_seasonal_amplitude_ppb': (float(np.max(b) - np.min(b)) if np.isfinite(b).all() else np.nan)})
    return pd.DataFrame(records)


# ============================================================
# === Statistical helpers
def get_statistics(values):
    v = np.asarray(values, dtype='f8')
    v = v[np.isfinite(v)]
    if not len(v):
        return (np.nan, np.nan, np.nan, np.nan)
    return (float(np.mean(v)), float(np.std(v)), float(np.min(v)), float(np.max(v)))


def get_best_lag(obs, mod):
    # --- Positive months means ACTM seasonality occurs later than OBS.
    #     Compare 12-month climatologies; remove means, no fitted phase shift.
    if not np.isfinite(obs).all() or not np.isfinite(mod).all():
        return np.nan, np.nan
    a = obs - np.mean(obs)
    b = mod - np.mean(mod)
    if np.std(a) < 0.1 or np.std(b) < 0.1:
        return np.nan, np.nan
    corr = []
    for lag in range(-3, 4):
        c = np.corrcoef(a, np.roll(b, -lag))[0, 1]
        corr.append((lag, c))
    return sorted(corr, key=lambda p: p[1], reverse=True)[0]


def get_filter_stats(inp, fits, sites):
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


def get_tdi_stats(nc, sites):
    numbers = nc['numbers']
    names = nc['names']
    data = nc['data']
    if list(numbers) != list(sites.station_id) or names != list(sites.site):
        raise ValueError(nc['filename'] + ': NetCDF station order/names differ from site list')
    valid = get_tdi_valid(data)
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
                        'n_eligible_for_inversion': int(np.sum(v & (data[:, i, 1] <= 20)))})
    return pd.DataFrame(summary)


def get_tdi_monthly(nc, sites, months):
    data = nc['data']
    valid = get_tdi_valid(data)
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
    v_actual = get_tdi_valid(actual)
    v_expected = get_tdi_valid(expected)
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
#     Same paired-days approach and 60-day minimum as get_filter_phase().
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
def get_primary_tdi(files):
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
    primary = get_primary_tdi(files)
    rec = np.asarray(monthly['fit_diff'][:, site_i], dtype='f8')
    ax[1, 1].plot(months, rec, lw=1, color='tab:orange', label='From filter fitted values')
    if primary is not None:
        obs_minus_mod = np.asarray(primary['data'][:, site_i, 0], dtype='f8').copy()
        obs_minus_mod[~get_tdi_valid(primary['data'])[:, site_i]] = np.nan
        ax[1, 1].plot(months, obs_minus_mod, '.', ms=2.5, color='black', label='TDI NetCDF')
    ax[1, 1].axhline(0, color='gray', lw=0.7)
    ax[1, 1].set_title('D. Monthly OBS - ACTM (filter reconstruction vs TDI)')
    ax[1, 1].set_ylabel('OBS - ACTM (ppb)')
    ax[1, 1].legend(fontsize=8)

    # --- E: Observation uncertainty by TDI inversion option
    uncertainties = get_uncertainty_files(files)
    colors = {'ctl': 'tab:green', 'ux2': 'tab:orange', 'ux4': 'tab:red'}
    for unc, nc in uncertainties.items():
        sigma = np.asarray(nc['data'][:, site_i, 1], dtype='f8').copy()
        sigma[~get_tdi_valid(nc['data'])[:, site_i]] = np.nan
        ax[2, 0].plot(months, sigma, lw=0.9, color=colors[unc], label=unc)
    ax[2, 0].axhline(20, ls='--', lw=0.9, color='gray', label='20 ppb')
    ax[2, 0].set_title('E. Monthly TDI uncertainty')
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
#     All statistics are derived, never used to modify D3 inputs.
def write_station_stats(outdir, sites, selected, f6_stats=None, phase=None, f7_files=None):
    f7_files = f7_files or {}
    primary = get_primary_tdi(f7_files)
    columns = (' ID  Site                                    Lat      Lon  Nobs  Nfit  NTDI  '
               'BiasTDI  RMSETDI SigmaTDI AmpOBS AmpACTM Lag')
    lines = ['CH4 D3 DIAGNOSTICS',
             'filter period: 1999-2025; TDI: monthly; station IDs from obsrvCH4_60sites.txt',
             'BiasTDI = mean(OBS-ACTM); RMSETDI = sqrt(mean((OBS-ACTM)^2)); units = ppb',
             'Lag > 0 means ACTM seasonal cycle is delayed relative to OBS (months)',
             'Nobs = days with raw OBS; Nfit = days with paired fitted OBS/ACTM; NTDI = valid months',
             columns, '-' * len(columns)]
    for i in selected:
        site = sites.iloc[i]
        rec = f6_stats.iloc[i] if f6_stats is not None else None
        ph = phase.iloc[i] if phase is not None else None
        nobs = int(rec.n_input_obs) if rec is not None else 0
        nfit = int(rec.n_fit_pairs) if rec is not None else 0
        nvalid, bias, rmse, sigma = 0, np.nan, np.nan, np.nan
        if primary is not None:
            values = primary['data'][:, i, :]
            mask = get_tdi_valid(primary['data'])[:, i]
            nvalid = int(mask.sum())
            if nvalid:
                diff = values[mask, 0]
                bias = float(np.mean(diff))
                rmse = float(np.sqrt(np.mean(diff * diff)))
                sigma = float(np.mean(values[mask, 1]))
        amp_obs = ph.obs_seasonal_amplitude_ppb if ph is not None else np.nan
        amp_mod = ph.actm_seasonal_amplitude_ppb if ph is not None else np.nan
        lag = ph.seasonal_lag_months if ph is not None else np.nan
        def num(x, n=7):
            return f'{x:{n}.1f}' if np.isfinite(x) else ' ' * (n - 3) + 'N/A'
        lag_string = f'{int(lag):3d}' if np.isfinite(lag) else 'N/A'
        lines.append(f'{int(site.station_id):3d}  {str(site.site)[:39]:39s} '
                     f'{float(site.lat):7.2f} {float(site.lon):8.2f} '
                     f'{nobs:5d} {nfit:5d} {nvalid:4d} '
                     f'{num(bias)} {num(rmse)} {num(sigma)} {num(amp_obs)} {num(amp_mod)} {lag_string:>3s}')
    filename = outdir / 'ch4_station_stats.txt'
    filename.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return filename


# ============================================================
# === Run read-only diagnostics for one CYC/INCA case
#     Default: plot every active station; --stations START END limits IDs.
def generate_d3_diagnostics(hcase, setup, outdir=None, plots=True, stations=None, check_f7=True,
                    skip_f6=False, skip_f7=False, endian='<'):
    hcase = hcase.upper()
    if hcase not in ['CYC', 'INCA']:
        raise ValueError('Expected CYC or INCA')
    if skip_f6 and skip_f7:
        raise ValueError('Cannot skip both filter and TDI')
    sites = get_sites(setup.sites_f)
    selected = get_plot_indices(sites, stations)
    days, months = get_dates()
    outdir = Path(outdir) if outdir else Path(setup.inv_wrk_dir) / 'diagnostics' / 'd3'
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
        print('\n\tfilter input  :', Path(setup.d3_filter_input_file).resolve())
        print('\tfilter output :', Path(setup.d3_filter_output_file).resolve())
        inp, fits = get_filter_files(setup.d3_filter_input_file, setup.d3_filter_output_file, endian)
        monthly = get_filter_monthly(fits, days, months)
        f6_stats = get_filter_stats(inp, fits, sites)
        phase = get_filter_phase(fits, sites, days)
        nbad = int(f6_stats.n_raw_preservation_errors.sum())
        print('\tfilter raw preservation mismatches:', nbad)
        nh_lag = phase.loc[phase.active & (phase.lat >= 50) & (phase.seasonal_lag_months == 1)]
        print('\tHigh-latitude NH active ACTM +1 mo lags:', len(nh_lag))

    f7_files = {}
    mismatch_examples = []
    mismatch_total = 0
    if not skip_f7:
        prefix = str(Path(setup.d3_observations_nc_prefix).expanduser().resolve())
        matches = sorted(Path(prefix).parent.glob(Path(prefix).name + '*.nc'))
        if not matches:
            print('\tWARNING: no TDI NetCDF found for', prefix)
        for filename in matches:
            if hcase.lower() not in filename.name.lower():
                continue
            nc = get_tdi_netcdf(filename)
            stat = get_tdi_stats(nc, sites)
            count_missing = int(stat.loc[stat.active, 'n_missing'].sum())
            count_high = int(stat.loc[stat.active, 'n_uncertainty_gt20'].sum())
            print('\tTDI', filename.name, '| missing:', count_missing, '| sigma>20:', count_high)
            if check_f7 and monthly is not None:
                unc_tag = filename.name[len(Path(prefix).name):-3][:3]
                if unc_tag in ('ctl', 'ux2', 'ux4'):
                    expected = get_expected_tdi(monthly, {'ctl': 1, 'ux2': 2, 'ux4': 4}[unc_tag], months)
                    result, examples = get_comparison(nc, expected, sites, months)
                    print('\tTDI vs filter mismatch months:', result['n_total_bad_months'], '(' + unc_tag + ')')
                    mismatch_total += result['n_total_bad_months']
                    mismatch_examples.extend(examples)
                else:
                    print('\tWARNING: TDI crosscheck skipped for unknown uncertainty:', filename.name)
            f7_files[filename.name] = nc

    stats_file = write_station_stats(outdir, sites, selected, f6_stats, phase, f7_files)
    print('\n\tStation statistics:', stats_file)
    if mismatch_total:
        print('\tWARNING: Filter -> TDI inconsistencies (first 10):')
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
            'comparison_errors': mismatch_total, 'statistics': stats_file}




# ============================================================
# === Run only diagnostics on existing D3 outputs (no Fortran reprocessing)
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='D3 CH4 processing and diagnostics')
    parser.add_argument('--case', choices=['INCA', 'CYC'], default='INCA')
    parser.add_argument('--diagnostics-only', action='store_true')
    parser.add_argument('--stations', nargs=2, type=int, metavar=('START', 'END'))
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args()
    os.environ['GCP_CASE'] = args.case
    os.environ['CH4_CASE'] = args.case
    obj = ObsACTM(case=args.case)
    if args.diagnostics_only:
        obj.diagnostics(stations=args.stations, plots=not args.no_plots)
    else:
        if args.stations or args.no_plots:
            parser.error('--stations/--no-plots require --diagnostics-only')
        obj.run()
