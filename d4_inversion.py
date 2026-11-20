# -*- coding: utf-8 -*-
"""D4: existing TDI Bayesian inversion and independent read-only diagnostics.

Numerical Fortran stays in ../inv_dir/tdi_inv/. No stage_impl, no helpers.
The diagnostics write ONLY to results2025/<CASE>/diagnostics/d4/.
"""
import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from _set_case import SetCase


# ============================================================
# === D4: TDI inversion
class Inversion(SetCase):

    CONTROL_FILES = ('sourcef', 'datf', 'datflist', 'srccorf', 'outputf',
                     'groupf', 'constf', 'recipef', 'greenf')
    STATIC_FILES = ('srccorf', 'outputf', 'groupf', 'constf', 'recipef', 'greenf')
    EXPECTED_RESULTS = ('prior.sources', 'predicted.sources')

    def __init__(self, case=None):
        super().__init__()
        if case is not None and self.hcase != str(case).upper():
            raise ValueError('SetCase mismatch: ' + self.hcase + ' != ' + str(case))
        self.case = self.hcase
        self.tdi_dir = Path(self.d4_tdi_dir).resolve()
        self.executable = Path(self.exe_inv).resolve()
        self.template = Path(self.inv_control_template).resolve()
        self.files_dir = Path(self.d4_files4run_dir).resolve()
        self.gf_dir = Path(self.d4_gf_dir).resolve()
        self.result_root = Path(self.inv_wrk_dir).resolve()
        self.diag_dir = self.result_root / 'diagnostics' / 'd4'
        self.priors = list(self.d4_priors)
        self.uncs = list(self.d4_uncs)
        self.models = list(self.d4_models)
        for values, allowed in ((self.priors, self.unpcases),
                                (self.uncs, self.unxcases),
                                (self.models, self.invcases)):
            if not values or set(values) - set(allowed):
                raise ValueError('Invalid D4 selection: ' + str(values))

    # --------------------------------------------------------
    # --- Selection, paths, and immutable experiment names
    def cases(self):
        for prior in self.priors:
            for model in self.models:
                for unc in self.uncs:
                    yield prior, unc, model

    def get_paths(self, prior, unc, model):
        if prior not in self.unpcases or unc not in self.unxcases or model not in self.invcases:
            raise ValueError('Invalid inversion experiment: ' + str((prior, unc, model)))
        prefix_f5 = Path(self.d1_uncertainty_output_prefix)
        prefix_f7 = Path(self.d3_observations_nc_prefix)
        source = Path(str(prefix_f5) + prior.removeprefix('p') + '.nc').resolve()
        obs = Path(str(prefix_f7) + unc + model + '.nc').resolve()
        run_dir = self.result_root / (self.icase + '_' + prior)
        run_dir = run_dir / (self.icase + '_' + unc + '_' + model)
        return source, obs, run_dir

    @staticmethod
    def require_file(path):
        path = Path(path)
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError('Missing or empty input: ' + str(path))

    # --------------------------------------------------------
    # --- True D1/D3 NetCDF validation, not just file existence
    def check_netcdf(self, source, obs):
        import xarray as xr
        info = {}
        with xr.open_dataset(obs, decode_times=False, mask_and_scale=False) as ds:
            if 'tracer_data' not in ds:
                raise ValueError('D3 observations missing tracer_data: ' + str(obs))
            v = ds['tracer_data']
            for dim, expected in (('tlen', self.ntime), ('data_type', 2)):
                if dim not in v.dims or ds.sizes[dim] != expected:
                    raise ValueError('D3 observations incorrect ' + dim + ': ' + str(dict(ds.sizes)))
            if 'station_count' not in v.dims or ds.sizes['station_count'] < 1:
                raise ValueError('D3 observations station dimension is missing or empty')
            data = np.asarray(v.transpose('tlen', 'station_count', 'data_type').values, dtype=np.float64)
            values, sigma = data[..., 0], data[..., 1]
            valid = np.isfinite(values) & np.isfinite(sigma)
            valid &= (values > -900.0) & (sigma > 0.0) & (sigma < 1.0e20)
            nvalid = int(valid.sum())
            if not nvalid:
                raise ValueError('No valid D3 station-month observations: ' + str(obs))
            info['observations_stations'] = int(ds.sizes['station_count'])
            info['observations_months'] = self.ntime
            info['observations_valid_station_months'] = nvalid
            info['observations_invalid_station_months'] = int(valid.size - nvalid)
        with xr.open_dataset(source, decode_times=False, mask_and_scale=False) as ds:
            if 'sources_land' not in ds:
                raise ValueError('D1 prior missing sources_land: ' + str(source))
            a = ds['sources_land']
            if 54 not in a.shape:
                raise ValueError('D1 prior expected 54 land regions: ' + str(a.shape))
            if 'month' not in ds.sizes or ds.sizes['month'] < self.ntime:
                raise ValueError('D1 prior has insufficient monthly data: ' + str(dict(ds.sizes)))
            if 'time' not in ds or ds['time'].size != ds.sizes['month']:
                raise ValueError('D1 prior missing decimal-year time coordinate')
            times = np.asarray(ds['time'].values, dtype=np.float64)
            if not np.isfinite(times).all() or not (times[0] < self.yr_s and times[-1] > self.yr_e):
                raise ValueError('D1 prior does not cover inversion years 1999-2025')
            info['prior_land_regions'] = 54
            info['prior_months'] = int(ds.sizes['month'])
        return info

    # --------------------------------------------------------
    # --- Validate required TDI assets, without changing anything
    def check_inputs(self, prior, unc, model, dry_run=False):
        source, obs, run_dir = self.get_paths(prior, unc, model)
        assets = [self.template, self.files_dir, self.gf_dir,
                  Path(self.sites_f), source, obs]
        if not dry_run:
            assets.append(self.executable)
        for path in assets:
            if not path.exists():
                raise FileNotFoundError('D4 input missing: ' + str(path))
        if not self.files_dir.is_dir() or not self.gf_dir.is_dir():
            raise NotADirectoryError('files4run or gmats is not a directory')
        if not list(self.gf_dir.iterdir()):
            raise ValueError('Empty Green function directory: ' + str(self.gf_dir))
        # --- tc3greens.f90 constructs ACTM.ch4j55.YYYY.L2.nc by concatenation.
        gf_year_files = sorted(self.gf_dir.glob('ACTM.ch4j55.*.L2.nc'))
        if not gf_year_files:
            raise FileNotFoundError('No ACTM.ch4j55.YYYY.L2.nc Green functions in ' + str(self.gf_dir))
        if not dry_run and not os.access(self.executable, os.X_OK):
            raise PermissionError('TDI executable is not executable: ' + str(self.executable))
        self.require_file(self.template)
        self.require_file(Path(self.sites_f))
        self.require_file(source)
        self.require_file(obs)
        nc_info = self.check_netcdf(source, obs)
        print('D4 input:', self.case, prior, unc, model)
        print('  D1 prior       :', source)
        print('  D3 observations:', obs)
        print('  work     :', run_dir)
        print('  greenf   :', str(self.gf_dir) + os.sep)
        print('  GF years :', len(gf_year_files), 'NetCDF files')
        print('  valid observations:', nc_info['observations_valid_station_months'])
        return source, obs, run_dir, nc_info

    # --------------------------------------------------------
    # --- Parse and update only file paths in legacy TDI control
    @staticmethod
    def get_control_value(contents, key):
        matches = list(re.finditer(r'(?im)^\s*' + re.escape(key) +
                                   r'\s*=\s*([\'\"])([^\'\"\r\n]*)\1', contents))
        if len(matches) != 1:
            raise ValueError('Expected exactly one TDI control key: ' + key)
        return matches[0].group(2)

    @staticmethod
    def replace_control_value(contents, key, value):
        if "'" in str(value) or '\n' in str(value) or '\r' in str(value):
            raise ValueError('Invalid TDI filename: ' + str(value))
        if not value or len(str(value)) > 80:
            raise ValueError('TDI character*80 path limit exceeded: ' + key + ' ' + str(value))
        p = re.compile(r'(?im)^(\s*' + re.escape(key) + r'\s*=\s*)([\'\"])[^\'\"\r\n]*\2')
        contents, count = p.subn(lambda m: m.group(1) + "'" + str(value) + "'", contents)
        if count != 1:
            raise ValueError('Missing/duplicated TDI control key: ' + key)
        return contents

    @staticmethod
    def _check_switch(contents, key, value):
        m = re.search(r'(?im)^\s*' + key + r'\s*=\s*(\.?\w+\.?)', contents)
        if not m or m.group(1).strip('.').lower() not in value:
            raise ValueError('TDI control has incompatible setting: ' + key)

    def localize_control_file(self, run_dir, key, value):
        """Copy small static inputs; retain Green functions in shared directory.

        tc3greens.f90 concatenates greenf(1:lastchar) + ACTM.ch4j55.*.L2.nc.
        Therefore greenf MUST end with '/' (unlike Path.resolve()/relpath()).
        Original nonzero s057 TDI control used an absolute gmats/ prefix.
        """
        if key == 'greenf':
            gf_dir = Path(self.gf_dir).resolve()
            if not gf_dir.is_dir():
                raise NotADirectoryError('D4 shared GF directory missing: ' + str(gf_dir))
            prefix = str(gf_dir) + os.sep
            if len(prefix) > 80:
                raise ValueError('TDI greenf character*80 limit exceeded: ' + prefix)
            return prefix
        p = Path(value)
        if not p.is_absolute() and (run_dir / p).exists():
            return str(p)
        candidates = [p] if p.is_absolute() else [self.tdi_dir / p,
                                                   self.files_dir / p,
                                                   self.files_dir / p.name,
                                                   Path(self.inv_dir) / p]
        matches = [item.resolve() for item in candidates if item.exists()]
        if not matches:
            raise FileNotFoundError('TDI control ' + key + ' unresolved: ' + value)
        src = matches[0]
        try:
            relative = src.relative_to(self.files_dir)
            return str(Path('files4run') / relative)
        except ValueError:
            pass
        try:
            relative = src.relative_to(self.gf_dir)
            # --- Shared GF remains outside run_dir, no copy and no symlink.
            return os.path.relpath(src, start=run_dir)
        except ValueError:
            pass
        alias = 'input_' + key
        dest = run_dir / alias
        if dest.exists() or dest.is_symlink():
            raise FileExistsError('Duplicate localized input: ' + str(dest))
        if src.is_dir():
            shutil.copytree(src, dest, symlinks=False)
        else:
            shutil.copy2(src, dest)
        return alias

    def get_control(self, run_dir):
        contents = self.template.read_text(encoding='utf-8', errors='replace')
        if '&control' not in contents.lower():
            raise ValueError('Not a TDI control namelist: ' + str(self.template))
        for key in ('firstsrc', 'lastsrc', 'firstdat', 'lastdat'):
            m = re.search(r'(?im)^\s*' + key + r'\s*=\s*(\d+)', contents)
            target = self.yr_s if key.startswith('first') else self.yr_e
            if not m or int(m.group(1)) != target:
                raise ValueError('TDI control incorrect period ' + key + ': expected ' + str(target))
        self._check_switch(contents, 'tc3tdi', ('t', 'true'))
        self._check_switch(contents, 'pseudodata', ('f', 'false'))
        for key in self.STATIC_FILES:
            filename = self.get_control_value(contents, key)
            localized = self.localize_control_file(run_dir, key, filename)
            contents = self.replace_control_value(contents, key, localized)
        for key, value in (('sourcef', 'sourcepriortdi.nc'),
                           ('datf', 'obsrvCH4_9925.nc'),
                           ('datflist', 'obsrvCH4_sites.txt')):
            contents = self.replace_control_value(contents, key, value)
        for key in self.CONTROL_FILES:
            value = self.get_control_value(contents, key)
            if not (run_dir / value).exists():
                raise FileNotFoundError('Control ' + key + ' is not accessible: ' + str(run_dir / value))
        return contents

    @staticmethod
    def sha256(path):
        h = hashlib.sha256()
        with Path(path).open('rb') as f:
            for block in iter(lambda: f.read(1024 * 1024), b''):
                h.update(block)
        return h.hexdigest()

    @staticmethod
    def write_json(path, payload):
        path = Path(path)
        temporary = path.with_name(path.name + '.tmp')
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        temporary.replace(path)

    # --------------------------------------------------------
    # --- Create unique TDI workdir. Never overwrite original results.
    def prepare_inversion(self, prior, unc, model, source, obs, run_dir):
        if run_dir.exists():
            raise FileExistsError('D4 run directory exists: ' + str(run_dir))
        run_dir.parent.mkdir(parents=True, exist_ok=True)
        run_dir.mkdir()
        try:
            shutil.copytree(self.files_dir, run_dir / 'files4run', symlinks=False)
            # --- Green functions are SHARED outside run_dir; never copy or symlink them.
            shutil.copy2(source, run_dir / 'sourcepriortdi.nc')
            shutil.copy2(obs, run_dir / 'obsrvCH4_9925.nc')
            shutil.copy2(Path(self.sites_f), run_dir / 'obsrvCH4_sites.txt')
            shutil.copy2(self.executable, run_dir / 'tdi_ch4')
            contents = self.get_control(run_dir)
            (run_dir / 'control').write_text(contents, encoding='utf-8')
            manifest = {'case': self.case, 'icase': self.icase, 'prior': prior,
                        'unc': unc, 'model': model,
                        'source_D1_prior': str(source), 'source_D3_observations': str(obs),
                        'source_sites': str(Path(self.sites_f).resolve()),
                        'source_executable': str(self.executable),
                        'source_control_template': str(self.template),
                        'source_files4run': str(self.files_dir),
                        'GF_path': str(self.gf_dir),
                        'D1_prior_sha256': self.sha256(source),
                        'D3_observations_sha256': self.sha256(obs),
                        'control_sha256': self.sha256(run_dir / 'control'),
                        'created_utc': datetime.now(timezone.utc).isoformat(),
                        'GF_note': 'Shared Green functions in inv_dir/gmats; no copy, no symlink',
                        'input_mode': 'LOCAL_COPIES_SHARED_GF_NO_SYMLINKS',
                        'local_D1_prior_sha256': self.sha256(run_dir / 'sourcepriortdi.nc'),
                        'local_D3_observations_sha256': self.sha256(run_dir / 'obsrvCH4_9925.nc'),
                        'local_TDI_sha256': self.sha256(run_dir / 'tdi_ch4')}
            (run_dir / 'RUN_INPUTS.txt').write_text(
                ''.join(str(key) + '=' + str(value) + '\n' for key, value in manifest.items()),
                encoding='utf-8')
            links = [str(path) for path in run_dir.rglob('*') if path.is_symlink()]
            if links:
                raise RuntimeError('D4 symlinks forbidden in run directory: ' + '; '.join(links))
            self.write_json(run_dir / 'RUN_STATUS.json', {'status': 'prepared',
                            'experiment': [self.case, prior, unc, model]})
        except Exception:
            # Preserve incomplete work directory for debugging; no silent cleanup.
            raise
        return run_dir

    # --------------------------------------------------------
    # --- Verify actual results; Fortran STOP can return exit status zero.
    def check_run(self, run_dir, strict=True):
        problems = []
        for name in self.EXPECTED_RESULTS:
            item = run_dir / name
            if not item.is_file() or item.stat().st_size == 0:
                problems.append('missing_or_empty_' + name)
        logs = sorted(run_dir.glob('tdi_*.out'))
        if not logs:
            problems.append('missing_tdi_log')
        else:
            log = logs[-1].read_text(encoding='utf-8', errors='replace')
            # Match fatal messages, not harmless words like 'error variance'.
            if re.search(r'(?im)^\s*(forrtl:|segmentation fault|fatal error:|cannot open file|error opening GF file)', log):
                problems.append('fatal_log_message')
        # --- TDI can exit with code zero and still write a full all-zero solution.
        # --- A decoded all-zero result is not a successful inversion.
        predicted_file = run_dir / 'predicted.sources'
        if predicted_file.is_file() and predicted_file.stat().st_size:
            try:
                coefficients = self.read_source_coefficients(predicted_file,
                                       months=self.ntime, first_year=self.yr_s)
            except Exception as exc:
                problems.append('predicted_sources_decoder_error: ' + str(exc))
            else:
                if coefficients is None:
                    problems.append('unrecognized_predicted_sources_format')
                elif not np.isfinite(coefficients).all():
                    problems.append('nonfinite_predicted_sources')
                elif not np.any(np.abs(coefficients) > 1e-12):
                    problems.append('all_zero_TDI_corrections_CHECK_SOLVER_AND_GREENS')
        status = 'PASS' if not problems else 'FAIL'
        if strict and problems:
            raise RuntimeError('TDI output check failed: ' + ', '.join(problems) + ' in ' + str(run_dir))
        return status, problems

    # --------------------------------------------------------
    # --- Verify modern AND original TDI RUN_INPUTS provenance.
    #     Older manifests store SHA256 values as '<hash>  local-file-name',
    #     not D1_prior_sha256=... / D3_observations_sha256=... key-value pairs.
    def check_provenance(self, prior, unc, model, source, obs, run_dir):
        manifest = run_dir / 'RUN_INPUTS.txt'
        if not manifest.is_file():
            return False, ['RUN_INPUTS.txt is missing']

        lines = manifest.read_text(encoding='utf-8', errors='replace').splitlines()
        metadata = {}
        legacy_hashes = {}
        for line in lines:
            match = re.fullmatch(r'\s*([0-9a-fA-F]{64})\s+\*?(.+?)\s*', line)
            if match:
                file_name = Path(match.group(2).strip().removeprefix('./')).as_posix()
                legacy_hashes[file_name] = match.group(1).lower()
            elif '=' in line and not line.lstrip().startswith('#'):
                key, value = line.split('=', 1)
                metadata[key.strip()] = value.strip()

        issues = []
        for key, expected in (('case', self.case), ('prior', prior),
                              ('unc', unc), ('model', model)):
            if metadata.get(key) != expected:
                issues.append(key + ' mismatch: recorded ' + repr(metadata.get(key)) +
                              ', expected ' + repr(expected))
        if 'icase' in metadata and metadata['icase'] != self.icase:
            issues.append('icase mismatch: recorded ' + metadata['icase'] +
                          ', expected ' + self.icase)

        for label, current, recorded_paths, digest_key, local_name in (
                ('D1 prior', Path(source), ('source_D1_prior',), 'D1_prior_sha256', 'sourcepriortdi.nc'),
                ('D3 observations', Path(obs), ('source_D3_observations',), 'D3_observations_sha256', 'obsrvCH4_9925.nc')):
            paths = [metadata[key] for key in recorded_paths if metadata.get(key)]
            if paths:
                for path in paths:
                    if Path(path).expanduser().resolve() != current.resolve():
                        issues.append(label + ' source path mismatch: recorded ' + path +
                                      ', expected ' + str(current))
            else:
                # --- Older manifests used differently named provenance keys.
                # --- Match the actual source basename, not historical stage labels.
                previous = [value for key, value in metadata.items()
                            if key.startswith('source_') and Path(value).name == current.name]
                if len(previous) == 1:
                    if Path(previous[0]).expanduser().resolve() != current.resolve():
                        issues.append(label + ' source path mismatch: recorded ' + previous[0] +
                                      ', expected ' + str(current))
                else:
                    issues.append(label + ' input path is absent from RUN_INPUTS.txt')

            # --- Prefer a checksum recorded at the time the inversion ran.
            recorded_hash = metadata.get(digest_key) or legacy_hashes.get(local_name)
            if recorded_hash is not None:
                actual_hash = self.sha256(current)
                if recorded_hash.lower() != actual_hash:
                    issues.append(label + ' SHA256 mismatch: recorded ' +
                                  recorded_hash[:12] + ', current ' + actual_hash[:12])
                # --- Independent local input must also match the run manifest.
                local_copy = run_dir / local_name
                if local_copy.is_file() and not local_copy.is_symlink():
                    if self.sha256(local_copy) != recorded_hash.lower():
                        issues.append(label + ' local copy SHA256 mismatch')
            else:
                # --- Original run may contain an independent COPY of the input.
                #     A symlink without a recorded hash cannot prove unchanged
                #     inputs, since its target might have been modified later.
                local = run_dir / local_name
                if local.is_file() and not local.is_symlink():
                    if self.sha256(local) != self.sha256(current):
                        issues.append(label + ' copied run input differs from current source')
                else:
                    issues.append(label + ' has no saved SHA256 or independent input copy; '
                                  'cannot safely verify the historical run')

        return not issues, issues

    def run_inversion(self, prior, unc, model, source, obs, run_dir):
        if run_dir.exists():
            # --- An older run may have been made with symlinks. Never silently
            # --- claim copy-only provenance or mutate a finished inversion.
            existing_links = [str(p.relative_to(run_dir)) for p in run_dir.rglob('*') if p.is_symlink()]
            if existing_links:
                raise FileExistsError('Existing D4 run contains symlinks (copy-only mode): ' +
                                      str(run_dir) + '\nLinks: ' + ', '.join(existing_links) +
                                      '\nMove or delete this old run before a fresh copy-only calculation.')
            status, output_issues = self.check_run(run_dir, strict=False)
            provenance_ok, provenance_issues = self.check_provenance(
                prior, unc, model, source, obs, run_dir)
            if status == 'PASS' and provenance_ok:
                print('D4: reuse verified run:', run_dir)
                return 'REUSED'
            issues = output_issues + provenance_issues
            raise FileExistsError('D4 run already exists: ' + str(run_dir) +
                                  '\nWill not overwrite existing inversion results.\nReasons:\n  - ' +
                                  '\n  - '.join(issues))
        self.prepare_inversion(prior, unc, model, source, obs, run_dir)
        log_file = run_dir / ('tdi_' + prior + '_' + model + '_' + unc + '.out')
        self.write_json(run_dir / 'RUN_STATUS.json', {'status': 'running',
                        'experiment': [self.case, prior, unc, model]})
        try:
            print('D4 run:', run_dir)
            with log_file.open('w', encoding='utf-8') as log:
                subprocess.run([str(run_dir / 'tdi_ch4')], cwd=str(run_dir),
                               stdout=log, stderr=subprocess.STDOUT, check=True)
            self.check_run(run_dir, strict=True)
            self.write_json(run_dir / 'RUN_STATUS.json', {'status': 'finished',
                            'experiment': [self.case, prior, unc, model]})
            print('D4 inversion PASS:', prior, unc, model)
            return 'PASS'
        except Exception as exc:
            self.write_json(run_dir / 'RUN_STATUS.json', {'status': 'failed',
                            'experiment': [self.case, prior, unc, model],
                            'error': str(exc)})
            raise

    def run(self, dry_run=False):
        self.require_file(self.template)
        for prior, unc, model in self.cases():
            source, obs, run_dir, _ = self.check_inputs(prior, unc, model, dry_run=dry_run)
            if dry_run:
                # Examine control without creating any file: check namelist fields/period.
                text = self.template.read_text(encoding='utf-8', errors='replace')
                self._check_switch(text, 'tc3tdi', ('t', 'true'))
                self._check_switch(text, 'pseudodata', ('f', 'false'))
                print('D4 DRY-RUN:', run_dir)
            else:
                self.run_inversion(prior, unc, model, source, obs, run_dir)
        print('\nD4 INVERSION: PASS' if not dry_run else '\nD4 INVERSION DRY-RUN: PASS')

    # ============================================================
    # === Part 2: read-only diagnostics
    def find_runs(self):
        root = self.result_root
        for prior_dir in sorted(root.glob(self.icase + '_p*')):
            if not prior_dir.is_dir():
                continue
            for run_dir in sorted(prior_dir.glob(self.icase + '_*_inv*')):
                if not run_dir.is_dir():
                    continue
                yield run_dir

    def diagnose_one(self, run_dir):
        name = run_dir.name
        prior = run_dir.parent.name.removeprefix(self.icase + '_')
        m = re.fullmatch(re.escape(self.icase) + r'_(ctl|ux2|ux4)_(inv\d+)', name)
        if not m:
            return None
        unc, model = m.groups()
        manifest = run_dir / 'RUN_INPUTS.txt'
        source = run_dir / 'sourcepriortdi.nc'
        obs = run_dir / 'obsrvCH4_9925.nc'
        row = {'case': self.case, 'icase': self.icase, 'prior': prior,
               'unc': unc, 'model': model, 'path': str(run_dir)}
        row['status'], issues = self.check_run(run_dir, strict=False)
        row['problems'] = ';'.join(issues)
        row['manifest'] = int(manifest.is_file())
        for label, item in (('prior_sources', run_dir / 'prior.sources'),
                            ('predicted_sources', run_dir / 'predicted.sources')):
            row[label + '_bytes'] = item.stat().st_size if item.is_file() else 0
        row['input_prior'] = int(source.is_file())
        row['input_observations'] = int(obs.is_file())
        row['control'] = int((run_dir / 'control').is_file())
        row['files4run'] = int((run_dir / 'files4run').is_dir())
        row['gf_shared'] = int(self.gf_dir.is_dir() and 
                                self.get_control_value((run_dir / 'control').read_text(), 'greenf') ==
                                str(self.gf_dir.resolve()) + os.sep) if (run_dir / 'control').is_file() else 0
        if manifest.is_file():
            meta = dict(line.split('=', 1) for line in manifest.read_text(encoding='utf-8').splitlines()
                        if '=' in line)
            for key in ('case', 'prior', 'unc', 'model'):
                if key in meta and meta[key] != row[key]:
                    row['status'] = 'FAIL'
                    row['problems'] += ';manifest_mismatch_' + key
        row['n_valid_obs'] = 0
        row['n_stations'] = 0
        row['mean_sigma_ppb'] = np.nan
        row['input_rmse_ppb'] = np.nan
        if source.is_file() and obs.is_file():
            try:
                obs_info = self.check_netcdf(source, obs)
                row['n_valid_obs'] = obs_info['observations_valid_station_months']
                row['n_stations'] = obs_info['observations_stations']
                data, sigma = self.get_observation_data(obs)
                row['mean_sigma_ppb'] = float(np.nanmean(sigma))
                row['input_rmse_ppb'] = float(np.sqrt(np.nanmean(data * data)))
            except Exception as exc:
                row['status'] = 'FAIL'
                row['problems'] += ';netcdf:' + str(exc)
        return row

    def get_observation_data(self, path):
        import xarray as xr
        with xr.open_dataset(path, decode_times=False, mask_and_scale=False) as ds:
            data = np.asarray(ds['tracer_data'].transpose('tlen', 'station_count', 'data_type').values,
                              dtype=np.float64)
        obs, sigma = data[..., 0], data[..., 1]
        good = np.isfinite(obs) & np.isfinite(sigma)
        good &= (obs > -900) & (sigma > 0) & (sigma < 1e20)
        return np.where(good, obs, np.nan), np.where(good, sigma, np.nan)

    def plot_observation_diagnostics(self, obs_path):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import pandas as pd

        data, sigma = self.get_observation_data(obs_path)
        if data.shape[0] != self.ntime:
            raise ValueError('Unexpected D4 observation months')
        x = pd.date_range(str(self.yr_s) + '-01-01', periods=self.ntime, freq='MS')
        n = np.isfinite(data).sum(axis=1)
        # --- Missing station-months are allowed; do not average empty months
        # --- Count is calculated from the same D3 observations validity mask as input checks
        means = np.full(self.ntime, np.nan, dtype=np.float64)
        sigmas = np.full(self.ntime, np.nan, dtype=np.float64)
        np.divide(np.nansum(data, axis=1), n, out=means, where=n > 0)
        nsigma = np.isfinite(sigma).sum(axis=1)
        np.divide(np.nansum(sigma, axis=1), nsigma, out=sigmas, where=nsigma > 0)
        out = self.diag_dir
        out.mkdir(parents=True, exist_ok=True)
        with (out / 'observation_monthly.csv').open('w', newline='') as f:
            w = csv.writer(f)
            w.writerow(('month', 'n_valid', 'mean_obs_minus_actm_ppb', 'mean_sigma_ppb'))
            for i in range(self.ntime):
                w.writerow((x[i].strftime('%Y-%m'), int(n[i]), float(means[i]), float(sigmas[i])))
        with (out / 'observation_stations.csv').open('w', newline='') as f:
            w = csv.writer(f)
            w.writerow(('station_index', 'valid_months', 'bias_obs_minus_actm_ppb', 'rmse_ppb', 'mean_sigma_ppb'))
            for j in range(data.shape[1]):
                yy = data[:, j]
                w.writerow((j + 1, int(np.isfinite(yy).sum()),
                            float(np.nanmean(yy)) if np.isfinite(yy).any() else np.nan,
                            float(np.sqrt(np.nanmean(yy * yy))) if np.isfinite(yy).any() else np.nan,
                            float(np.nanmean(sigma[:, j])) if np.isfinite(sigma[:, j]).any() else np.nan))

        fig, axes = plt.subplots(3, 1, figsize=(13, 10), sharex=True)
        axes[0].plot(x, n)
        axes[0].set_ylabel('Valid station-months')
        axes[1].plot(x, means)
        axes[1].set_ylabel('OBS - ACTM (ppb)')
        axes[2].plot(x, sigmas)
        axes[2].set_ylabel('Mean sigma (ppb)')
        for ax in axes:
            ax.axvspan(x[314], x[-1], color='gray', alpha=0.15)
        axes[0].text(0.99, 0.90, 'Mar-Dec 2025: OBS not yet available',
                     ha='right', va='top', transform=axes[0].transAxes, fontsize=8)
        fig.suptitle('D4 D3 observations inputs (not posterior fit): ' + self.case)
        fig.tight_layout()
        fig.savefig(out / 'observations_timeseries.png', dpi=180)
        plt.close(fig)

        # --- Show only active D3 observations rows, retaining source station IDs on the axis
        selected = np.flatnonzero(np.isfinite(data).any(axis=0))
        order = selected[np.argsort(np.isfinite(data[:, selected]).sum(axis=0))[::-1]]
        fig, ax = plt.subplots(figsize=(13, 6))
        ax.imshow(np.isfinite(data[:, order].T), origin='lower', aspect='auto',
                  interpolation='nearest', cmap='viridis', extent=[1999, 2026, 0.5, len(order) + 0.5])
        ticks = np.unique(np.linspace(0, len(order) - 1, min(12, len(order)), dtype=int))
        ax.set_yticks(ticks + 1)
        ax.set_yticklabels([str(order[k] + 1) for k in ticks])
        ax.set_xlabel('Year')
        ax.set_ylabel('D3 observations station index (active only)')
        ax.set_title('D3 observations coverage: ' + self.case + ' (' + str(len(order)) + ' active of ' +
                     str(data.shape[1]) + ' positions)')
        fig.tight_layout()
        fig.savefig(out / 'observations_coverage.png', dpi=180)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(9, 5))
        values = sigma[np.isfinite(sigma)]
        ax.hist(values, bins=50)
        ax.set_xlabel('Observation uncertainty (ppb)')
        ax.set_ylabel('Station-month count')
        ax.set_title('D3 observations uncertainty distribution')
        fig.tight_layout()
        fig.savefig(out / 'observation_uncertainty.png', dpi=180)
        plt.close(fig)

        return {'n_months': int(data.shape[0]), 'n_stations': int(data.shape[1]),
                'n_valid': int(np.isfinite(data).sum()),
                'n_empty_months': int(np.count_nonzero(n == 0)),
                'n_nonempty_stations': int(np.count_nonzero(np.isfinite(data).any(axis=0)))}

    # --------------------------------------------------------
    # --- D1 prior uncertainty plots. D1 prior fields are NOT posterior emissions.
    def plot_prior_inputs(self, rows):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import pandas as pd
        import xarray as xr
        obs_dates = pd.date_range(str(self.yr_s) + '-01-01', periods=self.ntime, freq='MS')
        self.diag_dir.mkdir(parents=True, exist_ok=True)
        curves = []
        shown = set()
        for row in rows:
            prior = row['prior']
            if prior in shown:
                continue
            shown.add(prior)
            source = Path(row['path']) / 'sourcepriortdi.nc'
            if not source.is_file():
                continue
            with xr.open_dataset(source, decode_times=False, mask_and_scale=False) as ds:
                if 'sources_land' not in ds:
                    continue
                dims = ('month', 'land_region', 'datatype_count')
                if not all(dim in ds['sources_land'].dims for dim in dims):
                    continue
                v = np.asarray(ds['sources_land'].transpose(*dims).values, dtype=np.float64)
                times = np.asarray(ds['time'].values, dtype=np.float64)
            sel = (times >= self.yr_s) & (times < self.yr_e + 1)
            if sel.sum() != self.ntime or v.shape[1] != 54 or v.shape[2] != 2:
                continue
            v = v[sel]
            self.diag_dir.mkdir(parents=True, exist_ok=True)
            curves.append((prior, np.sum(v[:, :, 1], axis=1)))
            fig, ax = plt.subplots(figsize=(12, 5))
            im = ax.imshow(v[:, :, 1].T, aspect='auto', origin='lower', interpolation='nearest')
            fig.colorbar(im, ax=ax, label='Prior uncertainty (Tg CH4/yr equivalent)')
            ax.set_xlabel('Month index, 1999-2025')
            ax.set_ylabel('Land region 1..54')
            ax.set_title('D1 prior uncertainty: ' + prior)
            fig.tight_layout()
            fig.savefig(self.diag_dir / ('prior_uncertainty_regions_' + prior + '.png'), dpi=180)
            plt.close(fig)
            with (self.diag_dir / ('prior_region_summary_' + prior + '.csv')).open('w', newline='') as f:
                w = csv.writer(f)
                w.writerow(('region', 'mean_prior_correction_Tg_per_year', 'mean_uncertainty_Tg_per_year'))
                for i in range(54):
                    w.writerow((i + 1, float(np.mean(v[:, i, 0])), float(np.mean(v[:, i, 1]))))
        # --- With one prior, show a monthly climatology rather than a misleading
        # --- one-series "sensitivity" graph (sum of sigmas is NOT global sigma).
        legacy = self.diag_dir / 'prior_uncertainty_sensitivity.png'
        if len(curves) == 1:
            legacy.unlink(missing_ok=True)
            prior, values = curves[0]
            climate = [float(np.mean(values[12:312][m::12])) for m in range(12)]
            fig, ax = plt.subplots(figsize=(9, 4.5))
            ax.plot(np.arange(1, 13), climate, '-o', ms=4)
            ax.set_xticks(np.arange(1, 13))
            ax.set_xlabel('Month of year')
            ax.set_ylabel('Sum of input regional sigma (NOT global 1σ)')
            ax.set_title('D1 prior uncertainty annual cycle: ' + prior + ', 2000-2024')
            fig.tight_layout()
            fig.savefig(self.diag_dir / ('prior_uncertainty_climatology_' + prior + '.png'), dpi=160)
            plt.close(fig)
        elif len(curves) > 1:
            fig, ax = plt.subplots(figsize=(10, 5))
            for label, values in curves:
                climate = [float(np.mean(values[12:312][m::12])) for m in range(12)]
                ax.plot(np.arange(1, 13), climate, '-o', ms=3, label=label)
            ax.set_xticks(np.arange(1, 13))
            ax.legend()
            ax.set_xlabel('Month of year')
            ax.set_ylabel('Sum of input regional sigma (NOT global 1σ)')
            ax.set_title('D1 prior uncertainty: monthly climatology, 2000-2024')
            fig.tight_layout()
            fig.savefig(legacy, dpi=160)
            plt.close(fig)
        return len(curves)

    def plot_experiment_comparison(self, rows):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        valid = [r for r in rows if r['n_valid_obs'] > 0]
        if len(valid) < 2:
            return
        labels = [r['prior'] + '_' + r['unc'] + '_' + r['model'] for r in valid]
        x = np.arange(len(labels))
        fig, axes = plt.subplots(2, 1, figsize=(max(10, len(labels) * 1.1), 8))
        axes[0].bar(x, [r['n_valid_obs'] for r in valid])
        axes[0].set_ylabel('Valid D3 observations station-months')
        axes[1].bar(x, [r['mean_sigma_ppb'] for r in valid])
        axes[1].set_ylabel('Mean D3 observations sigma (ppb)')
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(labels, rotation=45, ha='right')
        fig.suptitle('D4 experiment input comparison: ' + self.case)
        fig.tight_layout()
        fig.savefig(self.diag_dir / 'experiment_inputs.png', dpi=180)
        plt.close(fig)

    def write_log_diagnostics(self, rows):
        # --- Extract exact diagnostic log excerpts; do not invent chi-squared or convergence.
        key = re.compile(r'(?i)(chi.?square|converg|condition number|singular|svd|valid months|ntotdat|neqns|warning|error)')
        with (self.diag_dir / 'tdi_log_diagnostics.txt').open('w') as out:
            for row in rows:
                out.write('=== ' + row['prior'] + ' ' + row['unc'] + ' ' + row['model'] + ' ===\n')
                logs = sorted(Path(row['path']).glob('tdi_*.out'))
                if not logs:
                    out.write('No TDI log\n')
                    continue
                count = 0
                with logs[-1].open(encoding='utf-8', errors='replace') as f:
                    for line in f:
                        if key.search(line):
                            out.write(line[:300].rstrip() + '\n')
                            count += 1
                            if count >= 150:
                                out.write('(truncated)\n')
                                break
                if not count:
                    out.write('No recognizable convergence/SVD statistics in log\n')

    # --------------------------------------------------------
    # --- Strictly labelled, unit-agnostic source coefficient parsing
    @staticmethod
    def read_source_coefficients(path, months=324, regions=54, first_year=1999):
        """Read the regional TDI source format expected by D5.

        For each of 54 land regions, D5's read_pred_srcs() reads:
            land_01                  (A20)
            <header line>            (A)
            1999.041667 <coefficient> (list-directed)
            ... 324 monthly pairs ...

        Numeric-only matrices are retained for backward-compatible tests.
        The values are TDI regional coefficients; they are not posterior
        gridded CH4 fluxes (which are calculated later in D5).
        """
        path = Path(path)
        if not path.is_file() or path.stat().st_size > 30 * 1024 * 1024:
            return None
        try:
            lines = path.read_text(encoding='ascii').splitlines()
        except (UnicodeDecodeError, OSError):
            return None
        # --- Try the actual block-oriented TDI file first
        rows = [line.strip() for line in lines if line.strip()]
        if rows and re.fullmatch(r'land_0*1', rows[0], re.IGNORECASE):
            result = np.empty((months, regions), dtype=np.float64)
            reference_time = np.empty(months, dtype=np.float64)
            index = 0
            for region in range(regions):
                # --- A blank separator between regions is harmless
                while index < len(lines) and not lines[index].strip():
                    index += 1
                if index + 2 + months > len(lines):
                    return None
                name = lines[index].strip()
                match = re.fullmatch(r'land_(\d+)', name, re.IGNORECASE)
                if match is None or int(match.group(1)) != region + 1:
                    return None
                # --- Fortran consumes exactly one header record, even if empty
                index += 2
                for k in range(months):
                    fields = lines[index].replace('D', 'E').replace('d', 'e').split()
                    if len(fields) < 2:
                        return None
                    try:
                        time_value = float(fields[0])
                        coefficient = float(fields[1])
                    except ValueError:
                        return None
                    if not np.isfinite(time_value) or not np.isfinite(coefficient):
                        return None
                    expected = first_year + k // 12 + ((k % 12) + 0.5) / 12.0
                    if abs(time_value - expected) > 0.05:
                        return None
                    if region > 0 and abs(time_value - reference_time[k]) > 1e-5:
                        return None
                    result[k, region] = coefficient
                    if region == 0:
                        reference_time[k] = time_value
                    index += 1
            return result

        # --- Legacy reader for plain numeric matrices (synthetic testing)
        numeric = []
        for line in rows:
            if line.startswith(('#', '!', 'c ')):
                continue
            try:
                numeric.append([float(value.replace('D', 'E').replace('d', 'e'))
                                for value in line.split()])
            except ValueError:
                return None
        if not numeric or len(set(map(len, numeric))) != 1:
            return None
        arr = np.asarray(numeric, dtype=np.float64)
        if not np.isfinite(arr).all():
            return None
        if arr.shape == (regions, months):
            return arr.T
        if arr.shape == (months, regions):
            return arr
        return None

    # --------------------------------------------------------
    # --- Extract verified log indicators, not speculative convergence claims
    @staticmethod
    def read_tdi_log_metrics(run_dir):
        logs = sorted(Path(run_dir).glob('tdi_*.out'))
        result = {'nobs_tdi': None, 'ntotdat_tdi': None, 'excluded_tdi': None,
                  'neqns_tdi': None, 'nunks_tdi': None,
                  'dosvd': None, 'numsvdvecs': None}
        if not logs:
            return result
        text = logs[-1].read_text(encoding='utf-8', errors='replace')
        patterns = {
            'dosvd': r'(?im)^\s*DOSVD\s*=\s*([TF])',
            'numsvdvecs': r'(?im)^\s*NUMSVDVECS\s*=\s*(-?\d+)',
            'nobs_tdi': r'(?i)nobs,\s*valid\s+ntotdat,\s*excluded\s*=\s*(\d+)',
            'ntotdat_tdi': r'(?i)nobs,\s*valid\s+ntotdat,\s*excluded\s*=\s*\d+\s+(\d+)',
            'excluded_tdi': r'(?i)nobs,\s*valid\s+ntotdat,\s*excluded\s*=\s*\d+\s+\d+\s+(\d+)',
            'neqns_tdi': r'(?i)ndat,ntotdat,nconsts,neqns\s+\d+\s+\d+\s+\d+\s+(\d+)',
            'nunks_tdi': r'(?i)nunks\s*=\s*(\d+)'
        }
        for key, pattern in patterns.items():
            match = re.search(pattern, text)
            if match:
                result[key] = match.group(1) if key == 'dosvd' else int(match.group(1))
        return result

    # --------------------------------------------------------
    # --- Separate filesystem PASS from the numerical quality of TDI solution
    def evaluate_source_quality(self, row, predicted, prior=None):
        label = row['prior'] + '_' + row['unc'] + '_' + row['model']
        metrics = self.read_tdi_log_metrics(row['path'])
        finite = np.isfinite(predicted).all()
        zero = bool(finite and np.max(np.abs(predicted)) < 1e-8)
        n_nonzero = int(np.count_nonzero(np.abs(predicted) >= 1e-8)) if finite else 0
        flags = []
        if not finite:
            flags.append('NONFINITE_COEFFICIENTS')
        if zero:
            flags.append('ALL_ZERO_PREDICTED_SOURCES')
        # --- NUMSVDVECS may control SVD output only; do not call it a solver error.
        # --- Report the value but flag only independently demonstrated failures.
        if metrics['nobs_tdi'] is not None and metrics['nobs_tdi'] == 0:
            flags.append('TDI_SELECTED_NO_OBSERVATIONS')
        return {'label': label, 'case': self.case, 'status_files': row['status'],
                'status_numerical': 'CRITICAL' if zero or not finite else 'REVIEW' if flags else 'OK',
                'n_nonzero_coefficients': n_nonzero,
                'max_abs_coefficient': float(np.max(np.abs(predicted))) if finite else np.nan,
                'mean_abs_coefficient_2000_2024': float(np.mean(np.abs(predicted[12:312]))) if finite else np.nan,
                'n_regions': predicted.shape[1], 'n_months': predicted.shape[0],
                'flags': ';'.join(flags), 'nobs_tdi': metrics['nobs_tdi'],
                'ntotdat_tdi': metrics['ntotdat_tdi'], 'excluded_tdi': metrics['excluded_tdi'],
                'dosvd': metrics['dosvd'], 'numsvdvecs': metrics['numsvdvecs'],
                'nunks_tdi': metrics['nunks_tdi'],
                'prior_max_abs': float(np.max(np.abs(prior))) if prior is not None else np.nan}

    # --------------------------------------------------------
    # --- Quantitative annual and regional diagnostics for 2000-2024
    def write_source_statistics(self, label, predicted, dates):
        core = predicted[12:312, :]
        years = np.arange(self.yr_s, self.yr_e + 1)
        curve = np.sum(predicted, axis=1)
        with (self.diag_dir / ('source_correction_' + label + '.csv')).open('w', newline='') as f:
            w = csv.writer(f)
            w.writerow(('month', 'sum_54_tdi_regional_corrections_native'))
            for dt, value in zip(dates, curve):
                w.writerow((dt.strftime('%Y-%m'), float(value)))
        with (self.diag_dir / ('source_annual_' + label + '.csv')).open('w', newline='') as f:
            w = csv.writer(f)
            w.writerow(('year', 'mean_monthly_sum_native', 'min_monthly_sum_native',
                        'max_monthly_sum_native', 'interpretation'))
            for year in years:
                part = curve[(year - self.yr_s) * 12:(year - self.yr_s + 1) * 12]
                note = 'initial boundary' if year == self.yr_s else ('observations missing after Feb' if year == 2025 else 'analysis period')
                w.writerow((year, float(np.mean(part)), float(np.min(part)), float(np.max(part)), note))
        regional = []
        for i in range(predicted.shape[1]):
            xx = core[:, i]
            season = np.array([np.mean(xx[month::12]) for month in range(12)])
            annual = xx.reshape(25, 12).mean(axis=1)
            slope = float(np.polyfit(np.arange(2000, 2025), annual, 1)[0])
            regional.append({'region': i + 1, 'mean_2000_2024': float(np.mean(xx)),
                             'mean_abs_2000_2024': float(np.mean(np.abs(xx))),
                             'std_2000_2024': float(np.std(xx)),
                             'seasonal_amplitude_2000_2024': float(np.max(season) - np.min(season)),
                             'slope_native_per_year_2000_2024': slope,
                             'min_all': float(np.min(predicted[:, i])),
                             'max_all': float(np.max(predicted[:, i])),
                             'negative_months_2000_2024': int((xx < 0).sum()),
                             'nonzero_months_2000_2024': int((np.abs(xx) >= 1e-8).sum())})
        with (self.diag_dir / ('source_regional_' + label + '.csv')).open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(regional[0]))
            writer.writeheader()
            writer.writerows(regional)
        return regional, curve

    # --------------------------------------------------------
    # --- CYC vs INCA: use matching runs only; no TDI rerun
    def plot_case_comparison(self, row, predicted, dates):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        peer = 'INCA' if self.case == 'CYC' else 'CYC'
        other_dir = (self.result_root.parent / peer / (self.icase + '_' + row['prior']) /
                     (self.icase + '_' + row['unc'] + '_' + row['model']))
        other = self.read_source_coefficients(other_dir / 'predicted.sources',
                                              self.ntime, first_year=self.yr_s)
        if other is None or other.shape != predicted.shape:
            return False
        label = row['prior'] + '_' + row['unc'] + '_' + row['model']
        curves = {self.case: predicted.sum(axis=1), peer: other.sum(axis=1)}
        diff = (predicted - other)[12:312, :]
        fig, ax = plt.subplots(3, 1, figsize=(13, 10), constrained_layout=True)
        for name, values in curves.items():
            ax[0].plot(dates, values, label=name, lw=1.1)
        ax[0].legend()
        ax[0].axhline(0, lw=0.6, color='gray')
        ax[0].set_ylabel('Sum of TDI corrections')
        ax[0].set_title('Matching TDI experiments: ' + label)
        ax[1].plot(dates, curves.get('CYC') - curves.get('INCA'), lw=1)
        ax[1].axhline(0, lw=0.6, color='gray')
        ax[1].set_ylabel('CYC - INCA (native units)')
        ranking = np.argsort(np.mean(np.abs(diff), axis=0))[::-1][:15]
        ax[2].bar(np.arange(15), np.mean(np.abs(diff), axis=0)[ranking])
        ax[2].set_xticks(np.arange(15))
        ax[2].set_xticklabels([str(i + 1) for i in ranking])
        ax[2].set_ylabel('Mean |CYC - INCA|, 2000-2024')
        ax[2].set_xlabel('Land region')
        fig.savefig(self.diag_dir / ('compare_CYC_INCA_' + label + '.png'), dpi=150)
        plt.close(fig)
        return True

    # --------------------------------------------------------
    # --- Scientific one-page overview of the actual inversion coefficients
    def plot_inversion_overview(self, label, predicted, regional, dates, metrics, obs_path):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        from matplotlib.colors import TwoSlopeNorm

        core = predicted[12:312]
        months = dates[12:312]
        series = predicted.sum(axis=1)
        finite = np.isfinite(predicted).all()
        zero = bool(finite and np.max(np.abs(predicted)) < 1e-8)
        obs, _ = self.get_observation_data(obs_path) if obs_path.is_file() else (None, None)

        fig, axes = plt.subplots(2, 2, figsize=(16, 10), constrained_layout=True)
        fig.suptitle('D4 | ' + self.case + ' | ' + label + ' | TDI inversion quality', fontsize=16)

        # --- Monthly global sum; 1999 and pending 2025 are not in the core statistics
        ax = axes[0, 0]
        ax.plot(dates, series, lw=0.8, alpha=0.5, color='steelblue', label='Monthly sum')
        annual = series[12:312].reshape(25, 12).mean(axis=1)
        ax.plot(months[5::12], annual, '-o', ms=2.3, lw=1.5, color='darkblue', label='Annual mean')
        ax.axhline(0, lw=0.7, color='gray')
        ax.axvspan(dates[0], dates[11], alpha=0.13, color='gray')
        ax.axvspan(dates[314], dates[-1], alpha=0.13, color='gray')
        ax.set_title('Sum of 54 regional corrections (native TDI units)')
        ax.set_ylabel('Correction sum')
        ax.xaxis.set_major_locator(mdates.YearLocator(5))
        ax.legend(loc='upper right', fontsize=8)
        if zero:
            ax.text(0.5, 0.5, 'CRITICAL: ALL COEFFICIENTS ZERO',
                    transform=ax.transAxes, ha='center', va='center',
                    color='darkred', weight='bold', fontsize=12)

        # --- Heatmap in physical core period; symmetric limits around zero
        ax = axes[0, 1]
        if finite:
            vmax = float(np.nanpercentile(np.abs(core), 98))
            vmax = max(vmax, 1e-12)
            img = ax.imshow(core.T, origin='lower', aspect='auto', interpolation='nearest',
                            cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                            extent=(2000, 2025, 0.5, predicted.shape[1] + 0.5))
            fig.colorbar(img, ax=ax, shrink=0.82, label='Regional correction (native)')
        else:
            ax.text(0.5, 0.5, 'NON-FINITE CORRECTIONS', transform=ax.transAxes,
                    ha='center', color='darkred')
        ax.set_title('Regional corrections, 2000–2024')
        ax.set_xlabel('Year')
        ax.set_ylabel('Land region')

        # --- Rank the strongest regions, not 54 overlapping lines
        ax = axes[1, 0]
        ranked = sorted(regional, key=lambda r: r['mean_abs_2000_2024'], reverse=True)[:12]
        if zero or not finite:
            ax.axis('off')
            ax.text(0.5, 0.5, 'No meaningful regional ranking\n'
                    '(numerical inversion failed)', ha='center', va='center',
                    transform=ax.transAxes, fontsize=12, color='darkred')
        else:
            vals = [r['mean_2000_2024'] for r in ranked]
            ypos = np.arange(len(ranked))
            ax.barh(ypos, vals, color=['steelblue' if v >= 0 else 'indianred' for v in vals])
            ax.set_yticks(ypos, [f"Region {r['region']:02d}" for r in ranked])
            ax.invert_yaxis()
            ax.axvline(0, color='gray', lw=0.7)
            ax.set_xlabel('Mean signed correction (native)')
        ax.set_title('Top regions, 2000–2024')

        # --- Clarify that valid input ≠ equations selected by TDI
        ax = axes[1, 1]
        if obs is not None:
            count = np.isfinite(obs).sum(axis=1)
            ax.plot(dates, count, lw=0.9, color='darkgreen')
            ax.axvspan(dates[314], dates[-1], color='gray', alpha=0.13)
            ax.text(0.99, 0.95, 'Mar–Dec 2025: observations pending',
                    transform=ax.transAxes, ha='right', va='top', fontsize=8)
            ax.set_ylabel('Valid D3 input station-months')
            ax.set_xlabel('Year')
            ax.xaxis.set_major_locator(mdates.YearLocator(5))
            nvalid = int(np.isfinite(obs).sum())
        else:
            nvalid = 0
            ax.text(0.5, 0.5, 'No local observation NetCDF', transform=ax.transAxes,
                    ha='center')
        ax.set_title('Observational coverage (input, not posterior fit)')
        ax.text(0.01, 0.02,
                f"D3 valid: {nvalid:,}  |  TDI used: {metrics.get('ntotdat_tdi') or 'unknown'}",
                transform=ax.transAxes, fontsize=9, va='bottom', ha='left',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.75))

        out = self.diag_dir / ('inversion_overview_' + self.case + '_' + label + '.png')
        fig.savefig(out, dpi=170)
        plt.close(fig)
        return out

    # --------------------------------------------------------
    # --- Source graphics: real information, not 54 zero-line graphs
    def plot_source_sensitivity(self, rows):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.colors import TwoSlopeNorm
        import pandas as pd
        self.diag_dir.mkdir(parents=True, exist_ok=True)
        dates = pd.date_range(str(self.yr_s) + '-01-01', periods=self.ntime, freq='MS')
        quality = []
        curves = []
        compared = 0
        for row in rows:
            # --- Decode critical zero solutions even though check_run marks them FAIL.
            # --- Otherwise the diagnostic would silently suppress its own warning.
            if row['status'] != 'PASS' and 'all_zero_TDI_corrections_CHECK_SOLVER_AND_GREENS' not in row['problems']:
                continue
            folder = Path(row['path'])
            predicted = self.read_source_coefficients(folder / 'predicted.sources',
                                                      self.ntime, first_year=self.yr_s)
            prior = self.read_source_coefficients(folder / 'prior.sources',
                                                  self.ntime, first_year=self.yr_s)
            if predicted is None:
                print('D4: cannot decode TDI sources:', folder / 'predicted.sources')
                continue
            label = row['prior'] + '_' + row['unc'] + '_' + row['model']
            q = self.evaluate_source_quality(row, predicted, prior)
            quality.append(q)
            regional, curve = self.write_source_statistics(label, predicted, dates)
            self.plot_inversion_overview(label, predicted, regional, dates,
                                         self.read_tdi_log_metrics(folder),
                                         folder / 'obsrvCH4_9925.nc')
            curves.append((label, curve, q['status_numerical']))
            if self.plot_case_comparison(row, predicted, dates):
                compared += 1

            # --- An all-zero inversion is a numerical warning, not a physical result
            zero = q['status_numerical'] == 'CRITICAL'
            output = self.diag_dir / ('source_correction_regions_' + label + '.png')
            fig, ax = plt.subplots(figsize=(13, 6))
            if zero:
                ax.axis('off')
                ax.text(0.5, 0.65, 'ALL 54 REGIONAL CORRECTIONS ARE ZERO',
                        ha='center', va='center', transform=ax.transAxes,
                        weight='bold', fontsize=17)
                ax.text(0.5, 0.42, 'This is NOT a validated inversion result.\n'
                        'Inspect TDI SVD configuration and predicted.sources.',
                        ha='center', va='center', transform=ax.transAxes, fontsize=12)
            else:
                vmax = float(np.max(np.abs(predicted)))
                im = ax.imshow(predicted.T, aspect='auto', origin='lower',
                               cmap='RdBu_r', norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax),
                               interpolation='nearest', extent=[1999, 2026, 0.5, 54.5])
                fig.colorbar(im, ax=ax, label='TDI regional correction, native units')
                ax.axvspan(1999, 2000, color='grey', alpha=0.16)
                ax.axvspan(2025 + 2/12, 2026, color='grey', alpha=0.16)
                ax.set_xlabel('Year')
                ax.set_ylabel('Land region')
            ax.set_title(self.case + ' / ' + label + ' regional TDI corrections')
            fig.tight_layout()
            fig.savefig(output, dpi=160)
            plt.close(fig)

            # --- Previously generated 54 flat plots are misleading; clear ONLY this
            # --- diagnostic experiment subfolder. Never touch TDI run directories.
            region_dir = self.diag_dir / 'regional' / label
            region_dir.mkdir(parents=True, exist_ok=True)
            if zero:
                for png in region_dir.glob('region_??.png'):
                    png.unlink()
                (region_dir / 'CRITICAL_ZERO_SOLUTION.txt').write_text(
                    'All 54 predicted.sources corrections are zero for all 324 months.\n'
                    'These are not meaningful regional inversion results.\n'
                    'Inspect TDI logs, DOSVD and NUMSVDVECS.\n', encoding='utf-8')
                continue
            (region_dir / 'CRITICAL_ZERO_SOLUTION.txt').unlink(missing_ok=True)

            # --- 54 regions summarized as ranked signed average and variability
            by_abs = sorted(regional, key=lambda r: r['mean_abs_2000_2024'], reverse=True)
            ranked = by_abs[:20]
            fig, axs = plt.subplots(2, 1, figsize=(12, 9), constrained_layout=True)
            xx = np.arange(len(ranked))
            axs[0].bar(xx, [r['mean_2000_2024'] for r in ranked])
            axs[0].axhline(0, color='gray', lw=0.7)
            axs[0].set_xticks(xx)
            axs[0].set_xticklabels([str(r['region']) for r in ranked])
            axs[0].set_ylabel('Mean correction, 2000-2024')
            axs[0].set_title('Top 20 regions by mean absolute correction')
            axs[1].bar(xx, [r['std_2000_2024'] for r in ranked])
            axs[1].set_xticks(xx)
            axs[1].set_xticklabels([str(r['region']) for r in ranked])
            axs[1].set_xlabel('Land region')
            axs[1].set_ylabel('Monthly variability, 2000-2024')
            fig.savefig(region_dir / 'regional_ranking.png', dpi=150)
            plt.close(fig)

            # --- Default to the 12 most influential regions; optional full set of 54.
            # --- Always retain a 54-region heatmap and full CSV statistics.
            plot_all = getattr(self, 'd4_plot_all_regions', False)
            indices = list(range(54)) if plot_all else [r['region'] - 1 for r in by_abs[:12]]
            active_png = {'region_' + str(i + 1).zfill(2) + '.png' for i in indices}
            for png in region_dir.glob('region_??.png'):
                if png.name not in active_png:
                    png.unlink()
            for region in indices:
                values = predicted[:, region]
                fig, axes = plt.subplots(2, 1, figsize=(11, 6), constrained_layout=True)
                axes[0].plot(dates, values, lw=0.9, label='TDI correction')
                if prior is not None:
                    axes[0].plot(dates, prior[:, region], lw=0.8, alpha=0.6,
                                 label='prior.sources coefficient')
                axes[0].axhline(0, color='gray', lw=0.7)
                axes[0].axvspan(dates[0], dates[11], alpha=0.12, color='gray')
                axes[0].axvspan(dates[314], dates[-1], alpha=0.12, color='gray')
                axes[0].legend(loc='upper right', fontsize=8)
                axes[0].set_ylabel('TDI coefficients (native)')
                core = values[12:312]
                climate = [np.mean(core[mm::12]) for mm in range(12)]
                axes[1].plot(np.arange(1, 13), climate, '-o', ms=3)
                axes[1].axhline(0, color='gray', lw=0.7)
                axes[1].set_xticks(np.arange(1, 13))
                axes[1].set_xlabel('Month of year')
                axes[1].set_ylabel('Climatological correction, 2000-2024')
                axes[0].set_title(self.case + ' / ' + label + ' / region ' + str(region + 1).zfill(2))
                fig.savefig(region_dir / ('region_' + str(region + 1).zfill(2) + '.png'), dpi=120)
                plt.close(fig)

        with (self.diag_dir / 'source_quality.csv').open('w', newline='') as f:
            if quality:
                writer = csv.DictWriter(f, fieldnames=list(quality[0]))
                writer.writeheader()
                writer.writerows(quality)
            else:
                f.write('No decodable TDI source datasets\n')

        # --- One experiment already has its own overview; no duplicate line chart.
        if len(curves) > 1:
            fig, ax = plt.subplots(figsize=(13, 5))
            for label, curve, status in curves:
                if status == 'CRITICAL':
                    ax.plot(dates, curve, label=label + ': ALL ZERO (CRITICAL)', lw=1.8)
                else:
                    ax.plot(dates, curve, label=label, lw=1.1)
            ax.axhline(0, color='gray', lw=0.7)
            ax.axvspan(dates[0], dates[11], color='gray', alpha=0.12, label='1999 boundary')
            ax.axvspan(dates[314], dates[-1], color='gray', alpha=0.12,
                       label='Mar-Dec 2025: OBS unavailable')
            if all(s == 'CRITICAL' for _, _, s in curves):
                ax.text(0.5, 0.7, 'ALL PREDICTED TDI CORRECTIONS ARE ZERO',
                        ha='center', transform=ax.transAxes, fontsize=13, weight='bold')
            ax.legend(fontsize=8)
            ax.set_ylabel('Sum of land-region TDI corrections (native units)')
            ax.set_title('D4 regional corrections: ' + self.case)
            fig.tight_layout()
            fig.savefig(self.diag_dir / 'experiment_correction_comparison.png', dpi=160)
            plt.close(fig)
        else:
            (self.diag_dir / 'source_sensitivity.png').unlink(missing_ok=True)
            (self.diag_dir / 'experiment_correction_comparison.png').unlink(missing_ok=True)
        return quality, compared

    def diagnostics(self, dry_run=False):
        runs = list(self.find_runs())
        if not runs:
            raise FileNotFoundError('No ' + self.icase + ' inversion runs in ' + str(self.result_root))
        print('D4 diagnostics:', self.case, '| runs:', len(runs))
        if dry_run:
            for item in runs:
                print('  inspect:', item)
            return
        self.diag_dir.mkdir(parents=True, exist_ok=True)
        # --- Remove only obsolete stage-numbered plots/tables made by older D4.
        # --- Do not touch TDI run folders or user-created diagnostics.
        for file in self.diag_dir.iterdir():
            if file.is_file() and re.fullmatch(r'f[0-9]+_.*\.(png|csv)', file.name, re.I):
                file.unlink()
        rows = [self.diagnose_one(item) for item in runs]
        rows = [row for row in rows if row is not None]
        if not rows:
            raise ValueError('No recognizable D4 experiment directories')
        keys = list(rows[0])
        with (self.diag_dir / 'inversion_runs.csv').open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore')
            w.writeheader()
            w.writerows(rows)
        candidates = [Path(row['path']) / 'obsrvCH4_9925.nc' for row in rows if row['input_observations']]
        obs_info = self.plot_observation_diagnostics(candidates[0]) if candidates else None
        quality, n_compare = self.plot_source_sensitivity(rows)
        n_prior = self.plot_prior_inputs(rows)
        if len(rows) > 1:
            self.plot_experiment_comparison(rows)
        else:
            (self.diag_dir / 'experiment_inputs.png').unlink(missing_ok=True)
        self.write_log_diagnostics(rows)
        n_pass = sum(row['status'] == 'PASS' for row in rows)
        critical = [q for q in quality if q['status_numerical'] == 'CRITICAL']
        undecoded = [row for row in rows if 'unrecognized_predicted_sources_format' in row['problems']]
        lines = [
            'D4 | TDI inversion diagnostic',
            'Case: ' + self.case + '   Configuration: ' + self.icase,
            'Output: ' + str(self.diag_dir),
            'Experiments found: ' + str(len(rows)) +
            '   Run QC passed: ' + str(n_pass) + '/' + str(len(rows)),
            'Decodable regional datasets: ' + str(len(quality)),
            'Unrecognized TDI source files: ' + str(len(undecoded)),
            'NUMERICAL CRITICAL: ' + str(len(critical) + len(undecoded)),
            '',
            'INTERPRETATION WINDOWS',
            '  1999: initial inversion boundary (not used for core statistics)',
            '  2000-2024: primary evaluation interval',
            '  Mar-Dec 2025: observations not yet available (not a processing error)',
            '',
            'INPUT OBSERVATIONS',
        ]
        if obs_info:
            lines += [
                '  Valid D3 station-months: ' + str(obs_info['n_valid']),
                '  Stations with >=1 valid month: ' + str(obs_info['n_nonempty_stations']),
                '  Months with no data: ' + str(obs_info['n_empty_months']),
            ]
        else:
            lines += ['  Observation NetCDF unavailable']
        lines += ['', 'TDI NUMERICAL RESULTS']
        for q in quality:
            label = q['label']
            lines.append('  ' + label + ' : ' + q['status_numerical'])
            lines.append('    Nonzero coefficients: ' + str(q['n_nonzero_coefficients']) + '/' +
                         str(q['n_regions'] * q['n_months']))
            lines.append('    Max |regional correction|: ' + f"{q['max_abs_coefficient']:.5g}" +
                         ' (native TDI units)')
            lines.append('    TDI selected station series: ' + str(q['nobs_tdi']))
            lines.append('    TDI observation equations: ' + str(q['ntotdat_tdi']))
            lines.append('    TDI excluded months (among selected series): ' + str(q['excluded_tdi']))
            lines.append('    Solver switches: DOSVD=' + str(q['dosvd']) +
                         ' NUMSVDVECS=' + str(q['numsvdvecs']))
            if obs_info and q['nobs_tdi'] is not None:
                lines.append('    D3 contains ' + str(obs_info['n_nonempty_stations']) +
                             ' series with data; TDI uses ' + str(q['nobs_tdi']) +
                             '. Selection must be checked using station IDs.')
            if q['flags']:
                lines.append('    Flags: ' + q['flags'])
            regional_csv = self.diag_dir / ('source_regional_' + label + '.csv')
            if regional_csv.is_file() and q['status_numerical'] != 'CRITICAL':
                with regional_csv.open(newline='') as f:
                    regional = list(csv.DictReader(f))
                top = sorted(regional, key=lambda r: float(r['mean_abs_2000_2024']),
                             reverse=True)[:8]
                lines.append('    Top regions (2000-2024): region | mean signed | mean absolute')
                for r in top:
                    lines.append('      %02d | %+9.3f | %9.3f' %
                                 (int(r['region']), float(r['mean_2000_2024']),
                                  float(r['mean_abs_2000_2024'])))
            lines.append('    Figure: inversion_overview_' + self.case + '_' + label + '.png')
        lines += [
            '', 'OTHER CHECKS',
            '  D1 prior uncertainty variants: ' + str(n_prior),
            '  Matched CYC/INCA comparisons: ' + str(n_compare),
            '  Numerical status is separate from file completeness.',
            '  D3 OBS-ACTM residuals are inputs, NOT posterior model residuals.',
            '  predicted.sources corrections are NOT gridded posterior fluxes.',
            '  No TDI executable was run by this diagnostic.',
        ]
        report = '\n'.join(lines) + '\n'
        (self.diag_dir / 'diagnostic_report.txt').write_text(report, encoding='utf-8')
        print('D4 TDI DIAGNOSTICS |', self.case, '|', self.icase)
        print('  Experiments  :', len(rows))
        print('  Run QC passed:', str(n_pass) + '/' + str(len(rows)))
        print('  Numerical critical:', len(critical) + len(undecoded))
        if obs_info:
            print('  D3 valid station-months:', obs_info['n_valid'])
            print('  D3 stations with data :', obs_info['n_nonempty_stations'])
            print('  D3 months without data:', obs_info['n_empty_months'])
        for q in quality:
            print('  ', q['label'], ':', q['status_numerical'],
                  '| nonzero:', str(q['n_nonzero_coefficients']) + '/' +
                  str(q['n_regions'] * q['n_months']),
                  '| max abs:', f"{q['max_abs_coefficient']:.5g}",
                  '| TDI series:', q['nobs_tdi'],
                  '| TDI equations:', q['ntotdat_tdi'])
        print('  Full report  :', self.diag_dir / 'diagnostic_report.txt')
        print('  Overview     :', self.diag_dir / ('inversion_overview_' + self.case + '_*.png'))
        print('  Diagnostics  :', self.diag_dir)
        if critical or undecoded:
            raise RuntimeError('D4 NUMERICAL QUALITY FAILED: zero, nonfinite or unreadable inversion coefficients. '
                               'Inspect diagnostic_report.txt and source_quality.csv; TDI was not rerun.')


# ============================================================
# === Standard workflow entry points
def run_stage(case='INCA', dry_run=False):
    previous = os.environ.get('GCP_CASE')
    os.environ['GCP_CASE'] = str(case).upper()
    try:
        Inversion(case=case).run(dry_run=dry_run)
    finally:
        if previous is None:
            os.environ.pop('GCP_CASE', None)
        else:
            os.environ['GCP_CASE'] = previous


def run_diagnostics_stage(case='INCA', dry_run=False):
    previous = os.environ.get('GCP_CASE')
    os.environ['GCP_CASE'] = str(case).upper()
    try:
        Inversion(case=case).diagnostics(dry_run=dry_run)
    finally:
        if previous is None:
            os.environ.pop('GCP_CASE', None)
        else:
            os.environ['GCP_CASE'] = previous


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=('INCA', 'CYC'), required=True)
    parser.add_argument('--diagnostics-only', action='store_true')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    if args.diagnostics_only:
        run_diagnostics_stage(case=args.case, dry_run=args.dry_run)
    else:
        run_stage(case=args.case, dry_run=args.dry_run)
