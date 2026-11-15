# -*- coding: utf-8 -*-
"""Prepare and run TDI for CYC/INCA; detect the compacted valid-month reader.

Observation equations are compacted by the updated Fortran driver.
The Bayesian solver and its NetCDF inputs are unchanged.
"""
import hashlib
import os
import re
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# === The compacted-data Fortran reader must be positively identified

def supports_valid_months(executable):
    """Probe the executable in an empty directory so old binaries cannot run a job."""
    executable = Path(executable).expanduser().resolve() if executable else None
    if executable is None or not executable.is_file():
        return False
    with tempfile.TemporaryDirectory() as tmp:
        try:
            result = subprocess.run([str(executable), '--valid-months-version'],
                                    cwd=tmp, capture_output=True, text=True,
                                    timeout=5, check=False)
        except (OSError, subprocess.TimeoutExpired):
            return False
    return result.returncode == 0 and 'TDI_VALID_MONTHS_V1' in result.stdout


# ============================================================
# === Immutable input snapshots and manifest (no source symlinks)

def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def copy_checked(source, dest):
    # --- Verify that a file did not change while its private copy was made.
    source = Path(source)
    dest = Path(dest)
    if not source.is_file():
        raise FileNotFoundError(source)
    old_stat = source.stat()
    old_hash = sha256_file(source)
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, dest)
    new_stat = source.stat()
    new_hash = sha256_file(source)
    copy_hash = sha256_file(dest)
    if (old_stat.st_size != new_stat.st_size or
            old_stat.st_mtime_ns != new_stat.st_mtime_ns or
            old_hash != new_hash or copy_hash != old_hash):
        raise RuntimeError('Input changed during snapshot: ' + str(source))
    return copy_hash


def copy_master_files(source_dir, dest_dir):
    # --- Dereference symlinks: a run must not depend on mutable files4run links.
    source_dir = Path(source_dir)
    dest_dir = Path(dest_dir)
    if not source_dir.is_dir():
        raise FileNotFoundError('Inversion master directory: ' + str(source_dir))
    dest_dir.mkdir(parents=True, exist_ok=False)
    checksums = {}
    for src_root, dirs, files in os.walk(source_dir, followlinks=False):
        current = Path(src_root)
        relative = current.relative_to(source_dir)
        for dirname in dirs:
            if (current / dirname).is_symlink():
                raise ValueError('Symlink directory inside files4run is not supported: ' + str(current / dirname))
            (dest_dir / relative / dirname).mkdir(exist_ok=True)
        for filename in files:
            src = current / filename
            if not src.is_file():
                raise ValueError('Not a regular file inside files4run: ' + str(src))
            dst = dest_dir / relative / filename
            checksums[dst.relative_to(dest_dir.parent).as_posix()] = copy_checked(src, dst)
    if not checksums:
        raise ValueError('No files found in ' + str(source_dir))
    return checksums


def control_setting(contents, key):
    match = re.search(r"(?im)^\s*" + re.escape(key) + r"\s*=\s*['\"]([^'\"\r\n]+)['\"]", contents)
    if match is None:
        raise ValueError('Missing path in control: ' + key)
    return match.group(1)


def check_master_settings(contents):
    # --- All non-GF master inputs must come from the copied files4run folder.
    keys = ('srccorf', 'outputf', 'groupf', 'constf', 'recipef')
    for key in keys:
        value = control_setting(contents, key)
        path = Path(value)
        if path.is_absolute() or '..' in path.parts or not path.parts or path.parts[0] != 'files4run':
            raise ValueError('Master input must be inside local files4run: ' + key + '=' + value)


def gf_inventory(path):
    # --- External GF is NOT copied. Record its file metadata, not GBs of data.
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError('Green functions path: ' + str(path))
    if path.is_file():
        paths = [path]
        base = path.parent
    else:
        paths = sorted(p for p in path.rglob('*') if p.is_file())
        base = path
    digest = hashlib.sha256()
    total_bytes = 0
    for item in paths:
        stat = item.stat()
        total_bytes += stat.st_size
        line = '%s\t%d\t%d\n' % (item.relative_to(base).as_posix(), stat.st_size, stat.st_mtime_ns)
        digest.update(line.encode('utf-8'))
    return {'count': len(paths), 'bytes': total_bytes, 'metadata_sha256': digest.hexdigest()}


def write_run_manifest(run_dir, hcase, prior, unc, model, paths, checksums, gf_path, gf_info):
    run_dir = Path(run_dir)
    lines = ['# Immutable inputs for this inversion',
             'created_utc=' + datetime.now(timezone.utc).isoformat(),
             'case=' + hcase, 'icase=' + paths['icase'],
             'prior=' + prior, 'unc=' + unc, 'model=' + model,
             'source_F5=' + str(paths['F5']),
             'source_F7=' + str(paths['F7']),
             'source_sites=' + str(paths['sites']),
             'source_control_template=' + str(paths['control']),
             'source_files4run=' + str(paths['files4run']),
             'source_executable=' + str(paths['executable']),
             'vardump_enabled=false',
             'GF_mode=SHARED_NOT_COPIED',
             'GF_path=' + str(gf_path),
             'GF_file_count=' + str(gf_info['count']),
             'GF_total_bytes=' + str(gf_info['bytes']),
             'GF_metadata_sha256=' + gf_info['metadata_sha256'],
             'GF_warning=metadata fingerprint only; preserve GF directory unchanged',
             '', '# SHA256 hashes of COPIED local inputs']
    checksum_lines = []
    for name, checksum in sorted(checksums.items()):
        lines.append(checksum + '  ' + name)
        checksum_lines.append(checksum + '  ' + name)
    (run_dir / 'RUN_INPUTS.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    # --- sha256sum -c RUN_INPUTS.sha256, executed in the run directory.
    (run_dir / 'RUN_INPUTS.sha256').write_text('\n'.join(checksum_lines) + '\n', encoding='utf-8')


# ============================================================
# === Input NetCDF name produced by F7

def get_inversion_input(setup, hcase, unc='ctl', model='inv1'):
    if hcase not in ('CYC', 'INCA'):
        raise ValueError('Unknown ACTM case: ' + str(hcase))
    if unc not in ('ctl', 'ux2', 'ux4'):
        raise ValueError('Unknown data uncertainty: ' + str(unc))
    if model not in tuple(setup.invcases):
        raise ValueError('Unknown inversion model: ' + str(model))

    prefix = Path(setup.f7_output_nc_prefix).expanduser().resolve()
    if hcase.lower() not in prefix.name.lower():
        raise ValueError('F7 prefix does not match the selected case: ' + str(prefix))
    return Path(str(prefix) + unc + model + '.nc')


# ============================================================
# === Reconfigure only paths (all physical inversion switches unchanged)

def get_control_text(template, hcase, site_file, first_year, last_year):
    contents = template.read_text(encoding='utf-8')
    if '&control' not in contents.lower():
        raise ValueError('Not a valid control namelist: ' + str(template))

    settings = {'sourcef': 'sourcepriortdi.nc',
                'datf': 'obsrvCH4_9925.nc',
                'datflist': str(site_file)}

    for key, value in settings.items():
        if "'" in value or '\n' in value or '\r' in value:
            raise ValueError('Unsafe character in control value: ' + key)
        if len(value) > 80:
            raise ValueError(key + ' exceeds 80 characters in legacy Fortran: ' + value)
        pattern = re.compile(r"(?im)^(\s*" + re.escape(key) + r"\s*=\s*)(['\"])[^'\"\r\n]*\2")
        contents, count = pattern.subn(lambda match: match.group(1) + "'" + value + "'", contents)
        if count != 1:
            raise ValueError('Expected exactly one control key: ' + key)

    # --- Disable the enormous binary covariance dump for every new inversion.
    #     Apply this to the generated run control even when an older template
    #     still contains 'dumpmats = .true.'. The posterior solution is unaffected.
    pattern = re.compile(r'(?im)^(\s*dumpmats\s*=\s*)(\.?\w+\.?)\s*(?=$|[!#])')
    contents, count = pattern.subn(lambda match: match.group(1) + '.false.', contents)
    if count != 1:
        raise ValueError('Expected exactly one dumpmats switch in control')

    for key, expected in [('firstsrc', first_year), ('lastsrc', last_year),
                          ('firstdat', first_year), ('lastdat', last_year)]:
        match = re.search(r'(?im)^\s*' + key + r'\s*=\s*(\d+)\b', contents)
        if match is None or int(match.group(1)) != expected:
            raise ValueError('Inversion control has wrong ' + key + ', expected ' + str(expected))

    for key, expected in [('tc3tdi', True), ('pseudodata', False)]:
        match = re.search(r'(?im)^\s*' + key + r'\s*=\s*(\.?\w+\.?)', contents)
        if match is None:
            raise ValueError('Missing control switch: ' + key)
        value = match.group(1).strip('.').lower()
        if (value in ('t', 'true')) != expected:
            raise ValueError('Unexpected inversion control switch: ' + key)

    # --- Legacy character*80 applies to all file paths in control.
    keys = ('sourcef', 'datf', 'datflist', 'srccorf', 'outputf', 'groupf',
            'constf', 'recipef', 'greenf')
    for key in keys:
        match = re.search(r"(?im)^\s*" + key + r"\s*=\s*['\"]([^'\"\r\n]*)['\"]", contents)
        if match is None:
            raise ValueError('Missing control file path: ' + key)
        if len(match.group(1)) > 80:
            raise ValueError(key + ' exceeds 80 characters: ' + match.group(1))
        if not match.group(1):
            raise ValueError(key + ' is empty')
    return contents


# ============================================================
# === Build isolated run directory compatible with legacy results

def get_inversion_paths(setup, hcase, prior='p30', unc='ctl', model='inv1'):
    if prior not in tuple(setup.unpcases):
        raise ValueError('Unknown prior uncertainty: ' + str(prior))
    nc_obs = get_inversion_input(setup, hcase, unc, model)
    source_prefix = Path(setup.f5_output_prefix).expanduser().resolve()
    nc_source = Path(str(source_prefix) + prior.removeprefix('p') + '.nc')
    root = Path(setup.inv_wrk_dir).expanduser().resolve()
    run_dir = root / (setup.icase + '_' + prior) / (setup.icase + '_' + unc + '_' + model)
    return nc_obs, nc_source, run_dir


def check_inversion_inputs(obs_file, source_file, sites_file, allow_missing=False):
    from .check_inv_inputs import inspect_observations, inspect_sources
    obs = inspect_observations(obs_file, sites_file, strict=not allow_missing)
    src = inspect_sources(source_file)
    if not sites_file.is_file():
        raise FileNotFoundError(sites_file)
    if obs['n_stations'] < 61 or obs['n_months'] != 324:
        raise ValueError('Unexpected F7 observations: ' + str(obs))
    if src['n_regions'] != 54:
        raise ValueError('Unexpected F5 source regions: ' + str(src))
    return obs, src


# ============================================================
# === Validate selected F7/F5 input files without running TDI

def check_case_inversion(hcase, setup, priors=None, uncs=None, models=None):
    if hcase != setup.hcase:
        raise ValueError('SetCase hcase does not match: ' + hcase)
    priors = list(priors) if priors is not None else ['p30']
    uncs = list(uncs) if uncs is not None else ['ctl']
    models = list(models) if models is not None else ['inv1']
    sites = Path(setup.sites_f).expanduser().resolve()
    allow_missing = supports_valid_months(getattr(setup, 'exe_inv', None))
    print('\tINV compaction enabled:', allow_missing)
    for prior in priors:
        for model in models:
            for unc in uncs:
                obs, source, work = get_inversion_paths(setup, hcase, prior, unc, model)
                print('\n\tCase / prior / unc :', hcase, prior, unc)
                print('\tF7 NetCDF          :', obs)
                print('\tF5 NetCDF          :', source)
                obs_info, source_info = check_inversion_inputs(obs, source, sites, allow_missing=allow_missing)
                print('\tActive stations    :', obs_info['n_active_stations'])
                print('\tInvalid active data:', obs_info['n_invalid_active'])
                print('\tSource regions     :', source_info['n_regions'])
    print('\nInversion input checks: PASS')


# ============================================================
# === Run selected Bayesian inversion using local tdi_ch4

def run_inversion(hcase, setup, dry_run=False, priors=None, uncs=None, models=None, check_nc=True):
    if hcase != setup.hcase:
        raise ValueError('SetCase hcase does not match: ' + hcase)
    priors = list(priors) if priors is not None else ['p30']
    uncs = list(uncs) if uncs is not None else ['ctl']
    models = list(models) if models is not None else ['inv1']
    template = Path(setup.inv_control_template).expanduser().resolve()
    executable = Path(setup.exe_inv).expanduser().resolve() if setup.exe_inv else None
    sites = Path(setup.sites_f).expanduser().resolve()

    print('\n\tINV case           :', hcase)
    print('\tINV executable     :', executable or '(not set)')
    print('\tINV control source :', template)

    # --- Check configuration even for dry-run when template is accessible.
    if template.is_file():
        control_text = get_control_text(template, hcase, 'obsrvCH4_sites.txt', setup.yr_s, setup.yr_e)
        check_master_settings(control_text)
    elif dry_run:
        print('\tINV control        : template missing; content not checked')
        control_text = None
    else:
        raise FileNotFoundError(template)

    if not dry_run:
        if executable is None or not executable.is_file():
            raise FileNotFoundError('tdi_ch4 executable: ' + str(executable))
        if not os.access(executable, os.X_OK):
            raise PermissionError('tdi_ch4 is not executable: ' + str(executable))
        if control_text is None:
            raise FileNotFoundError(template)
        if not sites.is_file():
            raise FileNotFoundError(sites)

    allow_missing = supports_valid_months(executable) if not dry_run else False
    print('\tINV compaction enabled:', allow_missing)

    for prior in priors:
        for model in models:
            for unc in uncs:
                obs, source, run_dir = get_inversion_paths(setup, hcase, prior, unc, model)
                print('\n\tINV prior          :', prior)
                print('\tINV model          :', model)
                print('\tINV uncertainty    :', unc)
                print('\tINV observations   :', obs)
                print('\tINV source prior   :', source)
                print('\tINV work dir       :', run_dir)
                if dry_run:
                    continue

                if not obs.is_file():
                    raise FileNotFoundError(obs)
                if not source.is_file():
                    raise FileNotFoundError(source)
                if run_dir.exists():
                    raise FileExistsError('Inversion directory already exists. Refusing overwrite: ' + str(run_dir))
                if check_nc:
                    check_inversion_inputs(obs, source, sites, allow_missing=allow_missing)

                # --- New compacting binary: validate auxiliary master files
                #     before creating a results directory. Legacy/mocked
                #     executables preserve their previous behavior.
                if allow_missing:
                    for key in ('srccorf', 'outputf', 'groupf', 'constf', 'recipef'):
                        match = re.search(r"(?im)^\s*" + re.escape(key) + r"\s*=\s*['\"]([^'\"]+)['\"]", control_text)
                        if match is None:
                            raise ValueError('Missing inversion master file key: ' + key)
                        master = Path(match.group(1))
                        if not master.is_absolute():
                            master = template.parent / master
                        if not master.is_file():
                            raise FileNotFoundError('Missing inversion master file ' + key + ': ' + str(master))
                    match = re.search(r"(?im)^\s*greenf\s*=\s*['\"]([^'\"]+)['\"]", control_text)
                    if match is None or not Path(match.group(1)).exists():
                        raise FileNotFoundError('Missing inversion GF path: ' + (match.group(1) if match else '(unset)'))

                # --- Shared Green's functions remain external; record their metadata.
                gf_path = Path(control_setting(control_text, 'greenf'))
                if not gf_path.is_absolute():
                    raise ValueError('greenf must use an absolute shared GF path: ' + str(gf_path))
                gf_info = gf_inventory(gf_path)
                static_dir = template.parent / 'files4run'
                if not static_dir.is_dir():
                    raise FileNotFoundError('Missing inversion files4run: ' + str(static_dir))
                for key in ('srccorf', 'outputf', 'groupf', 'constf', 'recipef'):
                    master = Path(control_setting(control_text, key))
                    if not (static_dir.parent / master).is_file():
                        raise FileNotFoundError('Missing ' + key + ': ' + str(static_dir.parent / master))

                run_dir.parent.mkdir(parents=True, exist_ok=True)
                run_dir.mkdir(exist_ok=False)
                # --- All files below are real copies. Do not follow source changes.
                checksums = {}
                try:
                    checksums['obsrvCH4_9925.nc'] = copy_checked(obs, run_dir / 'obsrvCH4_9925.nc')
                    checksums['sourcepriortdi.nc'] = copy_checked(source, run_dir / 'sourcepriortdi.nc')
                    checksums['obsrvCH4_sites.txt'] = copy_checked(sites, run_dir / 'obsrvCH4_sites.txt')
                    checksums.update(copy_master_files(static_dir, run_dir / 'files4run'))
                    checksums['tdi_ch4'] = copy_checked(executable, run_dir / 'tdi_ch4')
                    # --- control uses only the frozen local filenames for data/sites.
                    (run_dir / 'control').write_text(control_text, encoding='utf-8')
                    checksums['control'] = sha256_file(run_dir / 'control')
                    write_run_manifest(run_dir, hcase, prior, unc, model,
                                       {'icase': setup.icase, 'F5': source, 'F7': obs,
                                        'sites': sites, 'control': template,
                                        'files4run': static_dir, 'executable': executable},
                                       checksums, gf_path, gf_info)
                    print('\tINV snapshot       :', run_dir / 'RUN_INPUTS.txt')
                except Exception as exc:
                    (run_dir / 'SNAPSHOT_FAILED.txt').write_text(str(exc) + '\n', encoding='utf-8')
                    raise

                logfile = run_dir / ('tdi_' + prior + '_' + model + '_' + unc + '.out')
                print('\tINV log            :', logfile)
                with logfile.open('w', encoding='utf-8') as log:
                    subprocess.run([str(run_dir / 'tdi_ch4')], cwd=run_dir,
                                   stdout=log, stderr=subprocess.STDOUT, check=True)
                # --- Warn if shared Green's functions changed during calculation.
                if gf_inventory(gf_path) != gf_info:
                    raise RuntimeError('Shared Green functions changed during run: ' + str(gf_path))
                # --- Some legacy STOP statements return exit code 0 even on failure.
                for artifact in ('prior.sources', 'predicted.sources'):
                    item = run_dir / artifact
                    if not item.is_file() or item.stat().st_size == 0:
                        raise RuntimeError('TDI did not produce ' + artifact + '; inspect ' + str(logfile))
                print('\tINV finished       :', run_dir)
