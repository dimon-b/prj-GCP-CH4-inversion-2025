# -*- coding: utf-8 -*-
"""Compile post-processing Fortran from a temporary patched copy."""

import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


# --- Project directories (do not depend on shell's current directory)
PROJECT_DIR = Path(__file__).resolve().parent
ROOT_DIR = PROJECT_DIR.parent
POST_SOURCE_DIR = PROJECT_DIR / 'fortran'


def full_path(value):
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = PROJECT_DIR / path
    return path.resolve()


def fortran_path(value, directory=False):
    path = full_path(value)
    if directory:
        return str(path) + os.sep
    return str(path)


def replace_parameter(source, name, value):
    # --- Replace a single named Fortran PARAMETER declaration, not a line number
    pattern = re.compile(
        r'(?m)^(\s*(?:character\s*\(\s*len\s*=\s*(?:\*|\d+)\s*\)|integer|logical)\s*,'
        r'\s*parameter\s*::\s*' + re.escape(name) + r'(?:\s*\([^)]*\))?\s*=)[^\n]*$',
        re.IGNORECASE,
    )
    def format_replacement(match):
        prefix = match.group(1) + ' '
        line = prefix + value
        if len(line) <= 120:
            return line
        if value.startswith("'") and value.endswith("'") and "'" not in value[1:-1]:
            chunks = [value[1:-1][i:i + 52] for i in range(0, len(value) - 2, 52)]
            if not chunks:
                return line
            return prefix + "'" + chunks[0] + "' // &\n" + "\n".join(
                "        & '" + chunk + "'" + (" // &" if i < len(chunks) - 2 else "")
                for i, chunk in enumerate(chunks[1:])
            )
        raise ValueError(f'Fortran PARAMETER {name} is too long; refusing a truncated line')

    source, count = pattern.subn(format_replacement, source)
    if count != 1:
        raise ValueError(f'Fortran PARAMETER {name}: expected 1 declaration, got {count}')
    return source


def replace_years(source, start, end):
    if end < start:
        raise ValueError('Invalid year interval')
    value = f'{start}, f_year = {end}, nyear = f_year-s_year+1'
    return replace_parameter(source, 's_year', value)


def replace_cases(source, model, prior, uncertainty):
    # --- These legacy conversion programs explicitly support one model at a time
    if model != 'inv1':
        raise ValueError('p3/p4 support inv1 only: update the Fortran recipe for other models')
    if not re.fullmatch(r'[A-Za-z0-9_]+', prior):
        raise ValueError('Invalid prior name')
    if not re.fullmatch(r'[A-Za-z0-9_]+', uncertainty):
        raise ValueError('Invalid uncertainty name')

    source = replace_parameter(source, 'ninv', '1, nunp = 1, nunx = 1')
    source = replace_parameter(source, 'invCAses', "['inv1']")
    source = replace_parameter(source, 'flxcases', "['fch4_gcp2025_prior_inv1']")
    source = replace_parameter(source, 'unpcases', f"['{prior}']")
    source = replace_parameter(source, 'unxcases', f"['{uncertainty}']")
    return source


def prepare_fortran(filename, setup, prior='p30', uncertainty='ctl', model='inv1', model_dir=None, lc_mode='file'):
    if lc_mode not in ('file', 'identity'):
        raise ValueError('LC mode must be file or identity')
    source_file = POST_SOURCE_DIR / filename
    if not source_file.is_file():
        raise FileNotFoundError(f'Fortran source missing: {source_file}')

    content = source_file.read_text(encoding='utf-8')
    content = replace_years(content, int(setup.yr_s), int(setup.yr_e))
    content = replace_parameter(content, 'work_dir', f"'{fortran_path(setup.inv_wrk_dir, directory=True)}'")

    if filename == 'pf1_ch4_burloss.f90':
        content = replace_parameter(content, 'burd_file', f"'{fortran_path(setup.inv_lsc_dir, directory=True)}frt_burd_{setup.icase}_{prior}_{uncertainty}_'")
        model_dir = model_dir or getattr(setup, 'inv_mod_dir', None)
        if not model_dir:
            raise ValueError('p1: set inv_mod_dir in _set_case.py or provide --model-dir')
        content = replace_parameter(content, 'ch4_dir', f"'{fortran_path(model_dir, directory=True)}'")
    elif filename in ('pf3_ch4_grd.f90', 'pf4_ch4_grd2gt3.f90'):
        content = replace_parameter(content, 'icase', f"'{setup.icase}'")
        content = replace_cases(content, model, prior, uncertainty)
        content = replace_parameter(content, 'flux_suffix', "'_raw'" if lc_mode == 'identity' else "''")
        if filename == 'pf3_ch4_grd.f90':
            content = replace_parameter(content, 'apply_lc', '.false.' if lc_mode == 'identity' else '.true.')
            content = replace_parameter(content, 'lsc_file', f"'gcp_LC_{setup.icase}_{prior}_{uncertainty}_'")
            content = replace_parameter(content, 'cgt_file', f"'check_GT_{setup.icase}_raw_'" if lc_mode == 'identity' else f"'check_GT_{setup.icase}_'")
            content = replace_parameter(content, 'flx_dir', f"'{fortran_path(setup.flx_apr_dir, directory=True)}'")
            content = replace_parameter(content, 'msk_file', f"'{fortran_path(Path(setup.inv_dir) / 'mask/mask80.bin')}'")
            content = replace_parameter(content, 'lnb_file', f"'{fortran_path(Path(setup.inv_dir) / 'mask/ch4_regn_pulse.grd')}'")
    else:
        raise ValueError(f'Unexpected Fortran file: {filename}')
    return content


def select_compiler():
    # --- Prefer the Intel compiler used by the production server
    requested = os.environ.get('GCP_FC', '').strip()
    if requested:
        if not shutil.which(requested):
            raise FileNotFoundError(f'GCP_FC compiler not found: {requested}')
        return requested
    for name in ['ifort', 'ifx']:
        if shutil.which(name):
            return name
    raise FileNotFoundError('Intel Fortran compiler not found (ifort/ifx). Load compiler module, or set GCP_FC=gfortran explicitly for testing.')


def check_inputs(paths):
    for filename in paths:
        path = full_path(filename)
        if not path.is_file():
            raise FileNotFoundError(f'Missing post-processing input: {path}')


def check_outputs(paths, overwrite=False):
    if overwrite:
        return
    existing = [str(full_path(name)) for name in paths if full_path(name).exists()]
    if existing:
        raise FileExistsError('Output exists (use --overwrite):\n' + '\n'.join(existing))


def run_fortran(filename, setup, prior='p30', uncertainty='ctl', model='inv1', model_dir=None,
                inputs=(), outputs=(), output_dirs=(), dry_run=False, overwrite=False, lc_mode='file'):
    content = prepare_fortran(filename, setup, prior=prior, uncertainty=uncertainty,
                              model=model, model_dir=model_dir, lc_mode=lc_mode)
    print('\n    POST Fortran        :', filename)
    print('    POST case           :', setup.hcase, setup.icase)
    print('    POST prior/unc      :', prior, uncertainty)
    print('    POST LC mode        :', lc_mode)
    print('    POST working dir    :', fortran_path(setup.inv_wrk_dir, directory=True))
    print('    POST inputs         :', len(inputs))
    print('    POST outputs        :', len(outputs))

    if dry_run:
        for name in inputs:
            print('    INPUT               :', full_path(name))
        for name in outputs:
            print('    OUTPUT              :', full_path(name))
        print('    DRY-RUN (no compilation or files written)')
        return

    check_inputs(inputs)
    check_outputs(outputs, overwrite=overwrite)
    for folder in output_dirs:
        full_path(folder).mkdir(parents=True, exist_ok=True)

    compiler = select_compiler()
    with tempfile.TemporaryDirectory(prefix='gcp_post_') as tmp:
        tmp_path = Path(tmp)
        src_path = tmp_path / filename
        src_path.write_text(content, encoding='utf-8')
        exe_path = tmp_path / filename.removesuffix('.f90')
        command = [compiler, '-O3', str(src_path), '-o', str(exe_path)]
        if Path(compiler).name == 'gfortran':
            command.insert(2, '-ffree-line-length-none')
        print('    POST compile        :', ' '.join(command))
        compile_result = subprocess.run(command, cwd=tmp_path, text=True, capture_output=True)
        if compile_result.returncode:
            raise RuntimeError('Fortran compilation failed:\n' + compile_result.stdout + compile_result.stderr)
        if compile_result.stderr.strip():
            print(compile_result.stderr)
        # --- Keep working-directory-relative auxiliary paths intact
        logfile = full_path(setup.inv_wrk_dir) / (filename.removesuffix('.f90') + '.log')
        logfile.parent.mkdir(parents=True, exist_ok=True)
        run_result = subprocess.run([str(exe_path)], cwd=PROJECT_DIR, text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        logfile.write_text(run_result.stdout, encoding='utf-8')
        print('    POST log            :', logfile)
        print(run_result.stdout[-5000:])
        if run_result.returncode or re.search(r'(?m)^\s*(?:ERROR\s+)?STOP(?:\s|$)', run_result.stdout):
            raise RuntimeError(f'Post-processing failed: {filename}, exit={run_result.returncode} (see {logfile})')
        for name in outputs:
            if not full_path(name).is_file():
                raise RuntimeError(f'Fortran completed without expected output: {full_path(name)}')
        print('    POST finished       :', filename)
