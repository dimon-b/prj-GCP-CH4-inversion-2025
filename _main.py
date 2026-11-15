# -*- coding: utf-8 -*-
"""
Project:    GCP/GMB 2025 project
Purpose:    Single entry point for CYC/INCA, without recompiling Fortran.
"""
import os
import sys
import argparse
import datetime
import subprocess
from pathlib import Path

import pandas as pd
import matplotlib as mpl


# ============================================================
# === Workflow switches (same default as the original _main.py)
RUN = {
    'flux': 0,
    'flux_gt3': 0,
    'f0': 0,
    'f2': 0,
    'f3': 0,
    'f4': 0,
    'f5': 0,
    'f6': 0,
    'f7': 0,
    'oh': 0,
    'actm_prior': 0,
    'obs': 0,
    'prepare': 0,
    'filter': 0,
    'inversion': 0,
    'inv_check': 0,
    'diag': 0,
    'server': 0,
    'actm_post': 0,
    'submission': 0,
    'validation': 1,
}

# --- Fortran preparation sequence
#     Names are logical identifiers, NOT assumed executable filenames.
FORTRAN_PREPARE = [
    'F2',  # f2_actm_ch4_at_obs_sites.f90
    'F3',  # f3_combch4_obs_actm.f90
    'F4',  # f4_ch4bu_mon_reg_t42.f90
    'F5',  # f5_write_prior_ch4.f90
    'F6',  # f6_nakazawa_filter.f90 (unchanged Nakazawa libraries)
    'F7',  # f7_write_tdi_obsrvCH4dif.f90
]


# ============================================================
# === Write runtime namelist for F6 (Nakazawa digital filter)
def write_f6_config(hcase, setup, dry_run=False):
    input_file = Path(setup.f6_input_file).expanduser().resolve()
    output_file = Path(setup.f6_output_file).expanduser().resolve()
    config_dir = Path(setup.inv_wrk_dir).expanduser().resolve() / 'config'
    config_file = config_dir / 'actm_f6.nml'

    if hcase.lower() not in input_file.name.lower() or hcase.lower() not in output_file.name.lower():
        raise ValueError('F6 paths do not match case ' + hcase)
    check_namelist_values([hcase, input_file, output_file])
    contents = ("&actm_f6\n"
                " hcase = '" + hcase + "'\n"
                " i_file = '" + str(input_file) + "'\n"
                " o_file = '" + str(output_file) + "'\n"
                "/\n")

    print('\n\tF6 case       :', hcase)
    print('\tF6 input      :', input_file)
    print('\tF6 output     :', output_file)
    print('\tF6 config     :', config_file)

    if not dry_run:
        if not input_file.is_file():
            raise FileNotFoundError(input_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        config_dir.mkdir(parents=True, exist_ok=True)
        config_file.write_text(contents, encoding='utf-8')

    return config_file


# ============================================================
# === Run Nakazawa filter (unmodified lsfbf4_main / lsfbf4_sub)
def run_f6(hcase, setup, dry_run=False):
    config_file = write_f6_config(hcase, setup, dry_run)
    # --- Nakazawa's legacy subroutines write auxiliary files (e.g. test2222.lst).
    #     Isolate them per case without touching the numerical routines.
    work_dir = Path(setup.inv_wrk_dir).expanduser().resolve() / 'f6_work'
    print('\tF6 working dir:', work_dir)
    run_configured_exe('F6', setup, config_file, dry_run, work_dir=work_dir)


# --- Backward-compatible function names for original filter stage
write_filter_config = write_f6_config
run_filter = run_f6


# ============================================================
# === Write runtime configuration for F2 (GT3 -> station .grd)
def write_f2_config(hcase, setup, dry_run=False):
    input_dir = Path(setup.f2_model_dir).expanduser().resolve()
    output_dir = Path(setup.f2_output_dir).expanduser().resolve()
    station_file = Path(setup.sites_f).expanduser().resolve()
    config_dir = Path(setup.inv_wrk_dir).expanduser().resolve() / 'config'
    config_file = config_dir / 'actm_f2.nml'
    output_prefix = setup.f2_output_prefix

    for value in [input_dir, output_dir, station_file, output_prefix, config_file]:
        if "'" in str(value):
            raise ValueError('Single quote in Fortran namelist value: ' + str(value))

    contents = ("&actm_f2\n"
                " hcase = '" + hcase + "'\n"
                " obs_file = '" + str(station_file) + "'\n"
                " i_dir = '" + str(input_dir) + os.sep + "'\n"
                " o_dir = '" + str(output_dir) + os.sep + "'\n"
                " o_fil = '" + output_prefix + "'\n"
                " s_year = " + str(setup.yr_s) + "\n"
                " f_year = " + str(setup.yr_e) + "\n"
                "/\n")

    print('\n\tF2 case        :', hcase)
    print('\tF2 site list   :', station_file)
    print('\tF2 model dir   :', input_dir)
    print('\tF2 output dir  :', output_dir)
    print('\tF2 config      :', config_file)

    if not dry_run:
        if not station_file.is_file():
            raise FileNotFoundError(station_file)
        if not input_dir.is_dir():
            raise FileNotFoundError(input_dir)
        # Protect against silently selecting the wrong case model directory.
        if hcase not in input_dir.name.upper():
            raise RuntimeError('Model directory does not identify ' + hcase + ': ' + str(input_dir))
        output_dir.mkdir(parents=True, exist_ok=True)
        config_dir.mkdir(parents=True, exist_ok=True)
        config_file.write_text(contents, encoding='utf-8')

    return config_file


# ============================================================
# === Run compiled F2 program (no recompilation)
def run_f2(hcase, setup, dry_run=False):
    config_file = write_f2_config(hcase, setup, dry_run)
    # --- The executable is configured in _set_case.py.
    #     CH4_EXE_F2 is an optional fallback for older configurations.
    exe_value = str(getattr(setup, 'exe_f2', '') or '').strip()
    if not exe_value:
        exe_value = os.environ.get('CH4_EXE_F2', '').strip()

    exe = Path(exe_value).expanduser().resolve() if exe_value else None
    print('\tF2 executable  :', exe if exe else '(not configured)')

    if not exe and not dry_run:
        raise RuntimeError('Set self.exe_f2 in _set_case.py or CH4_EXE_F2')

    if dry_run:
        return

    if not exe.is_file():
        raise FileNotFoundError(exe)

    subprocess.run([str(exe), str(config_file)], check=True)


# ============================================================
# === Write runtime configuration for F3 (station ACTM + observations)
def write_f3_config(hcase, setup, dry_run=False):
    input_dir = Path(setup.f3_model_dir).expanduser().resolve()
    obs_dir = Path(setup.f3_obs_dir).expanduser().resolve()
    output_dir = Path(setup.f3_output_dir).expanduser().resolve()
    station_file = Path(setup.sites_f).expanduser().resolve()
    config_dir = Path(setup.inv_wrk_dir).expanduser().resolve() / 'config'
    config_file = config_dir / 'actm_f3.nml'
    input_prefix = setup.f3_input_prefix
    output_prefix = setup.f3_output_prefix

    for value in [input_dir, obs_dir, output_dir, station_file, config_file, input_prefix, output_prefix]:
        if "'" in str(value):
            raise ValueError('Single quote in F3 namelist value: ' + str(value))

    if hcase.lower() not in input_prefix.lower() or hcase.lower() not in output_prefix.lower():
        raise ValueError('F3 filenames do not match case ' + hcase)

    contents = ("&actm_f3\n"
                " hcase = '" + hcase + "'\n"
                " obs_file = '" + str(station_file) + "'\n"
                " obs_dir = '" + str(obs_dir) + os.sep + "'\n"
                " i_dir = '" + str(input_dir) + os.sep + "'\n"
                " o_dir = '" + str(output_dir) + os.sep + "'\n"
                " i_fil = '" + input_prefix + "'\n"
                " o_yr = '" + output_prefix + "ch4yr.grd'\n"
                " o_mm = '" + output_prefix + "ch4mm.grd'\n"
                " o_dy = '" + output_prefix + "ch4dy.grd'\n"
                " s_year = " + str(setup.yr_s) + "\n"
                " f_year = " + str(setup.yr_e) + "\n"
                "/\n")

    print('\n\tF3 case        :', hcase)
    print('\tF3 site list   :', station_file)
    print('\tF3 obs dir     :', obs_dir)
    print('\tF3 model dir   :', input_dir)
    print('\tF3 output dir  :', output_dir)
    print('\tF3 config      :', config_file)

    if not dry_run:
        if not station_file.is_file():
            raise FileNotFoundError(station_file)
        if not obs_dir.is_dir():
            raise FileNotFoundError(obs_dir)
        if not input_dir.is_dir():
            raise FileNotFoundError(input_dir)
        missing = [input_dir / (input_prefix + str(y) + '.grd')
                   for y in range(setup.yr_s, setup.yr_e + 1)
                   if not (input_dir / (input_prefix + str(y) + '.grd')).is_file()]
        if missing:
            raise FileNotFoundError('Missing ' + str(len(missing)) + ' F2 input files; first: ' + str(missing[0]))
        output_dir.mkdir(parents=True, exist_ok=True)
        config_dir.mkdir(parents=True, exist_ok=True)
        config_file.write_text(contents, encoding='utf-8')

    return config_file


# ============================================================
# === Run F3 without recompiling between CYC and INCA
def run_f3(hcase, setup, dry_run=False):
    config_file = write_f3_config(hcase, setup, dry_run)
    exe_value = str(getattr(setup, 'exe_f3', '') or '').strip()
    if not exe_value:
        exe_value = os.environ.get('CH4_EXE_F3', '').strip()
    exe = Path(exe_value).expanduser().resolve() if exe_value else None
    print('\tF3 executable  :', exe if exe else '(not configured)')

    if dry_run:
        return
    if not exe:
        raise RuntimeError('Set self.exe_f3 in _set_case.py or CH4_EXE_F3')
    if not exe.is_file():
        raise FileNotFoundError(exe)
    subprocess.run([str(exe), str(config_file)], check=True)


# ============================================================
# === Paths and executables

def check_namelist_values(values):
    for value in values:
        if "'" in str(value) or '\n' in str(value):
            raise ValueError('Invalid Fortran namelist value: ' + str(value))


def write_case_namelist(stage, hcase, setup, fields, dry_run=False, required_files=None, output_paths=None):
    config_dir = Path(setup.inv_wrk_dir).expanduser().resolve() / 'config'
    config_file = config_dir / ('actm_' + stage.lower() + '.nml')
    values = [hcase] + list(fields.values())
    check_namelist_values(values)
    contents = '&actm_' + stage.lower() + '\n'
    contents += " hcase = '" + hcase + "'\n"
    for name, value in fields.items():
        if name in ['cfg_s_year', 'cfg_f_year']:
            contents += ' ' + name + ' = ' + str(int(value)) + '\n'
        else:
            contents += " " + name + " = '" + str(value) + "'\n"
    contents += '/\n'
    print('\n\t' + stage + ' case       :', hcase)
    for name, value in fields.items():
        print('\t' + stage + ' ' + name.ljust(11) + ':', value)
    print('\t' + stage + ' config     :', config_file)
    if not dry_run:
        for f in required_files or []:
            if not Path(f).is_file():
                raise FileNotFoundError(f)
        for out in output_paths or []:
            Path(out).parent.mkdir(parents=True, exist_ok=True)
        config_dir.mkdir(parents=True, exist_ok=True)
        config_file.write_text(contents, encoding='utf-8')
    return config_file


def run_configured_exe(stage, setup, config_file, dry_run=False, work_dir=None):
    exe_value = str(getattr(setup, 'exe_' + stage.lower(), '') or '').strip()
    if not exe_value:
        exe_value = os.environ.get('CH4_EXE_' + stage, '').strip()
    if stage == 'F6' and not exe_value:
        exe_value = os.environ.get('CH4_EXE_FILTER', '').strip()
    exe = Path(exe_value).expanduser().resolve() if exe_value else None
    print('\t' + stage + ' executable :', exe if exe else '(not configured)')
    if dry_run:
        return
    if exe is None:
        raise RuntimeError('Set self.exe_' + stage.lower() + ' in _set_case.py or CH4_EXE_' + stage)
    if not exe.is_file():
        raise FileNotFoundError(exe)
    if work_dir is not None:
        work_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run([str(exe), str(config_file)], cwd=str(work_dir) if work_dir is not None else None, check=True)


# ============================================================
# === F0: NetCDF prior to GT3 (the same prior for both CYC/INCA)
def write_f0_config(hcase, setup, dry_run=False):
    inp = Path(setup.f0_input_dir).expanduser().resolve()
    out = Path(setup.f0_output_dir).expanduser().resolve()
    input_nc = setup.f0_input_nc
    fields = {'inp_dir': str(inp) + os.sep, 'out_dir': str(out) + os.sep,
              'input_nc': input_nc, 'input_case': setup.f0_input_case,
              'cfg_s_year': setup.yr_s, 'cfg_f_year': setup.yr_e}
    return write_case_namelist('F0', hcase, setup, fields, dry_run, [inp / input_nc], [out/'fch4_gcp25_inv1_t42.gt3'])


def run_f0(hcase, setup, dry_run=False):
    config = write_f0_config(hcase, setup, dry_run)
    run_configured_exe('F0', setup, config, dry_run)


# ============================================================
# === F4: regional monthly prior CH4 fluxes (shared CYC/INCA)
def write_f4_config(hcase, setup, dry_run=False):
    flux = Path(setup.f4_flux_file).expanduser().resolve()
    mask = Path(setup.f4_mask_file).expanduser().resolve()
    out = Path(setup.f4_output_file).expanduser().resolve()
    fields = {'flux_file': str(flux), 'mask_file': str(mask), 'uncr_file': str(out),
              'cfg_s_year': setup.yr_s, 'cfg_f_year': setup.yr_e}
    return write_case_namelist('F4', hcase, setup, fields, dry_run, [flux, mask], [out])


def run_f4(hcase, setup, dry_run=False):
    config = write_f4_config(hcase, setup, dry_run)
    run_configured_exe('F4', setup, config, dry_run)


# ============================================================
# === F5: write NetCDF prior uncertainties (shared CYC/INCA)
def write_f5_config(hcase, setup, dry_run=False):
    inp = Path(setup.f5_input_file).expanduser().resolve()
    prefix = Path(setup.f5_output_prefix).expanduser().resolve()
    fields = {'un_file': str(inp), 'nc_file': str(prefix),
              'cfg_s_year': setup.yr_s, 'cfg_f_year': setup.yr_e}
    return write_case_namelist('F5', hcase, setup, fields, dry_run, [inp], [prefix])


def run_f5(hcase, setup, dry_run=False):
    config = write_f5_config(hcase, setup, dry_run)
    if not dry_run:
        for pc in [30, 50, 99]:
            target = Path(str(Path(setup.f5_output_prefix).expanduser().resolve()) + str(pc) + '.nc')
            if target.exists():
                raise FileExistsError('F5 uses nf_noclobber; existing output: ' + str(target))
    run_configured_exe('F5', setup, config, dry_run)


# ============================================================
# === F7: monthly OBS-ACTM differences and uncertainties
def write_f7_config(hcase, setup, dry_run=False):
    obs = Path(setup.sites_f).expanduser().resolve()
    daily = Path(setup.f7_input_file).expanduser().resolve()
    summary = Path(setup.f7_output_summary).expanduser().resolve()
    prior = Path(setup.f7_output_prior).expanduser().resolve()
    nc = Path(setup.f7_output_nc_prefix).expanduser().resolve()
    if any(hcase.lower() not in path.name.lower() for path in [daily, summary, prior, nc]):
        raise ValueError('F7 input/output paths do not match case ' + hcase)
    fields = {'obs_file': str(obs), 'o_obsrv': str(summary), 'i_infit': str(daily),
              'o_prior': str(prior), 'o_pr_nc': str(nc),
              'cfg_s_year': setup.yr_s, 'cfg_f_year': setup.yr_e}
    return write_case_namelist('F7', hcase, setup, fields, dry_run, [obs, daily], [summary, prior, nc])


def run_f7(hcase, setup, dry_run=False):
    config = write_f7_config(hcase, setup, dry_run)
    if not dry_run:
        # --- F7 uses nf_clobber: prevent overwriting existing case NetCDFs.
        prefix = Path(setup.f7_output_nc_prefix).expanduser().resolve()
        outputs = [Path(str(prefix) + unc + inv + '.nc')
                   for unc in setup.unxcases for inv in setup.invcases]
        outputs += [Path(setup.f7_output_summary).expanduser().resolve(),
                    Path(setup.f7_output_prior).expanduser().resolve()]
        existing = [path for path in outputs if path.exists()]
        if existing:
            raise FileExistsError('F7 would overwrite existing output: ' + str(existing[0]))
    run_configured_exe('F7', setup, config, dry_run)


# ============================================================
# === Run existing Fortran executable without recompilation
def run_fortran(stage, hcase, dry_run=False):

    # --- Every executable and work directory must be explicitly configured
    exe_name = 'CH4_EXE_' + stage
    cwd_name = 'CH4_WORK_' + stage
    exe_value = os.environ.get(exe_name, '').strip()
    cwd_value = os.environ.get(cwd_name, '').strip()

    if not exe_value or not cwd_value:
        raise RuntimeError('Set both ' + exe_name + ' and ' + cwd_name)

    exe = Path(exe_value).expanduser().resolve()
    work_dir = Path(cwd_value).expanduser().resolve()

    # --- Safety: old hard-coded Fortran paths are NOT changed by Python
    if os.environ.get('CH4_FORTRAN_CONFIG_READY', '0') != '1':
        raise RuntimeError('Fortran case configuration is not implemented/verified yet. '
                           'Do not run this stage with only CH4_CASE. '
                           'Add runtime path handling to the Fortran source first, '
                           'then set CH4_FORTRAN_CONFIG_READY=1.')

    if not exe.is_file():
        raise FileNotFoundError(exe)
    if not work_dir.is_dir():
        raise FileNotFoundError(work_dir)

    # --- Contract: adapted Fortran stage reads case from its own runtime config
    #     in this working directory; verify the file when stage is enabled.
    case_file = work_dir / 'actm_config.nml'
    if not case_file.is_file():
        raise FileNotFoundError(case_file)
    case_text = case_file.read_text(encoding='utf-8', errors='replace').upper()
    if "'" + hcase + "'" not in case_text and '"' + hcase + '"' not in case_text:
        raise RuntimeError('Config does not identify case ' + hcase + ': ' + str(case_file))

    cmd = [str(exe)]
    print('\n\tFORTRAN stage :', stage)
    print('\tCase          :', hcase)
    print('\tExecutable    :', exe)
    print('\tWork dir      :', work_dir)

    if not dry_run:
        subprocess.run(cmd, cwd=str(work_dir), check=True)


# ============================================================
# === Run ACTM model: optional existing submission command
def run_actm(stage, hcase, dry_run=False):
    env_name = 'CH4_ACTM_PRIOR_CMD' if stage == 'prior' else 'CH4_ACTM_POST_CMD'
    raise RuntimeError('ACTM submission is not configured. Set up a verified '
                       'case-specific submission workflow before enabling ' + env_name)


# ============================================================
# === Main
def main(hcase='CYC', selected=None, dry_run=False, inv_priors=None, inv_uncs=None, stations=None):

    # --- Set case BEFORE creating any SetCase-derived objects
    hcase = str(hcase).strip().upper()
    if hcase not in ['CYC', 'INCA']:
        raise ValueError('Unknown ACTM case: ' + hcase)
    os.environ['CH4_CASE'] = hcase

    # --- Settings
    pd.set_option('display.expand_frame_repr', False)
    mpl.rcParams['figure.dpi'] = 300

    switches = dict(RUN)
    if selected is not None:
        if selected == ['all']:
            raise ValueError('--steps all includes unconfigured ACTM submissions and is unsafe. Select stages explicitly.')
        unknown = set(selected) - set(switches)
        if unknown:
            raise ValueError('Unknown stages: ' + ', '.join(sorted(unknown)))
        # --- Aliases and grouped stages are expanded only once.
        normalized = set(selected)
        if 'flux_gt3' in normalized:
            normalized.discard('flux_gt3')
            normalized.add('f0')
        if 'filter' in normalized:
            normalized.discard('filter')
            normalized.add('f6')
        if 'prepare' in normalized:
            normalized.discard('prepare')
            normalized.update(('f2', 'f3', 'f4', 'f5', 'f6', 'f7'))
        switches = {key: int(key in normalized) for key in switches}

    # --- Check case consistency in SetCase
    import _set_case
    setup = _set_case.SetCase()
    assert setup.hcase == hcase

    print('\n\tACTM case        :', hcase)
    print('\tInversion output :', setup.inv_wrk_dir)
    print('\tStation list     :', setup.sites_f)
    print('\tEnabled stages   :', ', '.join(key for key, val in switches.items() if val))

    # ============================================================
    # === Convert input fluxes nc -> nc
    if switches['flux']:
        print('\n\t *-*-* Start c1_conflux_nc *-*-* ')
        if not dry_run:
            import c1_conflux_nc
            c1_conflux_nc.CnvFluxNc()

    # ============================================================
    # === Convert input fluxes nc -> gt3
    if switches['flux_gt3'] or switches['f0']:
        run_f0(hcase, setup, dry_run)

    # ============================================================
    # === INCA OH: CYC uses existing climatological OH
    if switches['oh']:
        if hcase == 'INCA':
            print('\n\t *-*-* Start i1_inca_OH *-*-* ')
            if not dry_run:
                import b1_inca_OH
                b1_inca_OH.CnvOHNc()
        else:
            print('\n\tSkip INCA OH conversion for CYC')

    # ============================================================
    # === Forward MIROC4-ACTM with prior
    if switches['actm_prior']:
        run_actm('prior', hcase, dry_run)

    # ============================================================
    # === Copy/Create observation files (shared stations in stage 1)
    if switches['obs']:
        print('\n\t *-*-* Start c3_obs_files *-*-* ')
        if not dry_run:
            import c3_obs_files
            c3_obs_files.ObsFilesTxt()

    # ============================================================
    # === Extract ACTM station concentrations (GT3 -> GRD)
    if switches['f2']:
        run_f2(hcase, setup, dry_run)

    # ============================================================
    # === Combine model and observations at stations
    if switches['f3']:
        run_f3(hcase, setup, dry_run)

    if switches['f4']:
        run_f4(hcase, setup, dry_run)

    if switches['f5']:
        run_f5(hcase, setup, dry_run)

    if switches['f6']:
        run_f6(hcase, setup, dry_run)

    if switches['f7']:
        run_f7(hcase, setup, dry_run)

    # ============================================================
    # === Prepare inversion inputs (existing Fortran programs)
    if switches['prepare']:
        for stage in FORTRAN_PREPARE:
            if stage == 'F2':
                run_f2(hcase, setup, dry_run)
            elif stage == 'F3':
                run_f3(hcase, setup, dry_run)
            elif stage == 'F4':
                run_f4(hcase, setup, dry_run)
            elif stage == 'F5':
                run_f5(hcase, setup, dry_run)
            elif stage == 'F6':
                run_f6(hcase, setup, dry_run)
            elif stage == 'F7':
                run_f7(hcase, setup, dry_run)
            else:
                run_fortran(stage, hcase, dry_run)

    # ============================================================
    # === Standalone daily observation/model filtering
    if switches['filter']:
        run_filter(hcase, setup, dry_run)

    # ============================================================
    # === Check case-specific F5/F7 NetCDF inputs
    if switches['inv_check']:
        from inversion.inversion_runner import check_case_inversion
        if not dry_run:
            check_case_inversion(hcase, setup, priors=inv_priors, uncs=inv_uncs)
        else:
            print('\n\tInput check skipped in --dry-run (needs existing NetCDF)')

    # ============================================================
    # === Read-only diagnostics of F6 and F7
    if switches['diag']:
        if not dry_run:
            import c4_diag_ch4
            c4_diag_ch4.run_diagnostics(hcase, setup, stations=stations)
        else:
            print('\n\tF6 input        :', setup.f6_input_file)
            print('\tF6 output       :', setup.f6_output_file)
            print('\tF7 NetCDF prefix:', setup.f7_output_nc_prefix)
            print('\tReport directory:', os.path.join(setup.inv_wrk_dir, 'diagnostics'))

    # ============================================================
    # === Bayesian inversion, old compiled executable
    if switches['inversion']:
        from inversion.inversion_runner import run_inversion
        run_inversion(hcase, setup, dry_run=dry_run, priors=inv_priors, uncs=inv_uncs)

    # ============================================================
    # === Loss/burden and conversion outputs (existing Python wrappers)
    if switches['server']:
        print('\n\t *-*-* Start Server *-*-* ')
        if not dry_run:
            import s0_server
            s0_server.ServerPart()

    # ============================================================
    # === Forward MIROC4-ACTM with posterior
    if switches['actm_post']:
        run_actm('post', hcase, dry_run)

    # ============================================================
    # === Submission
    if switches['submission']:
        print('\n\t *-*-* Start Submission *-*-* ')
        if not dry_run:
            import t0_submit
            t0_submit.Write2Submit()

    # ============================================================
    # === Validation using existing code
    if switches['validation']:
        print('\n\t *-*-* Start Validation *-*-* ')
        if not dry_run:
            import v3_check_hippo
            v3_check_hippo.CheckHippo()
            # import v1_check_3sites
            # v1_check_3sites.Check3Sites()
            # import v2_check_all_surf
            # v2_check_all_surf.CheckAllSurf()
            # import v4_check_GR
            # v4_check_GR.CheckGR()

    print('\nScript done!')


# ============================================================
# === Timing
def est_time(led):
    if led == 0:
        led = datetime.datetime.now()
        print('\nStart time       :', led, '\n')
        return led

    print('\nEnd time         :', datetime.datetime.now())
    print('Running time     :', datetime.datetime.now() - led)
    print('Main done')
    return led


# ============================================================
# === Command line
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='GCP/GMB2025 CYC/INCA workflow')
    parser.add_argument('--case', choices=['CYC', 'INCA'], default=os.environ.get('CH4_CASE', 'CYC').upper())
    parser.add_argument('--steps', nargs='+', choices=list(RUN) + ['all'], default=None)
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--stations', nargs=2, type=int, metavar=('START', 'END'), default=None, help='Inclusive station ID range for --steps diag')
    parser.add_argument('--inv-prior', nargs='+', choices=['p30', 'p50', 'p99', 'all'], default=['p30'])
    parser.add_argument('--inv-unc', nargs='+', choices=['ctl', 'ux2', 'ux4', 'all'], default=['ctl'])
    args = parser.parse_args()

    led = est_time(0)
    inv_priors = ['p30', 'p50', 'p99'] if 'all' in args.inv_prior else args.inv_prior
    inv_uncs = ['ctl', 'ux2', 'ux4'] if 'all' in args.inv_unc else args.inv_unc
    main(hcase=args.case, selected=args.steps, dry_run=args.dry_run, inv_priors=inv_priors, inv_uncs=inv_uncs, stations=args.stations)
    est_time(led)

