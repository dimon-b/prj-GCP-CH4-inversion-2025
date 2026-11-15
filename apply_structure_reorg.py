# -*- coding: utf-8 -*-
"""One-time reorganization of GCP/GMB2025 Python modules.

From c_gcpv3:
    python apply_structure_reorg.py           # preview
    python apply_structure_reorg.py --apply   # migrate, backup, verify
    python apply_structure_reorg.py --check   # verify an existing migration

Never edits numerical formulas, Fortran code, input datasets, or results.
"""

import argparse
import ast
import hashlib
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path


# ============================================================
# === Names, as agreed for the GCP/GMB 2025 project

ROOT = Path(__file__).resolve().parent

MOVES = {
    'b1_inca_OH.py': 'd0_inca_oh.py',
    'c1_conflux_nc.py': 'd1_conflux_nc.py',
    'c2_conflux_gt3.py': 'd2_conflux_gt3.py',
    'c3_obs_files.py': 'd3_obs_files.py',
    'c4_diag_ch4.py': 'd4_diag_ch4.py',
    's0_server.py': 'p0_postprocess.py',
    's1_Burden.py': 'p1_burden.py',
    's2_LossCorr.py': 'p2_loss_corr.py',
    's3_write_grd.py': 'p3_write_grd.py',
    's4_write_gt3.py': 'p4_write_gt3.py',
    's5_LossCorrPlt.py': 'p5_loss_corr_plot.py',
}

# --- _main.py in some older builds imported i1_inca_OH despite the actual
#     file being b1_inca_OH.py. This is an old, incorrect import, not the
#     current inversion module i1_inversion_runner.py.
ALIASES = {'i1_inca_OH': 'd0_inca_oh'}
OLD_MODULES = {a[:-3]: b[:-3] for a, b in MOVES.items()}
OLD_MODULES.update(ALIASES)


# ============================================================
# === Source handling

def project_files():
    scripts = list(ROOT.glob('*.py'))
    tests = list((ROOT / 'test').glob('*.py'))
    return sorted(set(scripts + tests))


def replace_refs(data):
    # --- Replace only complete Python module names, never strings such as
    #     s1_ch4_burloss.f90 or arbitrary 's1'/'c1' tokens.
    for old, new in sorted(OLD_MODULES.items(), key=lambda x: -len(x[0])):
        data = re.sub(r'(?<![A-Za-z0-9_])' + re.escape(old) + r'(?![A-Za-z0-9_])',
                      new, data)
    return data


def transform_main(source):
    # --- Add '--steps post' as an alias of the existing '--steps server'.
    #     Keeping 'server' fully backward-compatible avoids breaking scripts.
    marker = "    switches = dict(RUN)\n    if selected is not None:\n"
    if marker in source and "'server' if step == 'post' else step" not in source:
        source = source.replace(marker,
                                "    switches = dict(RUN)\n"
                                "    if selected is not None:\n"
                                "        selected = ['server' if step == 'post' else step for step in selected]\n")
    choice_old = "choices=list(RUN) + ['all']"
    if choice_old in source:
        source = source.replace(choice_old, "choices=list(RUN) + ['all', 'post']")
    source = source.replace("Start Server *-*-*", "Start Post-processing *-*-*")
    return source


def transform(filename, original):
    source = original.decode('utf-8')
    source = replace_refs(source)
    if filename == '_main.py':
        source = transform_main(source)
    ast.parse(source, filename=filename)
    return source.encode('utf-8')


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def all_refs_clean():
    problems = []
    for path in project_files():
        if path.name == 'apply_structure_reorg.py':
            continue
        data = path.read_text(encoding='utf-8')
        try:
            ast.parse(data, filename=str(path))
        except SyntaxError as exc:
            problems.append(str(path) + ': invalid syntax: ' + str(exc))
            continue
        for old in OLD_MODULES:
            if re.search(r'(?<![A-Za-z0-9_])' + re.escape(old) + r'(?![A-Za-z0-9_])', data):
                problems.append(str(path) + ': obsolete module ' + old)
    return problems


def verify():
    missing = []
    for old, new in MOVES.items():
        if not (ROOT / new).is_file():
            missing.append(new)
        if (ROOT / old).exists():
            missing.append('leftover legacy file ' + old)
    errors = missing + all_refs_clean()
    main_path = ROOT / '_main.py'
    if main_path.is_file():
        main = main_path.read_text(encoding='utf-8')
        if 'p0_postprocess' not in main:
            errors.append('_main.py does not refer to p0_postprocess')
    post_path = ROOT / 'p0_postprocess.py'
    if post_path.is_file():
        post = post_path.read_text(encoding='utf-8')
        for name in ['p1_burden', 'p2_loss_corr', 'p3_write_grd',
                     'p4_write_gt3', 'p5_loss_corr_plot']:
            if name not in post:
                errors.append('p0_postprocess.py missing expected import: ' + name)
    if errors:
        raise RuntimeError('STRUCTURE CHECK FAILED:\n  ' + '\n  '.join(errors))
    print('STRUCTURE CHECK: PASS')
    print('Python files: d0-d4, i1-i3 (unchanged), p0-p5, t*, v*')
    print('Fortran and inversion/: unchanged')


# ============================================================
# === Migration (sources are the live SERVER files, not old copies)

def migrate():
    old = {name: ROOT / name for name in MOVES}
    new = {name: ROOT / target for name, target in MOVES.items()}

    if all(not path.exists() for path in old.values()) and all(p.is_file() for p in new.values()):
        verify()
        print('Already migrated, nothing to change.')
        return

    missing_required = [name for name in MOVES if not old[name].is_file()]
    if missing_required:
        raise FileNotFoundError('Missing files (refusing partial migration): ' + ', '.join(missing_required))

    conflicts = [new[name].name for name in MOVES if new[name].exists()]
    if conflicts:
        raise FileExistsError('New filename already exists (refusing overwrite): ' + ', '.join(conflicts))

    # --- Read all intended changes and validate before creating backup.
    originals = {}
    staged = {}
    for name, old_path in old.items():
        target = new[name]
        original = old_path.read_bytes()
        originals[old_path] = original
        staged[target] = transform(target.name, original)

    for path in project_files():
        if path.name == 'apply_structure_reorg.py' or path in old.values():
            continue
        source = path.read_bytes()
        modified = transform(path.name, source)
        if modified != source:
            originals[path] = source
            staged[path] = modified

    readme = ROOT / 'README.md'
    if readme.is_file():
        source = readme.read_bytes()
        modified = replace_refs(source.decode('utf-8')).encode('utf-8')
        if modified != source:
            originals[readme] = source
            staged[readme] = modified

    now = datetime.now().strftime('%Y%m%d_%H%M%S')
    backup = ROOT.parent / ('backup_c_gcpv3_structure_' + now)
    if backup.exists():
        raise FileExistsError('Backup directory exists: ' + str(backup))
    backup.mkdir(parents=True)

    for path, data in originals.items():
        dest = backup / path.relative_to(ROOT)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(data)

    report = []
    report.append('GCP/GMB2025 Python module reorganization')
    report.append('Backup: ' + str(backup))
    report.append('')
    report.append('Old -> New:')
    for old_name, new_name in MOVES.items():
        report.append(old_name + ' -> ' + new_name)
    report.append('')
    report.append('Checksums of source files:')
    for path, data in sorted(originals.items()):
        report.append(sha256(data) + '  ' + str(path.relative_to(ROOT)))
    (backup / 'MIGRATION_REPORT.txt').write_text('\n'.join(report) + '\n', encoding='utf-8')

    created = []
    changed_existing = []
    try:
        for path, data in staged.items():
            path.write_bytes(data)
            if path in new.values():
                created.append(path)
            else:
                changed_existing.append(path)

        # --- Ensure migrated Python is syntactically valid without executing
        #     any scientific analysis or launching ACTM/TDI.
        for path in staged:
            if path.suffix == '.py':
                ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
        for old_name in MOVES:
            old_path = old[old_name]
            if old_path.is_file():
                old_path.unlink()
        verify()

    except Exception:
        for path, data in originals.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        for path in created:
            if path.exists():
                path.unlink()
        print('Migration FAILED. All old files restored. Backup retained: ' + str(backup))
        raise

    print('MIGRATION PASS')
    print('Backup:', backup)
    print('Old Python module names removed. No Fortran or data files were modified.')


def main():
    parser = argparse.ArgumentParser(description='Rename input and post-processing Python modules')
    parser.add_argument('--apply', action='store_true', help='Execute migration (with backup)')
    parser.add_argument('--check', action='store_true', help='Check migrated structure')
    args = parser.parse_args()
    if args.check:
        verify()
    elif args.apply:
        migrate()
    else:
        print('Project:', ROOT)
        print('Planned changes:')
        for old_name, new_name in MOVES.items():
            print(' ', old_name, '->', new_name)
        print('Use --apply to migrate, --check to verify. Fortran untouched.')


if __name__ == '__main__':
    main()
