# -*- coding: utf-8 -*-
"""Run configuration and binary smoke/regression tests; requires no real ACTM files."""
import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parent
project_dir = root.parent
scripts = ['test_config.py', 'test_f2.py', 'test_f3.py', 'test_filter.py',
           'test_f4_f7_config.py', 'test_user_merge.py', 'test_f4_binary.py', 'test_f6_regression.py', 'test_inversion.py',
           'test_c4_diag_ch4.py', 'test_inv_duplicate_names.py', 'test_read_globalview_guard.py',
           'test_valid_months.py']
for name in scripts:
    print('\n===', name, '===', flush=True)
    subprocess.run([sys.executable, str(root/name)], cwd=project_dir, check=True)
print('\nAll configured tests passed.')
