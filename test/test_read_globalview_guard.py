# -*- coding: utf-8 -*-
"""Check the NEW valid-month reader, not the retired read_globalview.f guard.

The numerical test is test_valid_months.py. No actual NetCDF or TDI run here.
"""
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'inversion' / 'read_globalview.f90'
DRIVER = ROOT / 'inversion' / 'tdi.f'


def main():
    if not SOURCE.is_file():
        raise FileNotFoundError('Missing updated Fortran reader: ' + str(SOURCE))
    if not DRIVER.is_file():
        raise FileNotFoundError('Missing updated TDI driver: ' + str(DRIVER))

    text = SOURCE.read_text(encoding='utf-8').lower()
    driver = DRIVER.read_text(encoding='utf-8').lower()
    for marker in ('module gv_valid_months', 'logical function gv_valid',
                   'subroutine count_globalview', 'subroutine read_globalview',
                   'ieee_is_finite', 'if (count_only) cycle'):
        assert marker in text, 'Missing expected reader code: ' + marker
    for marker in ('call count_globalview(', 'call read_globalview(',
                   'if (nread.ne.ntotdat)', 'tdi_valid_months_v1'):
        assert marker in driver, 'Missing TDI integration: ' + marker

    compiler = shutil.which('gfortran')
    if not compiler:
        print('Valid-month reader source: PASS; syntax SKIP (gfortran unavailable)')
        return
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        (tmpdir / 'netcdf.inc').write_text('integer, parameter :: nf_noerr=0\ninteger :: nf_open,nf_inq_dimid,nf_inq_dimlen,nf_inq_varid\ninteger :: nf_inq_vardimid,nf_get_var_text,nf_get_var_real,nf_close\ncharacter(len=80) :: nf_strerror\n', encoding='utf-8')
        subprocess.run([compiler, '-fsyntax-only', '-ffree-form',
                        '-ffree-line-length-none', '-J', str(tmpdir),
                        '-I', str(tmpdir), str(SOURCE)], check=True)
    print('Valid-month reader source and syntax: PASS')


if __name__ == '__main__':
    main()
