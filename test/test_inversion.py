# -*- coding: utf-8 -*-
"""Test inversion control wiring and NetCDF input schema without actual Bayesian inversion."""

# --- Import the project modules from c_gcpv3/ (one level above test/)
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import os
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.io import netcdf_file

from inversion.check_inv_inputs import inspect_observations, inspect_sources
from inversion.inversion_runner import get_control_text, get_inversion_paths, run_inversion, sha256_file


# ============================================================
# === Synthetic F7 monthly observations

def make_obs(path):
    with netcdf_file(str(path), 'w') as nc:
        nc.createDimension('station_count', 245)
        nc.createDimension('tlen', 324)
        nc.createDimension('data_type', 2)
        nc.createDimension('namelen', 45)
        names = nc.createVariable('site_name', 'c', ('station_count','namelen'))
        chars = np.full((245,45), b' ', dtype='S1')
        for i in range(245):
            name = ('site_%03d' % i).ljust(45).encode('ascii')
            chars[i] = np.frombuffer(name, dtype='S1')
        names[:] = chars
        ids = nc.createVariable('station_count', 'i', ('station_count',))
        ids[:] = np.arange(1, 246)
        time = nc.createVariable('tlen', 'f', ('tlen',))
        time[:] = 1999. + (np.arange(324) + .5) / 12.
        data = nc.createVariable('tracer_data', 'f', ('tlen','station_count','data_type'))
        values = np.empty((324,245,2), dtype='f4')
        values[:,:,0] = 12.5
        values[:,:,1] = 7.0
        data[:] = values


# ============================================================
# === Synthetic F5 regional prior

def make_source(path):
    with netcdf_file(str(path), 'w') as nc:
        for name, count in [('land_region',54),('month',324),('year',27),('presub_count',4),
                            ('presubp1',5),('datatype_count',2),('namelen',10)]:
            nc.createDimension(name,count)
        nc.createVariable('sources_land','f',('datatype_count','month','land_region'))[:] = np.ones((2,324,54),dtype='f4')
        nc.createVariable('sources_presub','f',('datatype_count','year','presub_count'))[:] = np.ones((2,27,4),dtype='f4')
        nc.createVariable('sources_offset','f',('datatype_count',))[:] = [0.,0.]
        nc.createVariable('time','f',('month',))[:] = 1999. + (np.arange(324) + .5) / 12.
        nc.createVariable('presub_name','c',('presubp1','namelen'))[:] = np.full((5,10), b' ',dtype='S1')


def run_tests():
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        inv = root / 'inv_dir'
        inv.mkdir()
        site_file = inv / 'obsrvCH4_60sites.txt'
        site_file.write_text('header\n245 no name\n1 9999 site_000 0 0 0\n2 9999 site_001 0 0 0\n', encoding='utf-8')

        template = root / 'control'
        old_control = Path(__file__).resolve().parents[1] / 'inversion' / 'control.example'
        template.write_text(old_control.read_text(encoding='utf-8'), encoding='utf-8')

        # --- All master inputs and Green's functions are available in this fixture.
        masters = root / 'files4run'
        masters.mkdir()
        for name in ('spatial-srccor', 'global-output', 'spatial-groups-tdi',
                     'spatial-consts-none', 'new-recipes'):
            (masters / name).write_text('frozen test master: ' + name + '\n')
        greens = root / 'gmats'
        greens.mkdir()
        (greens / 'g00.dat').write_bytes(b'green-matrix-test')
        old = template.read_text(encoding='utf-8')
        old = old.replace('/S/data01/G5070/y0715/prj_GCP_v25/inv_dir/gmats/', str(greens))
        # --- Simulate an older template; runner must force dumpmats off.
        template.write_text(old.replace('dumpmats = .false.', 'dumpmats = .true.'), encoding='utf-8')

        script = root / 'tdi_ch4'
        script.write_text('#!/bin/sh\n'
                          'if [ "$1" = "--valid-months-version" ]; then echo TDI_VALID_MONTHS_V1; exit 0; fi\n'
                          '[ -f obsrvCH4_9925.nc ] && [ ! -L obsrvCH4_9925.nc ] || exit 3\n'
                          '[ -f sourcepriortdi.nc ] && [ ! -L sourcepriortdi.nc ] || exit 4\n'
                          '[ -f obsrvCH4_sites.txt ] && [ ! -L obsrvCH4_sites.txt ] || exit 5\n'
                          '[ -d files4run ] && [ ! -L files4run ] || exit 6\n'
                          '[ -f control ] && [ -f tdi_ch4 ] || exit 7\n'
                          'grep -q "datflist = '"'"'obsrvCH4_sites.txt'"'"'" control || exit 8\n'
                          'grep -qi "dumpmats[[:space:]]*=[[:space:]]*\\.false\\." control || exit 9\n'
                          'printf "mock prior" > prior.sources\n'
                          'printf "mock predicted" > predicted.sources\n', encoding='utf-8')
        script.chmod(0o755)

        for case in ('CYC','INCA'):
            prefix = inv / ('obsrvCH4_' + case.lower() + '_9925')
            obs = Path(str(prefix) + 'ctlinv1.nc')
            source_prefix = inv / 'sourceprior' / 'sourcepriorch4_p'
            source_prefix.parent.mkdir(exist_ok=True)
            source = Path(str(source_prefix) + '30.nc')
            make_obs(obs)
            if not source.is_file():
                make_source(source)
            info = inspect_observations(obs,site_file,strict=True)
            assert info['n_active_stations'] == 2
            assert not info['n_invalid_active']
            assert inspect_sources(source)['n_regions'] == 54

            setup = SimpleNamespace(hcase=case,inv_control_template=str(template),exe_inv=str(script),sites_f=str(site_file),
                                    yr_s=1999,yr_e=2025, invcases=['inv1'],unpcases=['p30','p50','p99'],
                                    unxcases=['ctl','ux2','ux4'], f7_output_nc_prefix=str(prefix),
                                    f5_output_prefix=str(source_prefix), inv_wrk_dir=str(root/'results2025'/case),icase='s060')
            control = get_control_text(template,case,site_file,1999,2025)
            assert "datf     = 'obsrvCH4_9925.nc'" in control
            assert "sourcef  = 'sourcepriortdi.nc'" in control
            run_inversion(case, setup, dry_run=True)
            run_inversion(case, setup, check_nc=True)
            _, _, work_dir = get_inversion_paths(setup,case)
            assert (work_dir/'prior.sources').read_text() == 'mock prior'
            assert (work_dir/'predicted.sources').read_text() == 'mock predicted'
            assert (work_dir/'obsrvCH4_9925.nc').is_file()
            assert not (work_dir/'obsrvCH4_9925.nc').is_symlink()
            assert not (work_dir/'sourcepriortdi.nc').is_symlink()
            assert not (work_dir/'obsrvCH4_sites.txt').is_symlink()
            assert not (work_dir/'files4run').is_symlink()
            assert not (work_dir/'tdi_ch4').is_symlink()
            assert (work_dir/'control').read_text().count("'obsrvCH4_sites.txt'") == 1
            assert 'dumpmats = .false.' in (work_dir/'control').read_text()
            assert not (work_dir/'vardump').exists()
            assert 'vardump_enabled=false' in (work_dir/'RUN_INPUTS.txt').read_text()
            assert (work_dir/'RUN_INPUTS.txt').is_file()
            assert (work_dir/'RUN_INPUTS.sha256').is_file()
            for checksum_line in (work_dir/'RUN_INPUTS.sha256').read_text().splitlines():
                checksum, name = checksum_line.split('  ', 1)
                assert sha256_file(work_dir/name) == checksum
            assert 'GF_mode=SHARED_NOT_COPIED' in (work_dir/'RUN_INPUTS.txt').read_text()
            assert 'source_F7=' in (work_dir/'RUN_INPUTS.txt').read_text()
            assert (work_dir/'tdi_ch4').stat().st_mode & 0o111
            try:
                run_inversion(case, setup, check_nc=True)
            except FileExistsError:
                pass
            else:
                raise AssertionError('Existing result directory was overwritten')

            # --- Missing active observation months must stop production.
            with netcdf_file(str(obs), 'a') as nc:
                v = nc.variables['tracer_data']
                v[0,0,0] = -999.
            assert inspect_observations(work_dir/'obsrvCH4_9925.nc', site_file=None, strict=True)['n_invalid_all'] == 0
            try:
                inspect_observations(obs,site_file,strict=True)
            except ValueError as e:
                assert 'missing' in str(e)
            else:
                raise AssertionError('Invalid active station/month accepted')

        # --- Copy must survive later changes to the original files.
        frozen_dir = get_inversion_paths(setup, 'INCA')[2]
        orig_snapshot = (frozen_dir/'obsrvCH4_sites.txt').read_bytes()
        orig_master = (frozen_dir/'files4run/spatial-srccor').read_bytes()
        orig_src_hash = sha256_file(frozen_dir/'sourcepriortdi.nc')
        site_file.write_text('changed later\n', encoding='utf-8')
        (masters/'spatial-srccor').write_text('changed later\n', encoding='utf-8')
        source.write_bytes(b'changed later')
        assert (frozen_dir/'obsrvCH4_sites.txt').read_bytes() == orig_snapshot
        assert (frozen_dir/'files4run/spatial-srccor').read_bytes() == orig_master
        assert sha256_file(frozen_dir/'sourcepriortdi.nc') == orig_src_hash

        # --- Control must guard wrong year and overlong paths.
        broken = template.read_text().replace('lastdat = 2025','lastdat = 2024')
        template.write_text(broken)
        try:
            get_control_text(template,'CYC',site_file,1999,2025)
        except ValueError as e:
            assert 'lastdat' in str(e)
        else:
            raise AssertionError('Wrong inversion interval accepted')
        print('Inversion wiring and input NetCDF tests: PASS')


if __name__ == '__main__':
    run_tests()
