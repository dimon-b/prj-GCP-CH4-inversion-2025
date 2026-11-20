# -*- coding: utf-8 -*-
"""D1 (1998-2026) offline regression tests. Run from project root:
    python test/test_d1_fluxes.py
Does not read production inputs or modify project outputs.
Optional real-file validation: python test/test_d1_fluxes.py --production
"""
import argparse
import importlib.util
import os
import struct
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
SCRIPT = ROOT / 'd1_prior_fluxes.py'
if not SCRIPT.is_file():
    raise FileNotFoundError('Missing D1 module at ' + str(SCRIPT))

# The production SetCase class needs project-specific files. Offline unit tests
# instantiate D1 without calling __init__, so an optional stub is sufficient.
if '_set_case' not in sys.modules:
    try:
        import _set_case
    except ImportError:
        sys.modules['_set_case'] = types.SimpleNamespace(SetCase=object)
spec = importlib.util.spec_from_file_location('d1_prior_fluxes', SCRIPT)
d1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d1)

START, END = 1998, 2026
MONTHS = (END - START + 1) * 12


def runner():
    obj = object.__new__(d1.PriorFluxes)
    obj.yr_s, obj.yr_e, obj.years_nc = START, END, [2000, 2024]
    obj.R = 6371000.0
    return obj


def source_data():
    times = pd.date_range('2000-01-01', '2024-12-01', freq='MS')
    t = np.arange(300, dtype='float64')[:, None, None]
    wet = np.broadcast_to(12.0 + t / 4.0, (300, 2, 3)).copy()
    soil = np.broadcast_to(-1.0 - t / 300.0, (300, 2, 3)).copy()
    soil[:, 0, 0] = np.nan
    fossil = np.broadcast_to(3.0 + 0*t, (300, 2, 3)).copy()
    ds = xr.Dataset({'flux_ch4_wetlands': (('time', 'lat', 'lon'), wet),
                     'flux_ch4_soils': (('time', 'lat', 'lon'), soil),
                     'flux_ch4_fossil': (('time', 'lat', 'lon'), fossil)},
                    coords={'time': times, 'lat': [-0.5, 0.5], 'lon': [0.5, 1.5, 2.5]})
    ds['flux_ch4_soils'] = ds['flux_ch4_soils'].fillna(0.)
    sectors = ['flux_ch4_wetlands', 'flux_ch4_soils', 'flux_ch4_fossil']
    ds['flux_ch4_prior'] = sum(ds[v] for v in sectors)
    ds['flux_ch4_prior_soil0'] = ds['flux_ch4_prior'] - ds['flux_ch4_soils']
    return ds, sectors


def make_gt3(path, n=MONTHS, bad_record=None, bad_value=False):
    # 64 16-byte strings + 128x64 real32, Fortran big-endian sequential records.
    field = np.zeros((64, 128), dtype='>f4')
    with Path(path).open('wb') as f:
        for k in range(n):
            year, month = START+k//12, k % 12 + 1
            head = [''] * 64
            head[2], head[15], head[37] = 'QFLUX07', 'kg-CH4/m2/s', 'UR4'
            head[24], head[26], head[63] = str(k), f'{year:04d}{month:02d}16 120000', str(128*64)
            if k == bad_record:
                head[26] = '19000101 000000'
            h = ''.join(v.ljust(16)[:16] for v in head).encode('ascii')
            if k == bad_record and bad_value:
                field[0, 0] = np.nan
                head[26] = f'{year:04d}{month:02d}16 120000'
                h = ''.join(v.ljust(16)[:16] for v in head).encode('ascii')
            payload = field.tobytes()
            f.write(struct.pack('>i', len(h))); f.write(h); f.write(struct.pack('>i', len(h)))
            f.write(struct.pack('>i', len(payload))); f.write(payload); f.write(struct.pack('>i', len(payload)))
            if bad_value:
                field[0, 0] = 0.0


def write_regional_table(path, bad=False):
    with Path(path).open('w') as f:
        f.write('reg# ' + ' '.join(str(i) for i in range(1, 13)) + '\n')
        for i in range(1, 55):
            values = ['nan' if bad and i == 20 and m == 5 else '1.25' for m in range(12)]
            f.write(str(i) + ' ' + ' '.join(values) + '\n')


def write_uncertainty_nc(path, months=MONTHS, negative=False):
    values = np.zeros((2, months, 54), dtype=np.float32)
    values[1, :, :] = -1 if negative else 0.5
    time = np.array([year+(month-0.5)/12 for year in range(START, END+1)
                     for month in range(1, 13)], dtype=np.float32)[:months]
    ds = xr.Dataset({'sources_land': (('datatype_count', 'month', 'land_region'), values),
                     'time': ('month', time)},
                    coords={'year': np.arange(29), 'land_region': np.arange(54),
                            'datatype_count': [0, 1]})
    ds.to_netcdf(path)


class D1Tests(unittest.TestCase):
    def setUp(self):
        self.obj = runner()
        self.ds, self.sectors = source_data()
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)

    def test_01_sector_balance_negative_soil(self):
        self.obj.check_sector_balance(self.ds, self.sectors)
        self.assertLess(float(self.ds['flux_ch4_soils'].min()), 0)
        np.testing.assert_allclose(self.ds['flux_ch4_prior_soil0'],
                                   self.ds['flux_ch4_wetlands']+self.ds['flux_ch4_fossil'])

    def test_02_sector_mismatch_rejected(self):
        bad = self.ds.copy(deep=True)
        bad['flux_ch4_prior'].values[4, 0, 0] += 100
        with self.assertRaisesRegex(ValueError, 'sector sum'):
            self.obj.check_sector_balance(bad, self.sectors)

    def test_03_soil0_mismatch_rejected(self):
        bad = self.ds.copy(deep=True)
        bad['flux_ch4_prior_soil0'].values[1, 0, 1] = 9999
        with self.assertRaisesRegex(ValueError, 'soil0'):
            self.obj.check_sector_balance(bad, self.sectors)

    def test_04_nonfinite_sector_rejected(self):
        bad = self.ds.copy(deep=True)
        bad['flux_ch4_wetlands'].values[5, 0, 0] = np.nan
        with self.assertRaises(ValueError):
            self.obj.check_sector_balance(bad, self.sectors)

    def test_05_extension_all_fields(self):
        ext = self.obj.extend_months(self.ds)
        self.assertEqual(ext.sizes['time'], MONTHS)
        expected = pd.date_range('1998-01-01', '2026-12-01', freq='MS')
        np.testing.assert_array_equal(ext.time.values, expected.values)
        for dest, orig in [(1998, 2000), (1999, 2000), (2025, 2024), (2026, 2024)]:
            for v in self.ds.data_vars:
                np.testing.assert_array_equal(ext[v].sel(time=str(dest)).values,
                                              ext[v].sel(time=str(orig)).values)

    def test_06_missing_source_month_rejected(self):
        with self.assertRaisesRegex(ValueError, 'original monthly'):
            self.obj.extend_months(self.ds.isel(time=slice(1, None)))

    def test_07_missing_output_month_rejected(self):
        ext = self.obj.extend_months(self.ds).isel(time=slice(0, -1))
        with self.assertRaisesRegex(ValueError, '348 months'):
            self.obj.check_extended_ds(ext)

    def test_08_corrupt_boundary_copy_rejected(self):
        ext = self.obj.extend_months(self.ds)
        ext['flux_ch4_fossil'].values[0, 0, 0] += 10
        with self.assertRaisesRegex(ValueError, 'Edge copy differs'):
            self.obj.check_extended_ds(ext)

    def test_09_duplicate_output_month_rejected(self):
        ext = self.obj.extend_months(self.ds)
        dates = ext.time.values.copy(); dates[13] = dates[12]
        ext = ext.assign_coords(time=dates)
        with self.assertRaisesRegex(ValueError, 'time axis'):
            self.obj.check_extended_ds(ext)

    def test_10_gt3_348_records(self):
        p = self.folder / 'good.gt3'
        make_gt3(p)
        self.obj.check_gt3(p)

    def test_11_gt3_324_records_rejected(self):
        p = self.folder / 'old.gt3'
        make_gt3(p, n=324)
        with self.assertRaisesRegex(ValueError, '348'):
            self.obj.check_gt3(p)

    def test_12_gt3_bad_date_rejected(self):
        p = self.folder / 'date.gt3'
        make_gt3(p, bad_record=17)
        with self.assertRaisesRegex(ValueError, 'time header'):
            self.obj.check_gt3(p)

    def test_13_gt3_nan_rejected(self):
        p = self.folder / 'nan.gt3'
        make_gt3(p, bad_record=17, bad_value=True)
        with self.assertRaisesRegex(ValueError, 'flux values'):
            self.obj.check_gt3(p)

    def test_14_regional_valid(self):
        p = self.folder / 'f4.txt'; write_regional_table(p)
        self.obj.check_regional_table(p)

    def test_15_regional_nan_rejected(self):
        p = self.folder / 'f4_bad.txt'; write_regional_table(p, bad=True)
        with self.assertRaisesRegex(ValueError, 'non-finite'):
            self.obj.check_regional_table(p)

    def test_16_regional_missing_region_rejected(self):
        p = self.folder / 'f4_short.txt'; write_regional_table(p)
        p.write_text('\n'.join(p.read_text().splitlines()[:-1])+'\n')
        with self.assertRaisesRegex(ValueError, '54 region rows'):
            self.obj.check_regional_table(p)

    def test_17_uncertainty_valid(self):
        p = self.folder / 'p50.nc'; write_uncertainty_nc(p)
        self.obj.check_prior_uncertainty_nc(p)

    def test_18_uncertainty_old_324_rejected(self):
        p = self.folder / 'old_p50.nc'; write_uncertainty_nc(p, months=324)
        with self.assertRaisesRegex(ValueError, 'month/year'):
            self.obj.check_prior_uncertainty_nc(p)

    def test_19_uncertainty_negative_uncertainty_rejected(self):
        p = self.folder / 'bad_p50.nc'; write_uncertainty_nc(p, negative=True)
        with self.assertRaisesRegex(ValueError, 'negative prior uncertainty'):
            self.obj.check_prior_uncertainty_nc(p)

    def test_20_old_prior_nc_rejected_without_overwrite(self):
        # This regression reproduces the user's traceback: get_inp_flx sees
        # both old files and must reject them without changing a single byte.
        old = xr.Dataset({'flux_ch4_prior': (('time', 'lat', 'lon'),
                          np.zeros((324, 1, 1), dtype=np.float32))},
                         coords={'time': pd.date_range('1999-01-01', periods=324, freq='MS'),
                                 'lat': [0.], 'lon': [0.]})
        p = self.folder / 'gcp25_inv1.nc'
        q = self.folder / 'gcp25_inv1_full.nc'
        old.to_netcdf(p); old.to_netcdf(q)
        p0, q0 = p.read_bytes(), q.read_bytes()
        self.obj.flx_inp_dir, self.obj.flx_inp_nc = str(self.folder), 'missing.nc'
        self.obj.flx_apr_dir = str(self.folder)
        self.obj.flx_apr_nc, self.obj.flx_apr_nc_full = p.name, q.name
        with self.assertRaisesRegex(ValueError, '324.*348'):
            self.obj.get_inp_flx()
        self.assertEqual(p.read_bytes(), p0)
        self.assertEqual(q.read_bytes(), q0)

    def test_21_nc3_fortran_compatible_writer(self):
        p = self.folder / 'gcp25_inv1.nc'
        q = self.folder / 'gcp25_inv1_full.nc'
        self.obj.write_nc(self.ds, p, q)
        for name in (p, q):
            self.assertTrue(self.obj.is_fortran_nc(name))
            with name.open('rb') as f:
                self.assertEqual(f.read(4), b'CDF\x02')
            from scipy.io import netcdf_file
            with netcdf_file(name, 'r', mmap=False) as nc:
                self.assertTrue(nc.variables['time'].units.startswith(b'days since'))
                self.assertEqual(nc.variables['flux_ch4_prior'].units, b'kg-CH4/m2/s')
                self.assertEqual(nc.variables['time'].shape[0], 348)

    def test_22_rebuild_legacy_netcdf4_files(self):
        p = self.folder / 'gcp25_inv1.nc'
        q = self.folder / 'gcp25_inv1_full.nc'
        p.write_bytes(b'\x89HDFlegacy-prior')
        q.write_bytes(b'\x89HDFlegacy-full')
        inp = self.folder / 'GCP_Prior_CH4_fluxes.nc'
        self.ds.to_netcdf(inp)
        self.obj.flx_inp_dir = str(self.folder)
        self.obj.flx_inp_nc = inp.name
        self.obj.flx_apr_dir = str(self.folder)
        self.obj.flx_apr_nc = p.name
        self.obj.flx_apr_nc_full = q.name
        from unittest.mock import patch
        with patch.object(self.obj, 'check_prior_nc'), \
             patch.object(self.obj, 'manage_ds', return_value=self.ds), \
             patch.object(self.obj, 'check_total'):
            self.obj.get_inp_flx()
        self.assertEqual(p.read_bytes()[:4], b'CDF\x02')
        self.assertEqual(q.read_bytes()[:4], b'CDF\x02')
        self.assertEqual(p.with_name(p.name+'.netcdf4_backup').read_bytes(), b'\x89HDFlegacy-prior')
        self.assertEqual(q.with_name(q.name+'.netcdf4_backup').read_bytes(), b'\x89HDFlegacy-full')

    def test_23_fortran_periods_and_configs(self):
        # Cheap static preflight, not a substitute for compiler/real-data run.
        for filename, stage in [('df1_regrid_ch4_t42.f90', 'd1_regrid_ch4_t42'),
                                ('df1_aggregate_ch4_regions.f90', 'd1_aggregate_ch4_regions'),
                                ('df1_write_prior_uncertainty.f90', 'd1_write_prior_uncertainty')]:
            p = ROOT / 'fortran' / filename
            self.assertTrue(p.is_file(), 'Missing Fortran source: ' + filename)
            text = p.read_text().lower().replace(' ', '')
            self.assertIn('s_year=1998,f_year=2026', text, filename)
            self.assertIn('namelist/' + stage + '/', text, filename)


    def test_24_stale_regional_is_backed_up_and_recomputed(self):
        from unittest.mock import patch
        inp = self.folder / 'flux.gt3'
        mask = self.folder / 'mask.grd'
        out = self.folder / 'regional.dat'
        inp.write_bytes(b'new GT3 input')
        mask.write_bytes(b'mask')
        write_regional_table(out)
        old = out.read_bytes()
        self.obj.d1_gt3_file = str(inp)
        self.obj.d1_region_mask_file = str(mask)
        self.obj.d1_region_table_file = str(out)
        self.obj.write_case_namelist = lambda *a: 'test.nml'
        def run(stage, config):
            self.assertEqual(stage, 'regions')
            write_regional_table(out)
        with patch.object(self.obj, 'run_exe', side_effect=run):
            self.obj.get_regional_flux()
        backups = list(self.folder.glob('regional.dat.backup_*'))
        self.assertEqual(len(backups), 1)
        self.assertEqual(backups[0].read_bytes(), old)
        self.assertTrue(Path(str(out) + '.d1_1998-2026.sha256').is_file())
        with patch.object(self.obj, 'run_exe', side_effect=AssertionError('F4 reran')):
            self.obj.get_regional_flux()

    def test_25_stale_uncertainty_is_backed_up_and_recomputed(self):
        from unittest.mock import patch
        import time
        regional = self.folder / 'regional.dat'
        regional.write_text('new F4 output')
        prefix = self.folder / 'sourcepriorch4_p'
        output = [Path(str(prefix) + str(pc) + '.nc') for pc in (30, 50, 99)]
        for p in output:
            write_uncertainty_nc(p)
            os.utime(p, (1000000, 1000000))
        self.obj.d1_uncertainty_input_file = str(regional)
        self.obj.d1_uncertainty_output_prefix = str(prefix)
        self.obj.write_case_namelist = lambda *a: 'test.nml'
        def run(stage, config):
            self.assertEqual(stage, 'uncertainty')
            for p in output:
                write_uncertainty_nc(p)
        with patch.object(self.obj, 'run_exe', side_effect=run):
            self.obj.write_prior_uncertainty()
        for p in output:
            self.assertTrue(p.exists())
            self.assertEqual(len(list(self.folder.glob(p.name + '.backup_*'))), 1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--production', action='store_true')
    args = parser.parse_args()
    if args.production:
        import _set_case
        if not hasattr(_set_case, 'SetCase') or _set_case.SetCase is object:
            raise RuntimeError('Production mode requires real _set_case.py')
        case = os.environ.get('GCP_CASE', os.environ.get('CH4_CASE', 'INCA')).upper()
        os.environ['GCP_CASE'] = case
        os.environ['CH4_CASE'] = case
        obj = d1.PriorFluxes(case=case)
        prior = Path(obj.flx_apr_dir) / obj.flx_apr_nc
        full = Path(obj.flx_apr_dir) / obj.flx_apr_nc_full
        obj.check_prior_nc(prior)
        obj.check_prior_nc(full)
        inp = Path(obj.d1_regrid_input_dir) / obj.d1_regrid_input_nc
        gt3 = Path(obj.d1_gt3_file)
        obj.check_gt3(gt3)
        obj.check_flux_conservation(inp, gt3)
        obj.check_regional_table(Path(obj.d1_region_table_file))
        for pc in (30, 50, 99):
            obj.check_prior_uncertainty_nc(Path(str(obj.d1_uncertainty_output_prefix) + str(pc) + '.nc'))
        print('PRODUCTION D1 CHECKS PASSED')
    else:
        unittest.main(argv=[sys.argv[0]], verbosity=2)
