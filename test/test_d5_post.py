# -*- coding: utf-8 -*-
"""D5 integration tests: actual D1 NetCDF input -> posterior NC -> compiled Fortran GT3.

Run: python test/test_d5_post.py
"""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from d5_post import PostFlux


# ============================================================
# === Build real converter and synthesize an uncompressed D1 NetCDF
class TestD5Post(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):
            raise unittest.SkipTest('gfortran missing')
        cls.data = tempfile.TemporaryDirectory()
        cls.dir = Path(cls.data.name)
        cls.exe = cls.dir / 'df5_post_flux'
        subprocess.run(['gfortran', '-O2', '-ffree-line-length-none', '-fcheck=all',
                        '-o', str(cls.exe), str(ROOT / 'fortran' / 'df5_post_flux.f90')], check=True)
        cls.prior_dir = cls.dir / 'prior'
        cls.prior_dir.mkdir()
        cls.prior = cls.prior_dir / 'gcp25_inv1.nc'
        # --- Coordinates deliberately reversed/cyclic relative to 1-degree D5 grid
        with netcdf_file(str(cls.prior), mode='w', version=2) as nc:
            nc.createDimension('time', 348)
            nc.createDimension('lat', 180)
            nc.createDimension('lon', 360)
            tim = nc.createVariable('time', 'i', ('time',))
            tim.units = b'days since 1998-01-01 00:00:00'
            tim.calendar = b'proleptic_gregorian'
            dates = pd.date_range('1998-01-01', periods=348, freq='MS')
            tim[:] = [(x - dates[0]).days for x in dates]
            lat = nc.createVariable('lat', 'f', ('lat',))
            lon = nc.createVariable('lon', 'f', ('lon',))
            lat[:] = np.arange(89.5, -90, -1, dtype=np.float32)
            lon[:] = np.arange(0.5, 360, 1, dtype=np.float32)
            flux = nc.createVariable('flux_ch4_prior', 'f', ('time', 'lat', 'lon'))
            flux.units = b'kg-CH4/m2/s'
            flux[:] = np.float32(1.0e-13)
            lat_trend = np.arange(180, dtype=np.float32)[:, None] * np.float32(1.0e-15)
            lon_trend = np.arange(360, dtype=np.float32)[None, :] * np.float32(2.0e-16)
            flux[24, :, :] = np.float32(2.0e-13) + lat_trend + lon_trend
            flux[25, :, :] = np.float32(3.0e-13)
            flux[12, :, :] = np.float32(4.0e-13)

    @classmethod
    def tearDownClass(cls):
        cls.data.cleanup()

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.work = Path(self.temp.name)
        self.result_root = self.work / 'results2025' / 'CYC'
        self.inv_dir = self.work / 'inv_dir'
        mask_dir = self.inv_dir / 'mask'
        mask_dir.mkdir(parents=True)
        mask = np.zeros((180, 360), dtype=np.float32)
        mask[90:, :] = 1
        mask.tofile(mask_dir / 'mask80.bin')
        np.full((180, 360), 1e-13, dtype=np.float32).tofile(mask_dir / 'ch4_regn_pulse.grd')
        self.obj = object.__new__(PostFlux)
        self.obj.case = 'CYC'
        self.obj.hcase = 'CYC'
        self.obj.root = self.result_root
        self.obj.diag_root = self.result_root / 'diagnostics' / 'd5'
        self.obj.inv_dir = str(self.inv_dir) + '/'
        self.obj.flx_apr_dir = str(self.prior_dir) + '/'
        self.obj.flx_apr_nc = 'gcp25_inv1.nc'
        self.obj.icase = 's057'
        self.obj.unpcases = ['p30']
        self.obj.unxcases = ['ctl']
        self.obj.invcases = ['inv1']
        self.obj.executable = self.exe
        self.obj.source_code = ROOT / 'fortran' / 'df5_post_flux.f90'
        self.obj.yr_s = 2000
        self.obj.yr_e = 2000
        self.obj.ntime = 2
        self.source = self.obj.get_paths()['source']
        self.source.parent.mkdir(parents=True)
        self.make_sources()

    def tearDown(self):
        self.temp.cleanup()

    def make_sources(self, year=2000, nonzero=2.0):
        with self.source.open('w', encoding='ascii') as out:
            for r in range(1, 55):
                out.write('land_' + format(r, '02d') + '\n')
                out.write('    2 ,          2\n')
                for m in (1, 2):
                    stamp = year + (m - 0.5) / 12
                    coefficient = nonzero if r == 1 else 0
                    out.write(f'{stamp:.6f} {coefficient:.9e} 0.0\n')

    def test_direct_netcdf_prior_and_gt3(self):
        # --- This is a full run of the real Fortran converter and Python D5
        self.obj.run()
        paths = self.obj.get_paths()
        self.assertTrue(paths['nc'].is_file())
        self.assertTrue(paths['gt3'].is_file())
        self.assertFalse(any(self.result_root.rglob('*.grd')))
        self.assertEqual(paths['gt3'].stat().st_size, 2 * 33808)
        with paths['gt3'].open('rb') as gt3:
            self.assertEqual(int.from_bytes(gt3.read(4), 'big'), 1024)
            self.assertIn(b'kg-CH4/m2/s', gt3.read(1024))
            self.assertEqual(int.from_bytes(gt3.read(4), 'big'), 1024)
            self.assertEqual(int.from_bytes(gt3.read(4), 'big'), 128 * 64 * 4)
            field_t42 = np.frombuffer(gt3.read(128 * 64 * 4), dtype='>f4').reshape(64, 128)
            self.assertTrue(np.isfinite(field_t42).all())
            self.assertGreater(float(np.std(field_t42)), 0.0)
        with xr.open_dataset(paths['nc']) as nc:
            self.assertEqual(nc['flux_ch4_post'].shape, (2, 180, 360))
            self.assertEqual(pd.to_datetime(nc.time.values)[0], pd.Timestamp('2000-01-01'))
            post = nc['flux_ch4_post'].isel(time=0).values
            # --- For native D1 lat descending, lon 0-360, output is ascending lat/-180..180
            expected_sw = 2e-13 + 179e-15 + 180 * 2e-16
            expected_se = 2e-13 + 179e-15 + 179 * 2e-16
            np.testing.assert_allclose(post[0, 0], expected_sw, rtol=1.e-6)
            np.testing.assert_allclose(post[0, -1], expected_se, rtol=1.e-6)
            self.assertGreater(post[90, 0], 0)
        budget = pd.read_csv(paths['budget'])
        self.assertEqual(len(budget), 2)
        self.assertEqual(budget.month.tolist(), [1, 2])
        self.assertAlmostEqual(budget.correction_Tg_month.iloc[0], 2 * 31 / 365, places=7)
        self.assertAlmostEqual(budget.correction_Tg_month.iloc[1], 2 * 29 / 365, places=7)
        np.testing.assert_allclose(budget.prior_Tg_month + budget.correction_Tg_month,
                                   budget.post_Tg_month, rtol=3e-7)
        np.testing.assert_allclose(budget.post_Tg_month, budget.gt3_Tg_month, rtol=3e-7)
        self.assertTrue((paths['diag'] / 'post_overview.png').is_file())
        self.assertTrue((paths['diag'] / 'diagnostic_report.txt').is_file())
        self.assertTrue((paths['diag'] / 'annual_budget.csv').is_file())
        self.assertTrue((paths['diag'] / 'regional_2000_2024.csv').is_file())
        self.assertIn('NETCDF -> GT3', (paths['diag'] / 'fortran.log').read_text().upper())

    def test_annual_apr_corr_post_report_is_sum_not_december(self):
        # --- Legacy check_gtotal: APR + CORR = POST in annual mass units.
        import contextlib
        import io
        capture = io.StringIO()
        with contextlib.redirect_stdout(capture):
            self.obj.run()
        out = capture.getvalue()
        self.assertIn('Year         APR         CORR         POST', out)
        self.assertNotIn('Months', out)
        self.assertNotIn('2/12', out)
        self.assertNotIn('D5 processed: 2000 post [Tg/month]', out)
        annual = pd.read_csv(self.obj.get_paths()['diag'] / 'annual_budget.csv')
        monthly = pd.read_csv(self.obj.get_paths()['budget'])
        self.assertEqual(list(annual.year), [2000])
        self.assertNotIn('months', annual.columns)
        self.assertAlmostEqual(annual.APR_Tg_year.iloc[0], monthly.prior_Tg_month.sum(), places=7)
        self.assertAlmostEqual(annual.CORR_Tg_year.iloc[0], monthly.correction_Tg_month.sum(), places=7)
        self.assertAlmostEqual(annual.POST_Tg_year.iloc[0], monthly.post_Tg_month.sum(), places=7)
        self.assertAlmostEqual(annual.APR_Tg_year.iloc[0] + annual.CORR_Tg_year.iloc[0], annual.POST_Tg_year.iloc[0], places=6)
        report = (self.obj.get_paths()['diag'] / 'diagnostic_report.txt').read_text()
        self.assertIn('Year         APR         CORR         POST', report)
        self.assertNotIn('Months', report)

    def test_rerun_replaces_existing_posterior_without_flag(self):
        # --- Re-running D5 must replace NC/GT3 when D4 corrections change.
        self.obj.run()
        paths = self.obj.get_paths()
        with xr.open_dataset(paths['nc']) as ds:
            old_post = ds['flux_ch4_post'].isel(time=0).values.copy()
        old_gt3 = paths['gt3'].read_bytes()
        old_corr = pd.read_csv(paths['budget']).correction_Tg_month.iloc[0]

        self.make_sources(nonzero=4.0)
        self.obj.run()
        with xr.open_dataset(paths['nc']) as ds:
            new_post = ds['flux_ch4_post'].isel(time=0).values.copy()
        new_corr = pd.read_csv(paths['budget']).correction_Tg_month.iloc[0]

        self.assertGreater(float(np.max(np.abs(new_post - old_post))), 0.0)
        self.assertNotEqual(paths['gt3'].read_bytes(), old_gt3)
        self.assertAlmostEqual(new_corr, 2.0 * old_corr, places=7)
        self.assertNotIn('overwrite', str(PostFlux.run.__code__.co_varnames))

    def test_diagnostics_can_reprint_without_recreating_flux(self):
        # --- Re-running reports should not overwrite NetCDF or GT3.
        import contextlib
        import io
        self.obj.run()
        paths = self.obj.get_paths()
        original_gt3 = paths['gt3'].read_bytes()
        original_nc_size = paths['nc'].stat().st_size
        capture = io.StringIO()
        with contextlib.redirect_stdout(capture):
            self.obj.diagnostics(paths)
        self.assertIn('APR + CORR = POST', capture.getvalue())
        self.assertEqual(paths['gt3'].read_bytes(), original_gt3)
        self.assertEqual(paths['nc'].stat().st_size, original_nc_size)

    def test_full_calendar_year_corrected_total(self):
        # --- Each monthly TDI coefficient is Tg per 365-day year, not Tg/month.
        import contextlib
        import io
        self.obj.yr_s = 2001
        self.obj.yr_e = 2001
        self.obj.ntime = 12
        with self.source.open('w', encoding='ascii') as out:
            for r in range(1, 55):
                out.write(f'land_{r:02d}\n')
                out.write('12 , 2\n')
                for month in range(1, 13):
                    out.write(f'{2001 + (month - 0.5) / 12:.6f} {2.0 if r == 1 else 0.0:.9e} 0.0\n')
        capture = io.StringIO()
        with contextlib.redirect_stdout(capture):
            self.obj.run()
        annual = pd.read_csv(self.obj.get_paths()['diag'] / 'annual_budget.csv')
        self.assertNotIn('months', annual.columns)
        self.assertAlmostEqual(annual.CORR_Tg_year.iloc[0], 2.0, places=6)
        self.assertAlmostEqual(annual.POST_Tg_year.iloc[0], annual.APR_Tg_year.iloc[0] + 2., places=6)
        self.assertAlmostEqual(annual.POST_Tg_year.iloc[0], annual.GT3_Tg_year.iloc[0], places=6)
        self.assertIn('2001', capture.getvalue())
        self.assertIn('Core 2000-2024 (complete years only)', capture.getvalue())

    def test_coordinate_indices(self):
        with xr.open_dataset(self.prior) as ds:
            ilat, ilon = self.obj.check_prior_nc(ds)
            self.assertEqual(int(ilat[0]), 179)
            self.assertEqual(int(ilon[0]), 180)
            self.assertEqual(int(ilon[-1]), 179)

    def test_prior_index_1999(self):
        self.obj.yr_s = 1999
        self.obj.yr_e = 1999
        self.make_sources(year=1999)
        self.obj.run()
        with xr.open_dataset(self.obj.get_paths()['nc']) as ds:
            post = ds['flux_ch4_post'].isel(time=0).values
            np.testing.assert_allclose(post[:90, :], 4.e-13, rtol=1.e-6)

    def test_multiple_direct_access_grid_records(self):
        # --- Legacy Fortran reads REC=1, not the entire file.
        mask_path = self.obj.get_paths()['mask']
        basis_path = self.obj.get_paths()['basis']
        with mask_path.open('ab') as f:
            np.full((180, 360), 47.0, dtype='<f4').tofile(f)
        with basis_path.open('ab') as f:
            np.full((180, 360), -999.0, dtype='<f4').tofile(f)
        tdi = self.obj.read_sources(self.source)
        region, basis, area_lat, totals, counts = self.obj.read_regions(self.obj.get_paths(), tdi)
        self.assertEqual(int(counts[1]), 90 * 360)
        self.assertTrue(np.all(basis == np.float32(1.e-13)))
        self.obj.run()
        self.assertTrue(self.obj.get_paths()['gt3'].is_file())

    def test_mask80_non_inverted_ids_are_ignored(self):
        # --- The original df5_ch4_grd.f90 matches only k=1..54 and ignores
        #     all other IDs in mask80.bin, including 55..80.
        mask_path = self.obj.get_paths()['mask']
        mask = np.zeros((180, 360), dtype=np.float32)
        mask[:30, :] = 55
        mask[30:60, :] = 80
        mask[60:90, :] = -999
        mask[90:, :] = 1
        mask.tofile(mask_path)
        tdi = self.obj.read_sources(self.source)
        region, basis, area_lat, totals, counts = self.obj.read_regions(self.obj.get_paths(), tdi)
        self.assertEqual(int(counts[1]), 90 * 360)
        self.assertEqual(int(counts[0]), 90 * 360)
        self.assertTrue(np.all(region[:90, :] == 0))
        self.assertTrue(np.all(region[90:, :] == 1))
        self.obj.run()
        with xr.open_dataset(self.obj.get_paths()['nc']) as ds:
            post = ds['flux_ch4_post'].isel(time=0).values
            # --- Verify untouched prior for ignored mask IDs.
            with xr.open_dataset(self.prior) as prior_ds:
                ilat, ilon = self.obj.check_prior_nc(prior_ds)
                prior = prior_ds['flux_ch4_prior'].isel(time=24).values
                prior = prior[np.ix_(ilat, ilon)]
            np.testing.assert_array_equal(post[:90], prior[:90])
            self.assertGreater(float(np.mean(post[90:])), float(np.mean(prior[90:])))

    def test_mask_without_any_land_regions_rejected(self):
        mask_path = self.obj.get_paths()['mask']
        np.full((180, 360), 80, dtype='<f4').tofile(mask_path)
        tdi = self.obj.read_sources(self.source)
        with self.assertRaisesRegex(ValueError, 'no recognized regions'):
            self.obj.read_regions(self.obj.get_paths(), tdi)

    def test_reject_truncated_direct_access_grid(self):
        basis_path = self.obj.get_paths()['basis']
        with basis_path.open('ab') as f:
            f.write(b'bad')
        tdi = self.obj.read_sources(self.source)
        with self.assertRaisesRegex(ValueError, 'direct-access grid'):
            self.obj.read_regions(self.obj.get_paths(), tdi)

    def test_real_tdi_uploaded_block_format(self):
        # --- Source includes D4 additional sectors after the 54 land regions.
        #     Basic parser must tolerate names and date lines after land_54.
        with self.source.open('a', encoding='ascii') as f:
            f.write('  anth\n')
            f.write('       27 ,           2\n')
            f.write('1999.0000   0.00000E+00   1.00000E-02\n')
        self.assertEqual(self.obj.read_sources(self.source).shape, (2, 54))

    def test_trailing_tdi_records_are_ignored(self):
        # --- Real TDI files may contain records after the 54 land blocks.
        reference = self.obj.read_sources(self.source)
        with self.source.open('a', encoding='ascii') as out:
            out.write('extra_tdi_record\n')
            out.write('land_55\n')
            out.write('    2 ,          2\n')
            out.write('2000.041667 123.0\n')
        actual = self.obj.read_sources(self.source)
        np.testing.assert_array_equal(reference, actual)

    def test_missing_tdi_record_is_still_rejected(self):
        # --- Never hide truncation or corruption inside the 54 required regions.
        lines = self.source.read_text(encoding='ascii').splitlines()
        self.source.write_text('\n'.join(lines[:-1]) + '\n', encoding='ascii')
        with self.assertRaisesRegex(ValueError, 'Truncated D4 month'):
            self.obj.read_sources(self.source)

    def test_zero_predicted_rejected(self):
        self.make_sources(nonzero=0)
        with self.assertRaisesRegex(ValueError, 'ZERO predicted.sources'):
            self.obj.read_sources(self.source)

    def test_calendar_mismatch_rejected(self):
        self.make_sources(year=1999)
        with self.assertRaisesRegex(ValueError, 'calendar mismatch'):
            self.obj.read_sources(self.source)

    def test_missing_prior_does_not_create_outputs(self):
        self.obj.flx_apr_nc = 'does_not_exist.nc'
        with self.assertRaises(FileNotFoundError):
            self.obj.run()
        self.assertFalse(self.result_root.exists() and any(self.result_root.rglob('*.gt3')))

    def test_dry_run(self):
        self.obj.run(dry_run=True)
        self.assertFalse(self.obj.get_paths()['nc'].exists())

    def test_paths_have_netcdf_no_prior_grd(self):
        paths = self.obj.get_paths()
        self.assertEqual(paths['prior'].name, 'gcp25_inv1.nc')
        self.assertEqual(paths['nc'].name, 's057_ctl_inv1_noLC.nc')
        self.assertEqual(paths['gt3'].name, 's057_ctl_inv1_noLC.gt3')
        self.assertNotIn('grd', paths)
        self.assertNotIn('loss', paths)

    def test_reject_symlink(self):
        link = self.work / 'symlink.nc'
        link.symlink_to(self.prior)
        with self.assertRaisesRegex(ValueError, 'symlink'):
            self.obj.check_file(link)

    def test_namelist_quote(self):
        self.assertEqual(self.obj.quote('/tmp/a'), "'/tmp/a'")
        with self.assertRaises(ValueError):
            self.obj.quote("/tmp/it's")


if __name__ == '__main__':
    unittest.main(verbosity=2)
