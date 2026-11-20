# -*- coding: utf-8 -*-
"""D3 offline tests: python test/test_d3.py"""
import importlib.util
import subprocess
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
# In a standalone downloaded bundle, use a minimal SetCase double.
if not (ROOT / '_set_case.py').exists():
    import types
    settings = types.ModuleType('_set_case')
    class SetCase:
        def __init__(self):
            self.hcase = os.environ.get('GCP_CASE', 'INCA')
    settings.SetCase = SetCase
    sys.modules['_set_case'] = settings
spec = importlib.util.spec_from_file_location('d3_obs_actm', ROOT / 'd3_obs_actm.py')
d3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d3)


class D3Tests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.obj = object.__new__(d3.ObsACTM)
        self.obj.root = self.root
        self.obj.bin_dir = self.root / 'bin'
        self.obj.config_dir = self.root / 'config'
        self.obj.case = 'INCA'
        self.obj.yr_s, self.obj.yr_e = 1999, 2025
        self.obj.inv_wrk_dir = str(self.root / 'results')
        self.obj.d3_diagnostics_dir = str(self.root / 'results' / 'diagnostics')
        self.obj.d3_sample_output_prefix = 'miroc_ch4sfc_inca_'
        self.obj.d3_combine_input_prefix = self.obj.d3_sample_output_prefix
        self.obj.d3_combine_output_prefix = 'obsrv_actm_inca_'
        self.obj.unxcases = ['ctl', 'ux2', 'ux4']
        self.obj.invcases = ['inv1']
        self.obj.sites_f = str(self.root / 'sites.txt')
        Path(self.obj.sites_f).write_text('station\n')
        self.obj.d3_sample_model_dir = str(self.root / 'ACTM_INCA')
        Path(self.obj.d3_sample_model_dir).mkdir()
        self.obj.d3_sample_output_dir = str(self.root / 'actm')
        self.obj.d3_combine_model_dir = self.obj.d3_sample_output_dir
        self.obj.d3_combine_obs_dir = str(self.root / 'obs')
        Path(self.obj.d3_combine_obs_dir).mkdir()
        self.obj.d3_combine_output_dir = str(self.root / 'joined')
        self.obj.d3_filter_input_file = str(self.root / 'joined/obsrv_actm_inca_ch4dy.grd')
        self.obj.d3_filter_output_file = str(self.root / 'ch4obs_actm_inca_fits.bin')
        self.obj.d3_observations_input_file = self.obj.d3_filter_output_file
        self.obj.d3_observations_summary = str(self.root / 'obsrvCH4_2025_inca.txt')
        self.obj.d3_observations_prior = str(self.root / 'prior_data_inca.dif')
        self.obj.d3_observations_nc_prefix = str(self.root / 'obsrvCH4_inca_9925')

    def test_01_sample_namelist(self):
        cfg = self.obj.write_config('sample', {'s_year': 1999, 'f_year': 2025, 'o_fil': 'miroc_ch4sfc_inca_'})
        value = cfg.read_text()
        self.assertIn('&d3_sample_actm_sites', value)
        self.assertIn('s_year = 1999', value)
        self.assertIn("hcase = 'INCA'", value)

    def test_02_reject_quote(self):
        with self.assertRaises(ValueError):
            self.obj.write_config('combine', {'i_dir': "bad'path"})

    def test_03_all_stage_names(self):
        names = {
            'sample': 'd3_sample_actm_sites',
            'combine': 'd3_combine_obs_actm',
            'filter': 'd3_filter_ch4_nakazawa',
            'observations': 'd3_write_tdi_observations',
        }
        self.assertEqual(self.obj.NAMES, names)

    def nml_value(self, config, name):
        for line in config.read_text().splitlines():
            if line.strip().startswith(name + ' ='):
                return line.split('=', 1)[1].strip().strip("'")
        self.fail('Missing namelist field ' + name)

    def create_observations_products(self, config, content=b'new'):
        summary = Path(self.nml_value(config, 'o_obsrv'))
        prior = Path(self.nml_value(config, 'o_prior'))
        prefix = Path(self.nml_value(config, 'o_pr_nc'))
        for path in self.obj.observation_products(summary, prior, prefix):
            path.write_bytes(content)

    def test_04_sample_outputs(self):
        def fake(stage, config, work_dir=None):
            for year in range(1999, 2026):
                folder = Path(self.obj.d3_sample_output_dir)
                folder.mkdir(exist_ok=True)
                (folder / (self.obj.d3_sample_output_prefix + str(year) + '.grd')).write_bytes(b'123')
        with patch.object(self.obj, 'run_exe', side_effect=fake):
            self.obj.get_actm_at_sites()

    def test_05_sample_missing_year(self):
        with patch.object(self.obj, 'run_exe'):
            with self.assertRaises(FileNotFoundError):
                self.obj.get_actm_at_sites()

    def test_06_combine_outputs(self):
        folder = Path(self.obj.d3_sample_output_dir)
        folder.mkdir()
        for year in range(1999, 2026):
            (folder / (self.obj.d3_sample_output_prefix + str(year) + '.grd')).write_bytes(b'abc')
        def fake(stage, config, work_dir=None):
            out = Path(self.obj.d3_combine_output_dir)
            for suffix in ('ch4yr.grd', 'ch4mm.grd', 'ch4dy.grd'):
                (out / (self.obj.d3_combine_output_prefix + suffix)).write_bytes(b'abc')
        with patch.object(self.obj, 'run_exe', side_effect=fake):
            self.obj.combine_obs_actm()

    def test_07_filter_replaces_existing_with_backup(self):
        inp = Path(self.obj.d3_filter_input_file)
        inp.parent.mkdir(parents=True)
        inp.write_bytes(b'input')
        out = Path(self.obj.d3_filter_output_file)
        out.write_bytes(b'old')

        def fake(stage, config, work_dir=None):
            self.assertEqual(stage, 'filter')
            staged = Path(self.nml_value(config, 'o_file'))
            self.assertNotEqual(staged, out)
            staged.write_bytes(b'new')

        with patch.object(self.obj, 'run_exe', side_effect=fake), \
             patch.object(self.obj, 'check_filter_output'):
            self.obj.filter_ch4()
        self.assertEqual(out.read_bytes(), b'new')
        backups = list((self.root / 'results/backups').rglob('*.bin'))
        self.assertEqual(len(backups), 1)
        self.assertEqual(backups[0].read_bytes(), b'old')

    def test_08_observations_replace_existing_with_backup(self):
        Path(self.obj.d3_observations_input_file).write_bytes(b'abc')
        prior = Path(self.obj.d3_observations_prior)
        prior.write_bytes(b'old')

        def fake(stage, config, work_dir=None):
            self.assertEqual(stage, 'observations')
            self.create_observations_products(config)

        with patch.object(self.obj, 'run_exe', side_effect=fake):
            self.obj.write_tdi_observations()
        self.assertEqual(prior.read_bytes(), b'new')
        backups = list((self.root / 'results/backups').rglob('*.dif'))
        self.assertEqual(len(backups), 1)
        self.assertEqual(backups[0].read_bytes(), b'old')

    def test_09_filter_namelist(self):
        inp = Path(self.obj.d3_filter_input_file)
        inp.parent.mkdir(parents=True)
        inp.write_bytes(b'abc')

        def fake(stage, config, work_dir=None):
            self.assertIn('i_file =', config.read_text())
            self.assertIn('o_file =', config.read_text())
            Path(self.nml_value(config, 'o_file')).write_bytes(b'new')

        with patch.object(self.obj, 'run_exe', side_effect=fake), \
             patch.object(self.obj, 'check_filter_output'):
            self.obj.filter_ch4()
        self.assertEqual(Path(self.obj.d3_filter_output_file).read_bytes(), b'new')
        stable = self.obj.config_dir / 'd3_filter_ch4_nakazawa.nml'
        self.assertEqual(self.nml_value(stable, 'o_file'), self.obj.d3_filter_output_file)

    def test_10_observations_products(self):
        Path(self.obj.d3_observations_input_file).write_bytes(b'abc')

        def fake(stage, config, work_dir=None):
            value = config.read_text()
            self.assertIn('cfg_s_year = 1999', value)
            self.assertIn('cfg_f_year = 2025', value)
            self.create_observations_products(config, b'good')

        with patch.object(self.obj, 'run_exe', side_effect=fake):
            self.obj.write_tdi_observations()
        output = Path(self.obj.d3_observations_prior)
        self.assertEqual(output.read_bytes(), b'good')
        products = self.obj.observation_products(
            Path(self.obj.d3_observations_summary), output,
            Path(self.obj.d3_observations_nc_prefix))
        self.assertEqual(len(products), 9)
        self.assertTrue(all(p.is_file() for p in products))

    def test_10b_filter_failure_keeps_old_file(self):
        inp = Path(self.obj.d3_filter_input_file)
        inp.parent.mkdir(parents=True)
        inp.write_bytes(b'input')
        out = Path(self.obj.d3_filter_output_file)
        out.write_bytes(b'original')

        def fake(stage, config, work_dir=None):
            Path(self.nml_value(config, 'o_file')).write_bytes(b'broken')
            raise subprocess.CalledProcessError(1, ['fortran'])

        with patch.object(self.obj, 'run_exe', side_effect=fake):
            with self.assertRaises(subprocess.CalledProcessError):
                self.obj.filter_ch4()
        self.assertEqual(out.read_bytes(), b'original')

    def test_10c_missing_tdi_product_keeps_old_outputs(self):
        Path(self.obj.d3_observations_input_file).write_bytes(b'input')
        prior = Path(self.obj.d3_observations_prior)
        prior.write_bytes(b'original')

        def fake(stage, config, work_dir=None):
            Path(self.nml_value(config, 'o_prior')).write_bytes(b'partial')

        with patch.object(self.obj, 'run_exe', side_effect=fake):
            with self.assertRaises(FileNotFoundError):
                self.obj.write_tdi_observations()
        self.assertEqual(prior.read_bytes(), b'original')

    def test_10d_filter_binary_length(self):
        path = self.root / 'real_output.bin'
        with path.open('wb') as f:
            f.truncate(9862 * 245 * 2 * 3 * 4)
        self.obj.check_filter_output(path)
        with path.open('r+b') as f:
            f.truncate(128)
        with self.assertRaises(ValueError):
            self.obj.check_filter_output(path)

    def test_11_order(self):
        events = []
        with patch.object(self.obj, 'get_actm_at_sites', side_effect=lambda: events.append('sample')), \
             patch.object(self.obj, 'combine_obs_actm', side_effect=lambda: events.append('combine')), \
             patch.object(self.obj, 'filter_ch4', side_effect=lambda: events.append('filter')), \
             patch.object(self.obj, 'write_tdi_observations', side_effect=lambda: events.append('observations')), \
             patch.object(self.obj, 'diagnostics', side_effect=lambda: events.append('diagnostics')):
            self.obj.run()
        self.assertEqual(events, ['sample', 'combine', 'filter', 'observations', 'diagnostics'])

    def test_12_fortran_groups_match(self):
        for operation, name in self.obj.NAMES.items():
            source = ROOT / 'fortran' / ('df3_' + name[3:] + '.f90')
            self.assertTrue(source.is_file(), str(source))
            code = source.read_text()
            self.assertIn('namelist /' + name + '/', code.lower())
            self.assertNotIn('namelist /actm_f', code.lower())

    def test_13_filter_libraries(self):
        for name in ('lsfbf4_main.f', 'lsfbf4_sub.f'):
            self.assertTrue((ROOT / 'fortran' / name).is_file())

    def test_14_compile_command_sources(self):
        sources = [ROOT / 'fortran' / 'df3_filter_ch4_nakazawa.f90',
                   ROOT / 'fortran' / 'lsfbf4_main.f',
                   ROOT / 'fortran' / 'lsfbf4_sub.f']
        with patch.object(d3.subprocess, 'run') as subprocess_mock, \
             patch.object(d3.shutil, 'which', return_value='/usr/bin/gfortran'):
            self.obj.build_executable('filter', sources, self.root / 'bin' / 'df3_filter_ch4_nakazawa')
        args = subprocess_mock.call_args.args[0]
        for source in sources:
            self.assertIn(str(source), args)



    def test_15_diagnostics_path_and_call(self):
        self.obj.d3_filter_input_file = str(self.root / 'joined/raw.grd')
        self.obj.d3_filter_output_file = str(self.root / 'fits.bin')
        with patch.object(d3, 'generate_d3_diagnostics', return_value={'ok': True}) as func:
            self.assertEqual(self.obj.diagnostics(), {'ok': True})
        self.assertEqual(func.call_args.args[0], 'INCA')
        self.assertEqual(func.call_args.kwargs['outdir'], self.root / 'results/diagnostics')
        self.assertEqual(func.call_args.kwargs['plots'], True)

    def test_16_diagnostics_rejects_wrong_period(self):
        self.obj.yr_s = 2000
        with self.assertRaises(ValueError):
            self.obj.diagnostics()

    def test_17_diagnostic_plot_writes_png(self):
        import numpy as np
        import pandas as pd
        days = pd.date_range('2001-01-01', '2001-12-31', freq='D')
        months = pd.date_range('2001-01-01', '2001-12-01', freq='MS')
        count = len(days)
        inp = np.full((count, 2, 1), 1800.0, dtype='f4')
        inp[:, 1, :] = 1790.0
        fits = np.full((count, 2, 3, 1), 1800.0, dtype='f4')
        fits[:, 1, :, :] = 1790.0
        fits[:, :, 2, :] = 1790.0
        monthly = {'fit_diff': np.full((12, 1), 10.0),
                   'paired_raw_days': np.full((12, 1), 25),
                   'paired_fit_days': np.full((12, 1), 25)}
        sites = pd.DataFrame([{'station_id': 1, 'site': 'SYN', 'lat': 35., 'lon': 140.,
                               'active': True}])
        out = self.root / 'results' / 'diagnostics'
        out.mkdir(parents=True)
        path = d3.plot_station(out, 0, sites, days, months, inp, fits, monthly, {}, 'INCA')
        self.assertEqual(path.name, 'station_001_SYN.png')
        self.assertGreater(path.stat().st_size, 10000)

    def test_18_binary_reader_checks_file_size(self):
        out = self.root / 'binary.bin'
        out.write_bytes(b'01234567')
        data = d3.get_binary(out, (2,), endian='<')
        self.assertEqual(len(data), 2)
        with self.assertRaises(ValueError):
            d3.get_binary(out, (3,), endian='<')


    def test_19_full_diagnostics_synthetic(self):
        import numpy as np
        from scipy.io import netcdf_file
        from types import SimpleNamespace

        inv = self.root / 'inv_dir'
        inv.mkdir()
        site_path = inv / 'obsrvCH4_60sites.txt'
        site_path.write_text('header one\nheader two\n' + ''.join(
            f'{n:5d} {9999 if n <= 2 else 0:5d} ch4_test{n:03d} 50.0 10.0 100.0\n'
            for n in range(1, d3.N_SITE + 1)))
        input_file = inv / 'obsrv_actm_inca_ch4dy.grd'
        output_file = inv / 'ch4obs_actm_inca_fits.bin'
        raw = np.memmap(input_file, dtype='<f4', mode='w+', shape=(d3.N_DAY, 2, d3.N_SITE))
        fitted = np.memmap(output_file, dtype='<f4', mode='w+',
                           shape=(d3.N_DAY, 2, 3, d3.N_SITE))
        raw[:] = d3.FILL
        fitted[:] = d3.FILL
        raw[:, 0, 0], raw[:, 1, 0] = 1850.0, 1840.0
        fitted[:, 0, 0, 0], fitted[:, 1, 0, 0] = 1850.0, 1840.0
        fitted[:, 0, 1, 0], fitted[:, 1, 1, 0] = 1849.0, 1839.0
        fitted[:, 0, 2, 0], fitted[:, 1, 2, 0] = 1830.0, 1820.0
        raw.flush()
        fitted.flush()
        del raw, fitted

        days, months = d3.get_dates()
        _, fits = d3.get_filter_files(input_file, output_file)
        monthly = d3.get_filter_monthly(fits, days, months)
        expected = d3.get_expected_tdi(monthly, 1, months)
        self.assertTrue(np.allclose(expected[:, 0, 0], 10.0))
        self.assertTrue(np.allclose(expected[:, 0, 1], 6.0))
        sites = d3.get_sites(site_path)
        nc_path = inv / 'obsrvCH4_inca_9925ctlinv1.nc'
        with netcdf_file(str(nc_path), 'w') as nc:
            for name, size in [('station_count', d3.N_SITE), ('tlen', d3.N_MONTH),
                               ('data_type', 2), ('namelen', 45)]:
                nc.createDimension(name, size)
            nc.createVariable('station_count', 'i2', ('station_count',))[:] = sites.station_id.to_numpy(dtype='i2')
            labels = nc.createVariable('site_name', 'c', ('station_count', 'namelen'))
            raw_names = np.full((d3.N_SITE, 45), b' ', dtype='S1')
            for i, name in enumerate(sites.site):
                raw_names[i, :len(name)] = np.frombuffer(name.encode('ascii'), dtype='S1')
            labels[:] = raw_names
            nc.createVariable('tlen', 'f4', ('tlen',))[:] = np.asarray(
                1999 + (np.arange(d3.N_MONTH) + 0.5) / 12.0, dtype='f4')
            nc.createVariable('tracer_data', 'f4', ('tlen', 'station_count', 'data_type'))[:] = expected.astype('f4')
        nc = d3.get_tdi_netcdf(nc_path)
        stats = d3.get_tdi_stats(nc, sites)
        self.assertEqual(stats.iloc[0]['n_valid'], d3.N_MONTH)
        cmp, _ = d3.get_comparison(nc, expected, sites, months)
        self.assertEqual(cmp['n_total_bad_months'], 0)

        setup = SimpleNamespace(sites_f=str(site_path), inv_wrk_dir=str(self.root / 'results'),
                                d3_filter_input_file=str(input_file),
                                d3_filter_output_file=str(output_file),
                                d3_observations_nc_prefix=str(inv / 'obsrvCH4_inca_9925'))
        result = d3.generate_d3_diagnostics('INCA', setup, plots=True, stations=[1, 42])
        self.assertEqual(result['station_ids'], [1, 2])
        self.assertEqual(result['comparison_errors'], 0)
        images = list((self.root / 'results/diagnostics').glob('station_*.png'))
        self.assertEqual(len(images), 2)
        self.assertTrue((self.root / 'results/diagnostics/ch4_station_stats.txt').is_file())

if __name__ == '__main__':
    unittest.main(verbosity=2)
