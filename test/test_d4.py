# -*- coding: utf-8 -*-
"""Autonomous D4 Python tests with synthetic NetCDF and a fake TDI executable.

They test orchestration, not the numerical Fortran inversion.
"""
import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from d4_inversion import Inversion


# ============================================================
# === Create a minimal isolated TDI project with realistic filenames
class TestD4(unittest.TestCase):

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        previous = os.environ.get('GCP_CASE')
        os.environ['GCP_CASE'] = 'INCA'
        self.addCleanup(self.restore_env, previous)
        self.obj = Inversion(case='INCA')
        self.obj.d4_plot_all_regions = False
        a = self.obj
        a.result_root = self.root / 'results2025' / 'INCA'
        a.diag_dir = a.result_root / 'diagnostics' / 'd4'
        a.inv_dir = str(self.root / 'inv_dir') + '/'
        a.tdi_dir = self.root / 'inv_dir' / 'tdi_inv'
        a.files_dir = a.tdi_dir / 'files4run'
        a.gf_dir = self.root / 'inv_dir' / 'gmats'
        a.template = a.tdi_dir / 'control.example'
        a.executable = a.tdi_dir / 'tdi_ch4'
        a.d1_uncertainty_output_prefix = str(self.root / 'inv_dir' / 'sourceprior' / 'sourcepriorch4_p')
        a.d3_observations_nc_prefix = str(self.root / 'inv_dir' / 'obsrvCH4_inca_9925')
        a.sites_f = str(self.root / 'inv_dir' / 'obsrvCH4_60sites.txt')
        a.files_dir.mkdir(parents=True)
        a.gf_dir.mkdir(parents=True)
        (a.gf_dir / 'GF_01.dat').write_bytes(b'GF')
        (a.gf_dir / 'ACTM.ch4j55.2011.L2.nc').write_bytes(b'GF')
        Path(a.sites_f).write_text('test site\n')
        static = {'srccorf': 'spatial-srccor', 'outputf': 'global-output',
                  'groupf': 'spatial-groups', 'constf': 'spatial-consts-obs',
                  'recipef': 'new-recipes', 'greenf': 'global-block'}
        content = '&control\n'
        for name, value in static.items():
            (a.files_dir / value).write_text('test\n')
            content += ' ' + name + " = 'files4run/" + value + "'\n"
        content += " sourcef = 'old.nc'\n datf = 'old_obs.nc'\n datflist = 'old_sites.txt'\n"
        content += ' firstsrc=1999\n lastsrc=2025\n firstdat=1999\n lastdat=2025\n'
        content += ' tc3tdi=.true.\n pseudodata=.false.\n/\n'
        a.template.write_text(content)
        a.executable.write_text(
            '#!/bin/sh\n'
            + sys.executable + " - <<'PYCODE'\n"
            + "from pathlib import Path\n"
            + "Path('prior.sources').write_text((' '.join(['0.0'] * 54) + '\\n') * 324)\n"
            + "Path('predicted.sources').write_text((' '.join(['0.1'] * 54) + '\\n') * 324)\n"
            + "PYCODE\n"
            + "echo 'TDI finished'\n")
        a.executable.chmod(0o755)
        source, obs, self.run_dir = a.get_paths('p30', 'ctl', 'inv1')
        source.parent.mkdir(parents=True, exist_ok=True)
        obs.parent.mkdir(parents=True, exist_ok=True)
        values = np.zeros((324, 61, 2), dtype=np.float32)
        values[..., 0] = 2.0
        values[..., 1] = 8.0
        values[2, 3, 0] = -999.0
        xr.Dataset({'tracer_data': (('tlen', 'station_count', 'data_type'), values)}).to_netcdf(obs, engine='scipy')
        flux = np.ones((348, 54, 2), dtype=np.float32)
        xr.Dataset({'sources_land': (('month', 'land_region', 'datatype_count'), flux),
                    'time': ('month', np.linspace(1998.041667, 2026.958333, 348))}).to_netcdf(source, engine='scipy')
        self.source, self.obs = source, obs

    @staticmethod
    def restore_env(previous):
        if previous is None:
            os.environ.pop('GCP_CASE', None)
        else:
            os.environ['GCP_CASE'] = previous

    def test_greenf_shared_directory_with_trailing_separator(self):
        gf = self.obj.gf_dir / 'ACTM.ch4j55.2011.L2.nc'
        gf.write_bytes(b'GF')
        self.obj.prepare_inversion('p30', 'ctl', 'inv1', self.source, self.obs, self.run_dir)
        prepared = (self.run_dir / 'control').read_text()
        value = self.obj.get_control_value(prepared, 'greenf')
        self.assertEqual(value, str(self.obj.gf_dir.resolve()) + os.sep)
        self.assertTrue(value.endswith(os.sep))
        self.assertTrue(Path(value + 'ACTM.ch4j55.2011.L2.nc').samefile(gf))
        self.assertFalse((self.run_dir / 'gmats').exists())
        self.assertFalse(any(p.is_symlink() for p in self.run_dir.rglob('*')))

    def test_greenf_old_template_file_is_overridden_by_shared_directory(self):
        # --- The template may still use a legacy GF filename/global-block.
        self.obj.prepare_inversion('p30', 'ctl', 'inv1', self.source, self.obs, self.run_dir)
        prepared = (self.run_dir / 'control').read_text()
        value = self.obj.get_control_value(prepared, 'greenf')
        self.assertEqual(value, str(self.obj.gf_dir.resolve()) + os.sep)
        self.assertNotIn("greenf = 'files4run/global-block'", prepared)

    def test_greenf_over_80_characters_rejected(self):
        long_dir = self.root / ('a' * 75)
        long_dir.mkdir()
        self.obj.gf_dir = long_dir
        with self.assertRaisesRegex(ValueError, r'character\*80'):
            self.obj.localize_control_file(self.run_dir, 'greenf', 'dummy')

    def test_all_zero_decoded_predicted_sources_is_failure(self):
        self.obj.run()
        result = self.run_dir / 'predicted.sources'
        with result.open('w') as f:
            for r in range(1, 55):
                f.write(f'land_{r:02d}\n')
                f.write('decimal year, correction\n')
                for k in range(324):
                    time = 1999 + k // 12 + ((k % 12) + 0.5) / 12
                    f.write(f'{time:.8f} 0.0\n')
        status, problems = self.obj.check_run(self.run_dir, strict=False)
        self.assertEqual(status, 'FAIL')
        self.assertIn('all_zero_TDI_corrections_CHECK_SOLVER_AND_GREENS', problems)
        with self.assertRaisesRegex(RuntimeError, 'all_zero_TDI_corrections'):
            self.obj.check_run(self.run_dir, strict=True)

    def test_case(self):
        self.assertEqual(self.obj.case, 'INCA')

    def test_period(self):
        self.assertEqual(self.obj.ntime, 324)

    def test_paths(self):
        self.assertIn('s057_p30', str(self.run_dir))
        self.assertIn('s057_ctl_inv1', str(self.run_dir))

    def test_inputs(self):
        info = self.obj.check_inputs('p30', 'ctl', 'inv1')[-1]
        self.assertEqual(info['observations_valid_station_months'], 324 * 61 - 1)

    def test_control(self):
        self.obj.prepare_inversion('p30', 'ctl', 'inv1', self.source, self.obs, self.run_dir)
        text = (self.run_dir / 'control').read_text()
        self.assertIn("sourcef = 'sourcepriortdi.nc'", text)
        self.assertIn("greenf = '" + str(self.obj.gf_dir.resolve()) + os.sep + "'", text)
        self.assertFalse((self.run_dir / 'gmats').exists())
        self.assertTrue((self.run_dir / 'files4run').is_dir())
        self.assertFalse((self.run_dir / 'files4run').is_symlink())
        for name in ('sourcepriortdi.nc', 'obsrvCH4_9925.nc', 'obsrvCH4_sites.txt', 'tdi_ch4'):
            self.assertTrue((self.run_dir / name).is_file())
            self.assertFalse((self.run_dir / name).is_symlink())
        self.assertEqual(list(p for p in self.run_dir.rglob('*') if p.is_symlink()), [])

    def test_invalid_control_period(self):
        txt = self.obj.template.read_text().replace('firstdat=1999', 'firstdat=2000')
        self.obj.template.write_text(txt)
        with self.assertRaises(ValueError):
            self.obj.prepare_inversion('p30', 'ctl', 'inv1', self.source, self.obs, self.run_dir)
        # Invalid scientific periods are rejected, not silently replaced.

    def test_legacy_limit(self):
        with self.assertRaises(ValueError):
            self.obj.replace_control_value(self.obj.template.read_text(), 'sourcef', 'z' * 81)

    def test_missing_obs(self):
        self.obs.unlink()
        with self.assertRaises(FileNotFoundError):
            self.obj.check_inputs('p30', 'ctl', 'inv1')

    def test_empty_gf(self):
        (self.obj.gf_dir / 'GF_01.dat').unlink()
        (self.obj.gf_dir / 'ACTM.ch4j55.2011.L2.nc').unlink()
        with self.assertRaises(ValueError):
            self.obj.check_inputs('p30', 'ctl', 'inv1')

    def test_missing_tc3greens_expected_files_is_rejected(self):
        (self.obj.gf_dir / 'ACTM.ch4j55.2011.L2.nc').unlink()
        with self.assertRaisesRegex(FileNotFoundError, 'ACTM.ch4j55'):
            self.obj.check_inputs('p30', 'ctl', 'inv1')

    def test_old_run_with_symlinks_rejected(self):
        self.obj.run()
        (self.run_dir / 'files4run' / 'global-block').unlink()
        (self.run_dir / 'files4run' / 'global-block').symlink_to(self.obj.files_dir / 'global-block')
        with self.assertRaisesRegex(FileExistsError, 'contains symlinks'):
            self.obj.run_inversion('p30', 'ctl', 'inv1', self.source, self.obs, self.run_dir)

    def test_no_symlinks_and_copies_independent(self):
        self.obj.prepare_inversion('p30', 'ctl', 'inv1', self.source, self.obs, self.run_dir)
        assert not any(p.is_symlink() for p in self.run_dir.rglob('*'))
        source_bytes = (self.run_dir / 'sourcepriortdi.nc').read_bytes()
        obs_bytes = (self.run_dir / 'obsrvCH4_9925.nc').read_bytes()
        exe_bytes = (self.run_dir / 'tdi_ch4').read_bytes()
        self.source.unlink()
        self.obs.unlink()
        (self.obj.gf_dir / 'GF_01.dat').unlink()
        self.obj.executable.unlink()
        self.assertEqual((self.run_dir / 'sourcepriortdi.nc').read_bytes(), source_bytes)
        self.assertEqual((self.run_dir / 'obsrvCH4_9925.nc').read_bytes(), obs_bytes)
        self.assertFalse((self.run_dir / 'gmats').exists())
        self.assertEqual((self.run_dir / 'tdi_ch4').read_bytes(), exe_bytes)

    def test_execution(self):
        self.obj.run()
        self.assertEqual(self.obj.check_run(self.run_dir)[0], 'PASS')
        self.assertTrue((self.run_dir / 'RUN_INPUTS.txt').is_file())
        self.assertTrue((self.run_dir / 'RUN_STATUS.json').is_file())

    def test_idempotent(self):
        self.obj.run()
        self.assertEqual(self.obj.run_inversion('p30', 'ctl', 'inv1',
                         self.source, self.obs, self.run_dir), 'REUSED')

    def test_modified_input_refuses_reuse(self):
        self.obj.run()
        data = self.obs.read_bytes() + b'X'
        self.obs.write_bytes(data)
        with self.assertRaises(FileExistsError):
            self.obj.run_inversion('p30', 'ctl', 'inv1', self.source, self.obs, self.run_dir)

    def test_incomplete_refuses_overwrite(self):
        self.run_dir.mkdir(parents=True)
        with self.assertRaises(FileExistsError):
            self.obj.run_inversion('p30', 'ctl', 'inv1', self.source, self.obs, self.run_dir)

    def test_failed_run_preserves_workdir(self):
        self.obj.executable.write_text('#!/bin/sh\nexit 1\n')
        self.obj.executable.chmod(0o755)
        import subprocess
        with self.assertRaises(subprocess.CalledProcessError):
            self.obj.run()
        self.assertTrue((self.run_dir / 'RUN_STATUS.json').is_file())
        self.assertTrue((self.run_dir / 'control').is_file())

    def test_missing_output_fails(self):
        self.obj.run()
        (self.run_dir / 'predicted.sources').unlink()
        status, issues = self.obj.check_run(self.run_dir, strict=False)
        self.assertEqual(status, 'FAIL')
        self.assertIn('missing_or_empty_predicted.sources', issues)

    def test_unreadable_coefficients_fail_numerical_qc(self):
        self.obj.run()
        (self.run_dir / 'predicted.sources').write_text('not a regional inversion result\n')
        status, problems = self.obj.check_run(self.run_dir, strict=False)
        self.assertEqual(status, 'FAIL')
        self.assertIn('unrecognized_predicted_sources_format', problems)
        with self.assertRaisesRegex(RuntimeError, 'NUMERICAL QUALITY FAILED'):
            self.obj.diagnostics()

    def test_diagnostics(self):
        self.obj.run()
        self.obj.diagnostics()
        self.assertTrue((self.obj.diag_dir / 'inversion_runs.csv').is_file())
        self.assertTrue((self.obj.diag_dir / 'observations_timeseries.png').is_file())
        self.assertTrue((self.obj.diag_dir / 'observations_coverage.png').is_file())
        self.assertTrue((self.obj.diag_dir / 'observation_uncertainty.png').is_file())
        self.assertTrue((self.obj.diag_dir / 'prior_uncertainty_regions_p30.png').is_file())
        self.assertFalse((self.obj.diag_dir / 'experiment_inputs.png').exists())
        self.assertTrue((self.obj.diag_dir / 'tdi_log_diagnostics.txt').is_file())
        self.assertTrue((self.obj.diag_dir / 'inversion_overview_INCA_p30_ctl_inv1.png').is_file())

    def test_numeric_source_reader(self):
        arr = np.arange(54 * 324, dtype=float).reshape(324, 54)
        fn = self.root / 'sources.txt'
        np.savetxt(fn, arr)
        got = self.obj.read_source_coefficients(fn)
        self.assertEqual(got.shape, (324, 54))
        self.assertTrue(np.array_equal(got, arr))

    def test_source_plots(self):
        self.obj.run()
        prior = np.ones((324, 54))
        posterior = np.full((324, 54), 1.1)
        np.savetxt(self.run_dir / 'prior.sources', prior)
        np.savetxt(self.run_dir / 'predicted.sources', posterior)
        self.obj.diagnostics()
        self.assertFalse((self.obj.diag_dir / 'source_sensitivity.png').exists())
        self.assertTrue((self.obj.diag_dir / 'source_correction_regions_p30_ctl_inv1.png').is_file())
        self.assertTrue((self.obj.diag_dir / 'inversion_overview_INCA_p30_ctl_inv1.png').is_file())
        self.assertTrue((self.obj.diag_dir / 'regional' / 'p30_ctl_inv1' / 'region_01.png').is_file())

    def test_block_source_plots(self):
        self.obj.run()
        for filename, scale in [('prior.sources', 0.0), ('predicted.sources', 1.0)]:
            with (self.run_dir / filename).open('w') as out:
                for reg in range(1, 55):
                    out.write('land_' + str(reg).zfill(2) + '\n')
                    out.write('decimal year, CH4 correction\n')
                    for k in range(324):
                        t = 1999 + k // 12 + ((k % 12) + 0.5) / 12
                        out.write(f'{t:.8f} {scale * reg:.6f}\n')
        self.obj.diagnostics()
        report = (self.obj.diag_dir / 'diagnostic_report.txt').read_text()
        self.assertIn('Decodable regional datasets: 1', report)
        self.assertFalse((self.obj.diag_dir / 'source_sensitivity.png').exists())
        self.assertTrue((self.obj.diag_dir / 'source_annual_p30_ctl_inv1.csv').is_file())
        self.assertTrue((self.obj.diag_dir / 'regional' / 'p30_ctl_inv1' / 'region_54.png').is_file())


    def test_zero_solution_is_critical_and_no_flat_region_plots(self):
        self.obj.run()
        # --- A completed TDI run with all-zero 324 x 54 output must NOT pass QC.
        np.savetxt(self.run_dir / 'predicted.sources', np.zeros((324, 54)))
        np.savetxt(self.run_dir / 'prior.sources', np.zeros((324, 54)))
        log = next(self.run_dir.glob('tdi_*.out'))
        log.write_text(' DOSVD   = T,\n NUMSVDVECS      =           0,\n'
                       ' nobs, valid ntotdat, excluded =    57  17241  1227\n')
        with self.assertRaisesRegex(RuntimeError, 'NUMERICAL QUALITY FAILED'):
            self.obj.diagnostics()
        report = (self.obj.diag_dir / 'diagnostic_report.txt').read_text()
        self.assertIn('ALL_ZERO_PREDICTED_SOURCES', report)
        self.assertIn('NUMERICAL CRITICAL: 1', report)
        self.assertIn('TDI observation equations: 17241', report)
        self.assertTrue((self.obj.diag_dir / 'source_quality.csv').is_file())
        self.assertTrue((self.obj.diag_dir / 'regional' / 'p30_ctl_inv1' /
                         'CRITICAL_ZERO_SOLUTION.txt').exists())
        self.assertFalse((self.obj.diag_dir / 'regional' / 'p30_ctl_inv1' /
                          'region_01.png').exists())

    def test_cross_case_comparison(self):
        self.obj.run()
        prior = np.full((324, 54), 2.0)
        predicted = np.full((324, 54), 3.0)
        np.savetxt(self.run_dir / 'prior.sources', prior)
        np.savetxt(self.run_dir / 'predicted.sources', predicted)
        peer = (self.obj.result_root.parent / 'CYC' / 's057_p30' /
                's057_ctl_inv1' / 'predicted.sources')
        peer.parent.mkdir(parents=True)
        np.savetxt(peer, np.ones((324, 54)))
        self.obj.diagnostics()
        self.assertTrue((self.obj.diag_dir / 'compare_CYC_INCA_p30_ctl_inv1.png').is_file())
        report = (self.obj.diag_dir / 'diagnostic_report.txt').read_text()
        self.assertIn('Matched CYC/INCA comparisons: 1', report)

    def test_regional_stats_core_period(self):
        self.obj.run()
        predicted = np.zeros((324, 54))
        predicted[12:312, 36] = np.arange(300) / 100.0
        np.savetxt(self.run_dir / 'predicted.sources', predicted)
        self.obj.diagnostics()
        fn = self.obj.diag_dir / 'source_regional_p30_ctl_inv1.csv'
        text = fn.read_text()
        self.assertIn('mean_abs_2000_2024', text)
        self.assertIn('seasonal_amplitude_2000_2024', text)
        self.assertTrue((self.obj.diag_dir / 'regional' / 'p30_ctl_inv1' /
                         'regional_ranking.png').is_file())

    def write_legacy_manifest(self, include_hashes=True, copy_inputs=False):
        from d4_inversion import Inversion
        a = self.obj
        if copy_inputs:
            for fname in ('sourcepriortdi.nc', 'obsrvCH4_9925.nc'):
                local = self.run_dir / fname
                data = local.read_bytes()
                local.unlink()
                local.write_bytes(data)
        content = ('# Immutable inputs for this inversion\n'
                   'case=INCA\n'
                   'icase=s057\n'
                   'prior=p30\n'
                   'unc=ctl\n'
                   'model=inv1\n'
                   'source_D1_prior=' + str(self.source) + '\n'
                   'source_D3_observations=' + str(self.obs) + '\n')
        if include_hashes:
            content += ('\n# SHA256 hashes of COPIED local inputs\n'
                        + a.sha256(self.source) + '  sourcepriortdi.nc\n'
                        + a.sha256(self.obs) + '  obsrvCH4_9925.nc\n')
        (self.run_dir / 'RUN_INPUTS.txt').write_text(content)

    def test_legacy_manifest_sha256_inventory_reused(self):
        self.obj.run()
        self.write_legacy_manifest()
        self.assertEqual(self.obj.run_inversion('p30', 'ctl', 'inv1',
                         self.source, self.obs, self.run_dir), 'REUSED')

    def test_legacy_manifest_sha256_inventory_detects_changed_input(self):
        self.obj.run()
        self.write_legacy_manifest()
        self.obs.write_bytes(self.obs.read_bytes() + b'change')
        with self.assertRaisesRegex(FileExistsError, 'D3 observations SHA256 mismatch'):
            self.obj.run_inversion('p30', 'ctl', 'inv1',
                                   self.source, self.obs, self.run_dir)
        self.assertTrue((self.run_dir / 'predicted.sources').is_file())

    def test_legacy_manifest_independent_copies_reused(self):
        self.obj.run()
        self.write_legacy_manifest(include_hashes=False, copy_inputs=True)
        self.assertEqual(self.obj.run_inversion('p30', 'ctl', 'inv1',
                         self.source, self.obs, self.run_dir), 'REUSED')

    def test_legacy_manifest_missing_hashes_symlinks_refused(self):
        self.obj.run()
        self.write_legacy_manifest(include_hashes=False, copy_inputs=False)
        # --- Explicitly simulate an older run that used symlinks.
        for name, source in (('sourcepriortdi.nc', self.source), ('obsrvCH4_9925.nc', self.obs)):
            local = self.run_dir / name
            local.unlink()
            local.symlink_to(source)
        with self.assertRaisesRegex(FileExistsError, 'contains symlinks'):
            self.obj.run_inversion('p30', 'ctl', 'inv1',
                                   self.source, self.obs, self.run_dir)

    def test_existing_manifest_mismatch_shows_specific_reason(self):
        self.obj.run()
        manifest = self.run_dir / 'RUN_INPUTS.txt'
        data = manifest.read_text().replace('source_D3_observations=' + str(self.obs),
                                            'source_D3_observations=/incorrect/obs.nc')
        manifest.write_text(data)
        with self.assertRaisesRegex(FileExistsError, 'D3 observations source path mismatch') as ctx:
            self.obj.run_inversion('p30', 'ctl', 'inv1',
                                   self.source, self.obs, self.run_dir)
        self.assertNotIn('problems: []', str(ctx.exception))

    def test_legacy_manifest_experiment_mismatch_refused(self):
        self.obj.run()
        self.write_legacy_manifest()
        manifest = self.run_dir / 'RUN_INPUTS.txt'
        manifest.write_text(manifest.read_text().replace('case=INCA', 'case=CYC'))
        with self.assertRaisesRegex(FileExistsError, 'case mismatch'):
            self.obj.run_inversion('p30', 'ctl', 'inv1',
                                   self.source, self.obs, self.run_dir)

    def test_no_parse_mixed_source(self):
        fn = self.root / 'sources.txt'
        fn.write_text('land_01\n10.4 25.0\n')
        self.assertIsNone(self.obj.read_source_coefficients(fn))



# ============================================================
# === Tests for genuine TDI block format read by D5
    def test_older_manifest_source_key_names_reused_with_independent_copies(self):
        self.obj.run()
        manifest = self.run_dir / 'RUN_INPUTS.txt'
        text = manifest.read_text()
        # --- Older stage labels are read as metadata, never written by new D4.
        text = text.replace('source_D1_prior=', 'source_F' + '5=')
        text = text.replace('source_D3_observations=', 'source_F' + '7=')
        text = text.replace('D1_prior_sha256=', 'old_prior_hash=')
        text = text.replace('D3_observations_sha256=', 'old_obs_hash=')
        manifest.write_text(text)
        reused = self.obj.run_inversion('p30', 'ctl', 'inv1',
                                        self.source, self.obs, self.run_dir)
        self.assertEqual(reused, 'REUSED')
        self.assertFalse(any(path.is_symlink() for path in self.run_dir.rglob('*')))

    def test_stage_names_in_manifest_and_report(self):
        self.obj.run()
        manifest = (self.run_dir / 'RUN_INPUTS.txt').read_text()
        self.assertIn('source_D1_prior=', manifest)
        self.assertIn('source_D3_observations=', manifest)
        self.assertIn('D1_prior_sha256=', manifest)
        self.assertIn('D3_observations_sha256=', manifest)
        self.assertNotRegex(manifest, r'\bF[0-9]+\b')
        self.obj.diagnostics()
        report = (self.obj.diag_dir / 'diagnostic_report.txt').read_text()
        self.assertIn('D1 prior uncertainty variants:', report)
        self.assertIn('Valid D3 station-months:', report)

    def test_diagnostic_deletes_only_old_stage_numbered_outputs(self):
        self.obj.run()
        self.obj.diag_dir.mkdir(parents=True, exist_ok=True)
        old_fig = self.obj.diag_dir / 'f5_old_prior_plot.png'
        old_table = self.obj.diag_dir / 'f7_old_obs.csv'
        own_file = self.obj.diag_dir / 'my_notes.txt'
        for item in (old_fig, old_table, own_file):
            item.write_text('test')
        self.obj.diagnostics()
        self.assertFalse(old_fig.exists())
        self.assertFalse(old_table.exists())
        self.assertTrue(own_file.exists())

    def test_overview_zero_solution_is_diagnostic_not_result(self):
        self.obj.run()
        np.savetxt(self.run_dir / 'predicted.sources', np.zeros((324, 54)))
        with self.assertRaisesRegex(RuntimeError, 'NUMERICAL QUALITY FAILED'):
            self.obj.diagnostics()
        overview = self.obj.diag_dir / 'inversion_overview_INCA_p30_ctl_inv1.png'
        self.assertTrue(overview.is_file())
        self.assertGreater(overview.stat().st_size, 1000)
        report = (self.obj.diag_dir / 'diagnostic_report.txt').read_text()
        self.assertIn('NUMERICAL CRITICAL: 1', report)
        self.assertIn('Nonzero coefficients: 0/17496', report)


class TestD4BlockSource(unittest.TestCase):

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write_blocks(self, invalid_time=False, missing_month=False):
        file = self.root / 'predicted.sources'
        with file.open('w') as f:
            for region in range(1, 55):
                f.write('land_' + str(region).zfill(2) + '\n')
                f.write('time  correction\n')
                for k in range(324):
                    if missing_month and region == 54 and k == 323:
                        continue
                    t = 1999 + k // 12 + ((k % 12) + 0.5) / 12
                    if invalid_time and region == 2 and k == 10:
                        t += 1.0
                    v = float(region) * 0.01 + k * 0.001
                    f.write(f'{t:16.7f} {v:16.8f}\n')
        return file

    def test_actual_tdi_54_region_blocks(self):
        fn = self.write_blocks()
        got = Inversion.read_source_coefficients(fn)
        self.assertIsNotNone(got)
        self.assertEqual(got.shape, (324, 54))
        self.assertAlmostEqual(got[0, 0], 0.01)
        self.assertAlmostEqual(got[323, 53], 0.54 + 0.323)

    def test_tdi_block_wrong_time_rejected(self):
        fn = self.write_blocks(invalid_time=True)
        self.assertIsNone(Inversion.read_source_coefficients(fn))

    def test_tdi_block_missing_month_rejected(self):
        fn = self.write_blocks(missing_month=True)
        self.assertIsNone(Inversion.read_source_coefficients(fn))

    def test_empty_month_aggregation_no_warning(self):
        import warnings
        from unittest.mock import patch
        import numpy as np
        previous = os.environ.get('GCP_CASE')
        os.environ['GCP_CASE'] = 'INCA'
        try:
            inv = Inversion(case='INCA')
        finally:
            if previous is None:
                os.environ.pop('GCP_CASE', None)
            else:
                os.environ['GCP_CASE'] = previous
        inv.diag_dir = self.root / 'diagnostics'
        data = np.full((324, 5), np.nan)
        sig = np.full((324, 5), np.nan)
        data[1, 1] = 3.0
        sig[1, 1] = 9.0
        with patch.object(inv, 'get_observation_data', return_value=(data, sig)):
            with warnings.catch_warnings(record=True) as records:
                warnings.simplefilter('always', RuntimeWarning)
                out = inv.plot_observation_diagnostics(self.root / 'dummy.nc')
        self.assertEqual(out['n_empty_months'], 323)
        self.assertEqual(len([w for w in records if issubclass(w.category, RuntimeWarning)]), 0)
        csvfile = (inv.diag_dir / 'observation_monthly.csv').read_text()
        self.assertIn('1999-01,0,nan,nan', csvfile)


if __name__ == '__main__':
    unittest.main(verbosity=2)
