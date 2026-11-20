# -*- coding: utf-8 -*-
"""Read-only integration test of the three REAL D1 Fortran outputs.
Run after: python _main.py --case INCA --steps d1_prior
No compiler or synthetic stub is used. Does not modify any production output.
"""
import os
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from d1_prior_fluxes import PriorFluxes


class FortranIntegration(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        os.environ['GCP_CASE'] = 'INCA'
        os.environ['CH4_CASE'] = 'INCA'
        cls.job = PriorFluxes(case='INCA')

    def test_01_fortran_sources_and_namelists(self):
        names = {'regrid': 'd1_regrid_ch4_t42',
                 'regions': 'd1_aggregate_ch4_regions',
                 'uncertainty': 'd1_write_prior_uncertainty'}
        for stage, group in names.items():
            binary = Path(getattr(self.job, 'exe_d1_' + stage))
            source = Path(self.job.fort_dir) / (binary.name + '.f90')
            self.assertTrue(source.is_file(), 'Missing Fortran source: ' + str(source))
            self.assertIn('namelist /' + group + '/', source.read_text().lower())
            # --- The binary is only built if the corresponding output needs regeneration.

    def test_02_prior_netCDF(self):
        for name in (self.job.flx_apr_nc, self.job.flx_apr_nc_full):
            path = Path(self.job.flx_apr_dir) / name
            self.assertTrue(self.job.is_fortran_nc(path), 'D1 regrid needs NetCDF3 64-bit: ' + str(path))
            self.job.check_prior_nc(path)

    def test_03_regridded_gt3_and_conservation(self):
        gt3 = Path(self.job.d1_gt3_file)
        self.job.check_gt3(gt3)
        self.job.check_flux_conservation(Path(self.job.d1_regrid_input_dir) / self.job.d1_regrid_input_nc, gt3)

    def test_04_regions_and_provenance(self):
        import hashlib
        reg = Path(self.job.d1_region_table_file)
        self.job.check_regional_table(reg)
        marker = Path(str(reg) + '.d1_1998-2026.sha256')
        self.assertTrue(marker.is_file(), 'D1 regional provenance marker missing')
        with Path(self.job.d1_gt3_file).open('rb') as f:
            digest = hashlib.file_digest(f, 'sha256').hexdigest()
        self.assertEqual(marker.read_text().strip(), digest)

    def test_05_uncertainty(self):
        for pc in (30, 50, 99):
            p = Path(str(self.job.d1_uncertainty_output_prefix) + str(pc) + '.nc')
            self.job.check_prior_uncertainty_nc(p)
            self.assertGreaterEqual(p.stat().st_mtime, Path(self.job.d1_region_table_file).stat().st_mtime)


if __name__ == '__main__':
    unittest.main(verbosity=2)
