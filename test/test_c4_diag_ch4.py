# -*- coding: utf-8 -*-
"""Synthetic regression test for read-only F6 / F7 diagnostics."""
import os
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import c4_diag_ch4 as d


def main():
    with tempfile.TemporaryDirectory() as temp:
        base = Path(temp)
        inv = base / 'inv_dir'
        inv.mkdir()
        site_path = inv / 'obsrvCH4_60sites.txt'
        site_path.write_text('header one\nheader two\n' + ''.join(
            f'{n:5d} {9999 if n <= 2 else 0:5d} ch4_test{n:03d}_surface-flask_1 50.0 10.0 100.0\n'
            for n in range(1, d.N_SITE + 1)))

        inp_path = inv / 'obsrv_actm_cyc_ch4dy.grd'
        out_path = inv / 'ch4obs_actm_cyc_fits.bin'
        inp = np.memmap(inp_path, dtype='<f4', mode='w+', shape=(d.N_DAY, 2, d.N_SITE))
        fits = np.memmap(out_path, dtype='<f4', mode='w+', shape=(d.N_DAY, 2, 3, d.N_SITE))
        inp[:] = d.FILL
        fits[:] = d.FILL
        inp[:, 0, 0] = 1850.0
        inp[:, 1, 0] = 1840.0
        fits[:, 0, 0, 0] = 1850.0
        fits[:, 1, 0, 0] = 1840.0
        fits[:, 0, 1, 0] = 1849.0
        fits[:, 1, 1, 0] = 1839.0
        fits[:, 0, 2, 0] = 1830.0
        fits[:, 1, 2, 0] = 1820.0
        inp.flush(); fits.flush()
        del inp, fits

        days, months = d.get_dates()
        a, b = d.get_f6_files(inp_path, out_path)
        monthly = d.get_f6_monthly(b, days, months)
        expected = d.get_expected_f7(monthly, 1, months)
        assert np.allclose(expected[:, 0, 0], 10.0)
        assert np.allclose(expected[:, 0, 1], 6.0)
        assert np.all(expected[:, 1:, 0] == d.FILL)

        sites = d.get_sites(site_path)
        nc_path = inv / 'obsrvCH4_cyc_9925ctlinv1.nc'
        with netcdf_file(str(nc_path), 'w') as nc:
            for name, size in [('station_count', d.N_SITE), ('tlen', d.N_MONTH), ('data_type', 2), ('namelen', 45)]:
                nc.createDimension(name, size)
            numbers = nc.createVariable('station_count', 'i2', ('station_count',))
            names = nc.createVariable('site_name', 'c', ('station_count', 'namelen'))
            time = nc.createVariable('tlen', 'f4', ('tlen',))
            data = nc.createVariable('tracer_data', 'f4', ('tlen', 'station_count', 'data_type'))
            numbers[:] = sites.station_id.to_numpy(dtype='i2')
            raw_names = np.full((d.N_SITE, 45), b' ', dtype='S1')
            for i, name in enumerate(sites.site):
                raw_names[i, :len(name)] = np.frombuffer(name.encode('ascii'), dtype='S1')
            names[:] = raw_names
            time[:] = np.asarray(1999 + (np.arange(d.N_MONTH) + 0.5) / 12.0, dtype='f4')
            data[:] = expected.astype('f4')
        nc = d.get_f7_file(nc_path)
        stats = d.get_f7_stats(nc, sites)
        assert stats.iloc[0]['n_valid'] == d.N_MONTH
        assert stats.iloc[1]['n_valid'] == 0
        cmp, _ = d.get_comparison(nc, expected, sites, months)
        assert cmp['n_total_bad_months'] == 0, cmp

        # --- The actual inversion interval comes from control, not output years_nc.
        control = base / 'control'
        control.write_text('&control\n firstdat = 2000\n lastdat = 2024\n/\n')
        period_setup = SimpleNamespace(inv_control_template=str(control), inv_year_s=1999, inv_year_e=2025)
        first, last, source = d.get_inversion_period(period_setup)
        assert (first, last) == (2000, 2024)
        assert source.startswith('control:')
        invcheck = d.get_inversion_missing(nc, sites, months, first, last)
        assert invcheck['n_requested'] == 2
        assert invcheck['n_months'] == 300
        assert invcheck['n_invalid'] == 300  # station 2 has no data at all
        assert invcheck['by_id'][1]['missing'] == 0
        assert invcheck['by_id'][2]['fill'] == 300

        # --- Duplicated names: legacy read_globalview.f selects the LAST occurrence.
        duplicate_sites = sites.copy()
        duplicate_sites.loc[1, 'site'] = duplicate_sites.loc[0, 'site']
        duplicate_nc = dict(nc)
        duplicate_nc['names'] = list(nc['names'])
        duplicate_nc['names'][1] = duplicate_nc['names'][0]
        duplicated = d.get_duplicate_active_sites(duplicate_sites, duplicate_nc)
        assert len(duplicated) == 1
        assert duplicated[0]['ids'] == [1, 2]
        assert duplicated[0]['rows_identical'] is False
        duplicate_check = d.get_inversion_missing(duplicate_nc, duplicate_sites, months, first, last)
        assert duplicate_check['n_invalid'] == 600  # both active IDs read the missing second record
        assert duplicate_check['by_id'][1]['netcdf_station_id_used'] == 2

        # --- Numerically identical ACTIVE series with distinct station names.
        same_nc = dict(nc)
        same_data = np.array(nc['data'], copy=True)
        same_data[:, 1, :] = same_data[:, 0, :]
        same_nc['data'] = same_data
        pairs = d.get_identical_active_series(sites, same_nc)
        assert len(pairs) == 1
        assert pairs[0]['id_a'] == 1 and pairs[0]['id_b'] == 2
        assert pairs[0]['same_uncertainty'] is True
        assert pairs[0]['same_name'] is False
        same_data[10, 1, 1] += 0.01
        pairs = d.get_identical_active_series(sites, same_nc)
        assert len(pairs) == 1 and pairs[0]['same_uncertainty'] is False
        same_data[10, 1, 0] += 0.01
        assert d.get_identical_active_series(sites, same_nc) == []

        # --- Phase profile and margin: ACTM one month later, without assuming confidence.
        x = np.sin(np.arange(12) * np.pi / 6.0)
        lag, corr = d.get_best_lag(x, np.roll(x, 1))
        profile = d.get_lag_correlations(x, np.roll(x, 1))
        assert lag == 1 and np.isclose(corr, 1.0)
        assert len(profile) == 7
        assert profile[1] > profile[0]

        # --- F7 monthly missing raw OBS doubles sigma while fitted data stay valid.
        monthly_no_raw = dict(monthly)
        monthly_no_raw['obs_raw_mean'] = monthly['obs_raw_mean'].copy()
        monthly_no_raw['obs_raw_mean'][0, 0] = np.nan
        expected_no_raw = d.get_expected_f7(monthly_no_raw, 1, months)
        assert np.isclose(expected_no_raw[0, 0, 1], 12.0)
        assert np.isclose(expected_no_raw[1, 0, 1], 6.0)

        setup = SimpleNamespace(sites_f=str(site_path), inv_wrk_dir=str(base/'results2025'/'CYC'),
                                f6_input_file=str(inp_path), f6_output_file=str(out_path),
                                f7_output_nc_prefix=str(inv / 'obsrvCH4_cyc_9925'))
        result = d.run_diagnostics('CYC', setup, plots=True, stations=[1, 42])
        assert result['station_ids'] == [1, 2], result
        assert d.get_plot_indices(sites, [1, 42]) == [0, 1]
        assert d.get_plot_indices(sites) == [0, 1]
        assert d.get_plot_indices(sites, [3, 42]) == []
        try:
            d.get_plot_indices(sites, [42, 1])
            raise AssertionError('Reversed station range not rejected')
        except ValueError:
            pass
        diagdir = base/'results2025'/'CYC'/'diagnostics'
        assert (diagdir/'ch4_station_stats.txt').is_file()
        assert not list(diagdir.glob('*.csv'))
        report = (diagdir/'ch4_station_stats.txt').read_text()
        for heading in ('FirstF7', 'LastF7', 'NoRaw', 'InvBad', 'SEASONAL LAG CORRELATIONS',
                        'DUPLICATE ACTIVE STATION NAMES', 'IDENTICAL ACTIVE F7 SERIES', 'INVERSION INPUT VALIDITY',
                        'literal -999'):
            assert heading in report, heading
        # --- Verify figures are generated for all selected ACTIVE sites (1, 2)
        figures = sorted(diagdir.glob('station_*.png'))
        assert len(figures) == 2, figures
        assert all(img.stat().st_size > 10000 for img in figures)
        print('F6/F7 diagnostics synthetic test: PASS')


if __name__ == '__main__':
    main()
