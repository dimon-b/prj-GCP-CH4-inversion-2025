# -*- coding: utf-8 -*-
"""Synthetic unit and end-to-end checks; no production data changed."""

import importlib.util
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import xarray as xr


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('d0_inca_oh', ROOT / 'd0_inca_oh.py')
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def check_vertical():
    # --- Source levels in deliberately descending pressure order
    ps = np.full((2, 3), 100000.0)
    p = np.broadcast_to(np.array([90000.0, 50000.0, 10000.0])[:, None, None], (3, 2, 3)).copy()
    v = 1e6 + 10.0 * p
    v[:, 0, 0] = [1.9e6, 1.5e6, 1.1e6]
    sigma = np.array([0.95, 0.5, 0.3, 0.1, 0.05])
    a, outside = mod.interp_vertical(v, p, ps, sigma, max_outside_pct=90)
    assert np.isfinite(a).all()
    assert a.shape == (5, 2, 3)
    assert abs(a[0, 1, 1] - 1.9e6) < 1
    assert abs(a[1, 1, 1] - 1.5e6) < 1
    assert a[2, 1, 1] > 1.1e6 and a[2, 1, 1] < 1.5e6
    assert abs(a[3, 1, 1] - 1.1e6) < 1
    assert abs(a[4, 1, 1] - 1.1e6) < 1
    assert outside > 0

    # --- One source upper level NaN -> no missing output
    v[2, :, :] = np.nan
    a, _ = mod.interp_vertical(v, p, ps, sigma, max_outside_pct=90)
    assert np.isfinite(a).all()
    assert a.min() >= 0

    # --- Cannot silently invent surface pressure
    ds = xr.Dataset({'OH': (('time_counter', 'presnivs', 'lat', 'lon'), np.ones((1, 3, 2, 3)), {'units': 'molecules/m3'}),
                     'pmid': (('time_counter', 'presnivs', 'lat', 'lon'), np.ones((1, 3, 2, 3)), {'units': 'Pa'})},
                    coords={'time_counter': pd.date_range('2000-01-01', periods=1),
                            'presnivs': [1, 2, 3], 'lat': [-20, 20], 'lon': [0, 120, 240]})
    try:
        mod.get_source_metadata(ds, SimpleNamespace(ps_var=None, assume_oh_unit=None, assume_pressure_unit=None))
    except ValueError as e:
        assert 'requires pmid' in str(e)
    else:
        raise AssertionError('Missing ps was accepted')


def check_periodic():
    # --- Longitude 0 must use neighbors from both sides
    lat = np.array([-30., 30.])
    lon = np.array([45., 135., 225., 315.])
    ds = xr.Dataset({'OH': (('z', 'y', 'x'), np.ones((3, 2, 4))),
                     'pmid': (('z', 'y', 'x'), np.ones((3, 2, 4))),
                     'ps': (('y', 'x'), np.ones((2, 4)))},
                    coords={'z': [1, 2, 3], 'y': lat, 'x': lon})
    cfg = dict(z='z', y='y', x='x')
    out = mod.horizontal_regrid(ds, cfg, lat, np.array([0., 90., 180., 270.]))
    assert np.isfinite(out['OH']).all()
    assert out.sizes['x'] == 4


def check_monthly_conversion():
    # --- Full daily coverage in leap year tests daily -> monthly interpolation pipeline
    dates = pd.date_range('2000-01-01', '2000-12-31', freq='D')
    lat = np.array([-30., 30.])
    lon = np.array([45., 135., 225., 315.])
    pressure = np.array([90000., 50000., 10000.])
    p = np.broadcast_to(pressure[None, :, None, None], (len(dates), 3, 2, 4)).copy()
    v = 1e6 + 10.0 * p
    ps = np.full((len(dates), 2, 4), 100000.0)
    ds = xr.Dataset({
        'OH': (('time_counter', 'presnivs', 'lat', 'lon'), v, {'units': 'molecules/m3'}),
        'pmid': (('time_counter', 'presnivs', 'lat', 'lon'), p, {'units': 'Pa'}),
        'ps': (('time_counter', 'lat', 'lon'), ps, {'units': 'Pa'}),
    }, coords={'time_counter': dates, 'presnivs': [1, 2, 3], 'lat': lat, 'lon': lon})
    args = SimpleNamespace(ps_var=None, assume_oh_unit=None, assume_pressure_unit=None,
                           min_days=20, max_outside_pct=1.0, allow_incomplete_months=False)
    cfg = mod.get_source_metadata(ds, args)
    out = mod.build_months(ds, cfg, 2000, lat, np.array([0., 90., 180., 270.]),
                           np.array([0.9, 0.5, 0.1]), args)
    assert len(out) == 12
    assert all(x.shape == (3, 2, 4) for x in out)
    assert np.isfinite(np.stack(out)).all()
    assert abs(float(out[0][0, 0, 0]) - 1.9) < 1e-5


def check_end_to_end():
    # --- Requires netCDF4 to be available; the server already has this
    try:
        import netCDF4  # noqa: F401
    except ImportError:
        print('NetCDF4 output test: SKIPPED (package not installed in this test environment)')
        return

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        src = tmp / 'daily'
        src.mkdir()
        ref = tmp / 'actm.nc'
        out = tmp / 'v2.nc'
        lat = [-30., 30.]
        lon = [45., 135., 225., 315.]
        sigma = [0.9, 0.5, 0.1]
        xr.Dataset(coords={'level': ('level', sigma), 'y': ('y', lat), 'x': ('x', lon)}).to_netcdf(ref)
        year = 2000
        dates = pd.date_range(f'{year}-01-01', f'{year}-12-31', freq='D')
        psrc = np.array([90000., 50000., 10000.])
        pmid = np.broadcast_to(psrc[None, :, None, None], (len(dates), 3, 2, 4)).copy()
        oh = 1e6 + 10.0 * pmid
        ps = np.full((len(dates), 2, 4), 100000.)
        ds = xr.Dataset(
            data_vars={
                'OH': (('time_counter', 'presnivs', 'lat', 'lon'), oh, {'units': 'molecules/m3'}),
                'pmid': (('time_counter', 'presnivs', 'lat', 'lon'), pmid, {'units': 'Pa'}),
                'ps': (('time_counter', 'lat', 'lon'), ps, {'units': 'Pa'}),
            }, coords={'time_counter': dates, 'presnivs': [1, 2, 3], 'lat': lat, 'lon': lon})
        ds.to_netcdf(src / 'INCA_DM_2000_OH_scaled.nc')
        args = SimpleNamespace(input_dir=str(src), reference=str(ref), output=str(out),
                               start_year=2000, end_year=2000, min_days=20,
                               max_outside_pct=1.0, extend_edges=False,
                               ps_var=None, assume_oh_unit=None, assume_pressure_unit=None)
        mod.run(args)
        with xr.open_dataset(out) as new:
            assert new.sizes == {'time': 12, 'level': 3, 'y': 2, 'x': 4}
            assert np.isfinite(new.OH.values).all()
            assert abs(float(new.OH.isel(time=0, level=0, y=0, x=0)) - 1.9) < 1e-5
        try:
            mod.run(args)
        except FileExistsError:
            pass
        else:
            raise AssertionError('Overwrite guard failed')



def check_nominal_mode():
    # --- OH-only LMDZ-style file: presnivs is pseudo-pressure at reference Ps
    dates = pd.date_range('2000-01-01', '2000-12-31', freq='D')
    lat = np.array([-30., 30.])
    lon = np.array([45., 135., 225., 315.])
    pvals = np.array([90000., 50000., 10000.])
    input_oh = np.broadcast_to((1e6 + 10 * pvals)[None, :, None, None],
                               (len(dates), 3, 2, 4)).copy()
    ds = xr.Dataset({'OH': (('time_counter', 'presnivs', 'lat', 'lon'),
                            input_oh, {'units': 'molec. m-3'})},
                    coords={'time_counter': dates, 'presnivs': ('presnivs', pvals, {'units': 'Pa'}),
                            'lat': lat, 'lon': lon})
    args = SimpleNamespace(vertical_mode='nominal', assume_oh_unit=None,
                           assume_pressure_unit=None, nominal_reference_pressure_pa=100000.,
                           min_days=20, max_outside_pct=1.0, allow_incomplete_months=False)
    cfg = mod.get_source_metadata(ds, args)
    assert cfg['mode'] == 'nominal'
    assert cfg['oh_factor'] == 1e-6
    out = mod.build_months(ds, cfg, 2000, lat, np.array([0., 90., 180., 270.]),
                           np.array([0.9, 0.5, 0.1]), args)
    assert np.isfinite(np.stack(out)).all()
    assert abs(float(out[0][0, 0, 0]) - 1.9) < 1e-5

    # --- Real units are compulsory, no silent default
    ds['presnivs'].attrs.clear()
    try:
        mod.get_source_metadata(ds, args)
    except ValueError as e:
        assert 'Pressure units not confirmed' in str(e)
    else:
        raise AssertionError('Nominal mode silently assumed presnivs units')


def check_nominal_end_to_end():
    try:
        import netCDF4
    except ImportError:
        print('Nominal end-to-end: SKIPPED (netCDF4 absent)')
        return
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        src = tmp / 'daily'
        src.mkdir()
        ref = tmp / 'actm.nc'
        out = tmp / 'nominal.nc'
        lat = [-30., 30.]
        lon = [45., 135., 225., 315.]
        sigma = [0.9, 0.5, 0.1]
        xr.Dataset(coords={'level': ('level', sigma), 'y': ('y', lat),
                           'x': ('x', [0., 90., 180., 270.])}).to_netcdf(ref)
        dates = pd.date_range('2000-01-01', '2000-12-31', freq='D')
        pvals = np.array([90000., 50000., 10000.])
        oh = np.broadcast_to((1e6 + 10*pvals)[None, :, None, None],
                             (len(dates), 3, 2, 4)).copy()
        ds = xr.Dataset({'OH': (('time_counter', 'presnivs', 'lat', 'lon'),
                                oh, {'units': 'molec. m-3'})},
                        coords={'time_counter': dates, 'presnivs': ('presnivs', pvals, {'units': 'Pa'}),
                                'lat': lat, 'lon': lon})
        ds.to_netcdf(src / 'INCA_DM_2000_OH_scaled.nc')
        args = SimpleNamespace(input_dir=str(src), reference=str(ref), output=str(out),
                               start_year=2000, end_year=2000, min_days=20,
                               max_outside_pct=1., extend_edges=False, vertical_mode='nominal',
                               nominal_reference_pressure_pa=100000., ps_var=None,
                               assume_oh_unit=None, assume_pressure_unit=None)
        assert mod.inspect_daily(src / 'INCA_DM_2000_OH_scaled.nc', args)
        mod.run(args)
        with xr.open_dataset(out) as result:
            assert result.attrs['vertical_mapping'] == 'nominal'
            assert result.attrs['nominal_reference_pressure_pa'] == 100000.
            assert result.sizes['time'] == 12
            assert np.isfinite(result.OH.values).all()


def check_daily_coverage():
    # --- Leap year: all 366 daily timestamps required, including February 29
    dates = pd.date_range('2000-01-01', '2000-12-31', freq='D')
    pres = np.array([90000., 50000., 10000.])
    lat = np.array([-30., 30.])
    lon = np.array([45., 135., 225., 315.])
    oh = np.broadcast_to((1e6 + 10*pres)[None, :, None, None],
                         (len(dates), 3, 2, 4)).copy()
    ds = xr.Dataset({'OH': (('time_counter', 'presnivs', 'lat', 'lon'),
                            oh, {'units': 'molec. m-3'})},
                    coords={'time_counter': dates,
                            'presnivs': ('presnivs', pres, {'units': 'Pa'}),
                            'lat': lat, 'lon': lon})
    args = SimpleNamespace(vertical_mode='nominal', assume_oh_unit=None,
                           assume_pressure_unit=None, nominal_reference_pressure_pa=100000.,
                           min_days=20, max_outside_pct=1., allow_incomplete_months=False)
    cfg = mod.get_source_metadata(ds, args)
    result = mod.build_months(ds, cfg, 2000, lat, np.array([0., 90., 180., 270.]),
                              np.array([0.9, 0.5, 0.1]), args)
    assert len(result) == 12
    assert np.isfinite(np.stack(result)).all()

    # --- Missing a date must block a full-period production conversion
    gap = ds.drop_sel(time_counter=pd.Timestamp('2000-02-29'))
    try:
        mod.build_months(gap, cfg, 2000, lat, np.array([0., 90., 180., 270.]),
                         np.array([0.9, 0.5, 0.1]), args)
    except ValueError as err:
        assert '366 daily records' in str(err) or 'missing daily records' in str(err)
    else:
        raise AssertionError('Missing February 29 was ignored')

    # --- Duplicate a day must block too
    doubled = xr.concat([ds, ds.isel(time_counter=[0])], dim='time_counter')
    try:
        mod.build_months(doubled, cfg, 2000, lat, np.array([0., 90., 180., 270.]),
                         np.array([0.9, 0.5, 0.1]), args)
    except ValueError as err:
        assert 'duplicate daily dates' in str(err)
    else:
        raise AssertionError('Duplicate day was ignored')


def check_full_period():
    # --- 25 source years + copied 1999/2025 edges = 324 months.
    # --- Every synthetic source year includes all days; no bypass of production validation.
    try:
        import netCDF4
    except ImportError:
        print('Full-period 324-month test: SKIPPED (netCDF4 absent)')
        return
    with tempfile.TemporaryDirectory() as dirname:
        tmp = Path(dirname)
        src = tmp / 'daily'
        src.mkdir()
        ref = tmp / 'actm.nc'
        out = tmp / '324months.nc'
        lat = [-30., 30.]
        lon = [45., 135., 225., 315.]
        sigma = [0.9, 0.5, 0.1]
        xr.Dataset(coords={'level': ('level', sigma), 'y': ('y', lat),
                           'x': ('x', [0., 90., 180., 270.])}).to_netcdf(ref)
        pressure = [90000., 50000., 10000.]
        for year in range(2000, 2025):
            dates = pd.date_range(f'{year}-01-01', f'{year}-12-31', freq='D')
            value = np.broadcast_to(np.array(pressure)[None, :, None, None] * 10.0 +
                                    float(year), (len(dates), 3, 2, 4)).copy()
            ds = xr.Dataset({'OH': (('time_counter', 'presnivs', 'lat', 'lon'), value,
                                    {'units': 'molec. m-3'})},
                            coords={'time_counter': dates,
                                    'presnivs': ('presnivs', pressure, {'units': 'Pa'}),
                                    'lat': lat, 'lon': lon})
            ds.to_netcdf(src/f'INCA_DM_{year}_OH_scaled.nc')
        args = SimpleNamespace(input_dir=str(src), reference=str(ref), output=str(out),
                               start_year=2000, end_year=2024, vertical_mode='nominal',
                               nominal_reference_pressure_pa=100000., ps_var=None,
                               assume_oh_unit=None, assume_pressure_unit=None,
                               min_days=20, allow_incomplete_months=False,
                               max_outside_pct=1., extend_edges=True)
        mod.run(args)
        with xr.open_dataset(out) as result:
            assert result.sizes['time'] == 324
            times = pd.DatetimeIndex(result.time.values)
            assert (times[0].year, times[-1].year) == (1999, 2025)
            assert np.array_equal(result.OH.isel(time=0).values, result.OH.isel(time=12).values)
            assert np.array_equal(result.OH.isel(time=-1).values, result.OH.isel(time=-13).values)
            assert np.isfinite(result.OH.values).all()
        print('INCA OH 324-month edge-extension test: PASS')


def main():
    check_vertical()
    check_periodic()
    check_monthly_conversion()
    check_end_to_end()
    check_nominal_mode()
    check_nominal_end_to_end()
    check_daily_coverage()
    check_full_period()
    print('INCA OH v2 physical and nominal synthetic checks: PASS')


if __name__ == '__main__':
    main()
