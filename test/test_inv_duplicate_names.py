# -*- coding: utf-8 -*-
"""Test that NetCDF inspection permits duplicate site names, not missing data."""
import sys
import tempfile
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'inversion'))
from check_inv_inputs import inspect_observations


def main():
    with tempfile.TemporaryDirectory() as temp:
        folder = Path(temp)
        fname = folder / 'obsrvCH4_cyc_9925ctlinv1.nc'
        sitelist = folder / 'obsrvCH4_60sites.txt'
        sitelist.write_text('header one\nheader two\n' + ''.join(
            f'{i:3d} {9999 if i <= 2 else 0:4d} station_{1 if i in (1, 2) else i:03d} 10 20 30\n'
            for i in range(1, 246)))
        with netcdf_file(str(fname), 'w') as nc:
            for name, size in (('station_count', 245), ('tlen', 324), ('data_type', 2), ('namelen', 45)):
                nc.createDimension(name, size)
            nums = nc.createVariable('station_count', 'i2', ('station_count',))
            names = nc.createVariable('site_name', 'c', ('station_count', 'namelen'))
            years = nc.createVariable('tlen', 'f4', ('tlen',))
            data = nc.createVariable('tracer_data', 'f4', ('tlen', 'station_count', 'data_type'))
            nums[:] = np.arange(1, 246, dtype='i2')
            raw = np.full((245, 45), b' ', dtype='S1')
            for i in range(245):
                name = 'station_001' if i <= 1 else f'station_{i+1:03d}'
                raw[i, :len(name)] = np.frombuffer(name.encode(), dtype='S1')
            names[:] = raw
            years[:] = np.asarray(1999 + (np.arange(324) + 0.5)/12., dtype='f4')
            data[:] = np.array([1., 10.], dtype='f4')
        info = inspect_observations(fname, sitelist, strict=True)
        assert info['n_active_stations'] == 2, info
        assert info['n_invalid_active'] == 0, info
        with netcdf_file(str(fname), 'a') as nc:
            dat = np.array(nc.variables['tracer_data'][:], copy=True)
            dat[1, 1, :] = -999.
            nc.variables['tracer_data'][:] = dat
        try:
            inspect_observations(fname, sitelist, strict=True)
            raise AssertionError('Missing active month was not rejected')
        except ValueError as error:
            assert 'missing/invalid' in str(error), error
        print('Inversion input duplicate-name and missing-value guard: PASS')


if __name__ == '__main__':
    main()
