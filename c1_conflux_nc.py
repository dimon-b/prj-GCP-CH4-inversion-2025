# -*- coding: utf-8 -*-
"""
Finalized:  26/xx/xx
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib import gridspec
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import regionmask
import _set_case


class CnvFluxNc(_set_case.SetCase):
    def __init__(self):
        super().__init__()

        print('Read Input Fluxes from dir: ')
        print('\t\t', self.flx_inp_dir)

        # ---
        self.get_inp_flx()

    # --- get_inp_flx
    def get_inp_flx(self):

        # --- get_inp_nc
        def get_inp_nc():
            path = self.flx_inp_dir + self.flx_inp_nc
            ds_ = xr.open_dataset(path, engine="h5netcdf", decode_times=False)
            return ds_

        # --- manage_ds
        def manage_ds(ds_):
            # - get lengths
            ntime = ds_.sizes["time"]
            nseason = ds_.sizes["time_climato"]

            # - make seasonal variables match the full time axis
            expanded = {}
            for v in [v for v in ds_.data_vars if "time_climato" in ds_[v].dims]:
                # tile the seasonal data to match the time dimension
                arr = ds_[v]
                reps = ntime // nseason
                expanded[v] = xr.DataArray(
                    np.tile(arr.values, (reps, 1, 1)),  # repeat along time axis
                    dims=("time", "lat", "lon"),
                    coords={"time": ds_.time, "lat": ds_.lat, "lon": ds_.lon},
                    attrs=arr.attrs)

            # - keep the original "time" variables as they are
            time_vars = {v: ds_[v] for v in ds_.data_vars if "time" in ds_[v].dims}

            # - join
            ds_expanded = xr.Dataset({**time_vars, **expanded})

            # - Select sector fluxes only
            flux_vars = [v for v in ds_expanded.data_vars if
                         v.startswith('flux_ch4_') and v not in ['flux_ch4_prior', 'flux_ch4_prior_soil0']]
            for v in flux_vars:
                print(v, 'NaN cells:', int(ds_expanded[v].isnull().sum().item()))

            # - Treat masked soil cells as zero uptake before summing
            ds_expanded['flux_ch4_soils'] = ds_expanded['flux_ch4_soils'].fillna(0.0)

            # - Add total and total without soil
            ds_expanded['flux_ch4_prior'] = sum(ds_expanded[v] for v in flux_vars)
            ds_expanded['flux_ch4_prior_soil0'] = ds_expanded['flux_ch4_prior'] - ds_expanded['flux_ch4_soils']

            # - Make sure time is datetime
            if not np.issubdtype(ds_expanded.time.dtype, np.datetime64):
                ds_expanded['time'] = xr.decode_cf(ds_expanded[['time']]).time

            return ds_expanded

        # --- check_total
        def check_total(ds):

            # - Grid-cell area [m2]
            dlat = np.deg2rad(ds.lat[1] - ds.lat[0])
            dlon = np.deg2rad(ds.lon[1] - ds.lon[0])
            lat_rad = np.deg2rad(ds.lat)
            # cell_area_2d = (self.R ** 2 * dlat * dlon * np.cos(lat_rad))
            cell_area_2d = xr.DataArray(self.R ** 2 * dlat * dlon * np.cos(lat_rad),
                                        coords={"lat": ds.lat},
                                        dims=["lat"])
            cell_area = cell_area_2d.broadcast_like(ds["flux_ch4_wetlands"])

            # - sanity check: total Earth surface area
            earth_area = cell_area.isel(time=0).sum(dim=("lat", "lon"))
            print("")
            print(f"Check total Earth area = {earth_area.item():.2e} m2")
            print("Expected ≈ 5.10e14 m2")
            print("")

            # - Seconds per month (time-aware, xarray-safe)
            time = pd.to_datetime(ds.time.values)
            seconds_per_month = xr.DataArray(time.days_in_month * 24 * 3600,
                                             coords={"time": ds.time},
                                             dims=("time",),
                                             name="seconds_per_month",
                                             )

            # - Flux variables to integrate
            flux_vars = [v for v in ds.data_vars if v.startswith("flux_ch4_")]

            # - Integration
            records = []
            for var in flux_vars:
                # kg m-2 s-1 -> kg per month per grid cell
                monthly_flux = ds[var] * cell_area * seconds_per_month
                # global monthly sum [kg]
                global_monthly = monthly_flux.sum(dim=("lat", "lon"))
                # annual totals [Tg yr-1]
                global_annual = global_monthly.groupby("time.year").sum(dim="time") / 1e9
                for year, value in zip(global_annual.year.values, global_annual.values):
                    records.append(
                        {"year": int(year), "variable": var.replace("flux_ch4_", ""), "Tg_CH4_yr": float(value), })

            # - Output table
            df = pd.DataFrame(records)
            table = df.pivot(index="year", columns="variable", values="Tg_CH4_yr", )

            # - Put totals at the end
            total_cols = [c for c in table.columns if c.lower().startswith("prior")]
            other_cols = [c for c in table.columns if c not in total_cols]
            table = table[other_cols + total_cols]

            print(table.round(4))

        # - write 'flux_ch4_total' to NetCDF
        def write_nc_0(ds_i):
            ds_ = ds_i.copy()

            # - Copy full year 2000 → 1999
            year_1999 = ds_.sel(time=slice("2000-01", "2000-12"))
            year_1999 = year_1999.assign_coords(time=pd.date_range("1999-01-01", "1999-12-01", freq="MS"))

            # - Copy full year 2024 → 2025
            year_2025 = ds_.sel(time=slice("2024-01", "2024-12"))
            year_2025 = year_2025.assign_coords(time=pd.date_range("2025-01-01", "2025-12-01", freq="MS"))

            # - Concatenate with original dataset and sort by time
            ds_ = xr.concat([year_1999, ds_, year_2025], dim="time").sortby("time")

            # - Check time records
            print(f"\nTime coverage after expansion: {len(ds_.time.values)}")

            # - rm Attributes
            ds_.attrs = {}
            for v in ds_.data_vars:
                ds_[v].attrs = {}

            # - Save to NetCDF flux_ch4_prior only
            ds_["flux_ch4_prior"].fillna(0.0).to_netcdf(self.flx_apr_dir + self.flx_apr_nc)

            # - Save to NetCDF ds full
            ds_.fillna(0.0).to_netcdf(self.flx_apr_dir + self.flx_apr_nc_full)

            return ds_

        # --- write_nc
        def write_nc(ds_i):
            ds_ = ds_i.copy()

            # - Copy full year 2000 to 1999
            year_1999 = ds_.sel(time=slice('2000-01', '2000-12'))
            year_1999 = year_1999.assign_coords(time=pd.date_range('1999-01-01', '1999-12-01', freq='MS'))

            # - Copy full year 2024 to 2025
            year_2025 = ds_.sel(time=slice('2024-01', '2024-12'))
            year_2025 = year_2025.assign_coords(time=pd.date_range('2025-01-01', '2025-12-01', freq='MS'))

            # - Concatenate with original dataset and sort by time
            ds_ = xr.concat([year_1999, ds_, year_2025], dim='time').sortby('time')
            print(f'\nTime coverage after expansion: {ds_.sizes["time"]}')

            # - Remove attributes and inherited encoding
            ds_.attrs = {}
            for v in ds_.variables:
                ds_[v].attrs = {}
                ds_[v].encoding = {}

            # - Time encoding compatible with NetCDF Fortran
            encoding = {'time': {'dtype': 'int32', 'units': 'days since 1999-01-01 00:00:00', 'calendar': 'proleptic_gregorian'}}

            # - Save prior flux only
            ds_['flux_ch4_prior'].fillna(0.0).to_netcdf(self.flx_apr_dir + self.flx_apr_nc, engine='scipy', format='NETCDF3_64BIT', encoding=encoding)

            # - Save all components
            ds_.fillna(0.0).to_netcdf(self.flx_apr_dir + self.flx_apr_nc_full, engine='scipy', format='NETCDF3_64BIT', encoding=encoding)

            return ds_

        # --- write_nc
        def write_nc_0(ds_i):
            from pathlib import Path

            ds_ = ds_i.copy()

            # - Copy full year 2000 to 1999
            year_1999 = ds_.sel(time=slice('2000-01', '2000-12'))
            year_1999 = year_1999.assign_coords(time=pd.date_range('1999-01-01', '1999-12-01', freq='MS'))

            # - Copy full year 2024 to 2025
            year_2025 = ds_.sel(time=slice('2024-01', '2024-12'))
            year_2025 = year_2025.assign_coords(time=pd.date_range('2025-01-01', '2025-12-01', freq='MS'))

            # - Concatenate and sort
            ds_ = xr.concat([year_1999, ds_, year_2025], dim='time').sortby('time')
            print(f'\nTime coverage after expansion: {ds_.sizes["time"]}')

            # - Remove attributes
            ds_.attrs = {}
            for v in ds_.data_vars:
                ds_[v].attrs = {}

            # - Prepare output directory
            out_dir = Path(self.flx_apr_dir)
            out_dir.mkdir(parents=True, exist_ok=True)
            prior_path = out_dir / self.flx_apr_nc
            full_path = out_dir / self.flx_apr_nc_full

            # - Write NetCDF using h5netcdf
            print(f'Write prior: {prior_path}')
            ds_['flux_ch4_prior'].fillna(0.0).to_netcdf(str(prior_path), engine='h5netcdf', mode='w')

            print(f'Write all sectors: {full_path}')
            ds_.fillna(0.0).to_netcdf(str(full_path), engine='h5netcdf', mode='w')

            return ds_

        # - write 'flux_ch4_prior' to grd
        def write_grd(ds_i, dtype=np.float32):
            da = ds_i["flux_ch4_prior"].transpose("time", "lat", "lon")
            ntime, nlat, nlon = da.shape
            assert ntime % self.nmonth == 0, "time is not divisible by 12"
            nyear = ntime // self.nmonth
            print(f'Write *.grd for {nyear} years {nyear * 12} recs to {self.flx_apr_grd} ')

            with open(self.flx_apr_grd, "wb") as f:
                for iy in range(nyear):
                    for im in range(self.nmonth):
                        t = iy * self.nmonth + im
                        slab = da.isel(time=t).values.astype(dtype, copy=False)
                        slab = np.nan_to_num(slab, nan=0.0)
                        f.write(slab.tobytes(order="C"))

        # =======================================
        ds_1 = get_inp_nc()
        ds_2 = manage_ds(ds_1)
        check_total(ds_2)
        ds_nc = write_nc(ds_2)
        # - no need now
        # write_grd(ds_nc, dtype=np.float32)
