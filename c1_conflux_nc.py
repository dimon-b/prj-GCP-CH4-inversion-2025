# -*- coding: utf-8 -*-
"""
Created/Corrected:    25/09/09 09:09
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import numpy as np
import pandas as pd
import xarray as xr
import _set_case


class CnvFluxNc(_set_case.SetCase):
    def __init__(self):
        super().__init__()

        print('Read Input Fluxes from dir: ')
        print('\t\t', self.inp_dir)

        # ---
        self.get_inp_flx()

    # --- get_inp_flx
    def get_inp_flx(self):
        path = self.inp_dir + self.inp_flx
        ds = xr.open_dataset(path, decode_times=False)
        print(ds)

        # - get lengths
        ntime = ds.dims["time"]
        nseason = ds.dims["time_climato"]

        # - make seasonal variables match the full time axis
        expanded = {}
        for v in [v for v in ds.data_vars if "time_climato" in ds[v].dims]:
            # tile the seasonal data to match the time dimension
            arr = ds[v]
            reps = ntime // nseason
            expanded[v] = xr.DataArray(
                np.tile(arr.values, (reps, 1, 1)),  # repeat along time axis
                dims=("time", "lat", "lon"),
                coords={"time": ds.time, "lat": ds.lat, "lon": ds.lon},
                attrs=arr.attrs)

        # - keep the original "time" variables as they are
        time_vars = {v: ds[v] for v in ds.data_vars if "time" in ds[v].dims}

        # - join
        ds_expanded = xr.Dataset({**time_vars, **expanded})

        # - add total
        sum_all = sum([ds_expanded[v] for v in ds_expanded.data_vars])
        ds_expanded["flux_ch4_total"] = sum_all

        # - Make sure time is datetime
        if not np.issubdtype(ds_expanded.time.dtype, np.datetime64):
            ds_expanded['time'] = xr.decode_cf(ds_expanded[['time']]).time

        # - Compute grid cell area (m²)
        dlat = np.deg2rad(ds_expanded.lat[1] - ds_expanded.lat[0])
        dlon = np.deg2rad(ds_expanded.lon[1] - ds_expanded.lon[0])
        lat_rad = np.deg2rad(ds_expanded.lat)
        cell_area = (self.R ** 2 * dlat * dlon * np.cos(lat_rad)).broadcast_like(ds_expanded.flux_ch4_soils)

        # - check cell_area_total
        cell_area_total = cell_area.sum(dim=("lat", "lon"))
        print('')
        print(f"Check Total Cell Area = {cell_area_total.values[0]:.2e}")
        print('')

        # - Compute seconds per month
        time = pd.to_datetime(ds_expanded.time.values)
        days_in_month = time.days_in_month.values
        seconds_per_month = days_in_month * 24 * 3600

        # - print results in table
        results = []
        for var in ds_expanded.data_vars:
            monthly_flux_per_cell = ds_expanded[var] * cell_area * seconds_per_month[:, np.newaxis, np.newaxis]
            global_monthly_flux = monthly_flux_per_cell.sum(dim=("lat", "lon"))
            global_annual_Tg = global_monthly_flux.groupby("time.year").sum(dim="time") / 1e9  # convert to Tg/year
            for year, value in zip(global_annual_Tg.year.values, global_annual_Tg.values):
                results.append({"variable": var[9:], "year": int(year), "global_annual_Tg": float(value)})

        df_all = pd.DataFrame(results)
        df_table = df_all.pivot(index="year", columns="variable", values="global_annual_Tg")
        cols = [c for c in df_table.columns if c != "total"] + ["total"]
        df_table = df_table[cols]
        print(df_table.round(4))

        # - write 'flux_ch4_total' to NetCDF
        ds_expanded_ = ds_expanded.rename({"flux_ch4_total": "fch4"})

        # Copy full year 2000 → 1999
        year_1999 = ds_expanded_.sel(time=slice("2000-01", "2000-12"))
        year_1999 = year_1999.assign_coords(time=pd.date_range("1999-01-01", "1999-12-01", freq="MS")                                                )

        # Copy full year 2024 → 2025
        year_2025 = ds_expanded_.sel(time=slice("2024-01", "2024-12"))
        year_2025 = year_2025.assign_coords(time=pd.date_range("2025-01-01", "2025-12-01", freq="MS")                                                )

        # Concatenate with original dataset and sort by time
        ds_expanded_ = xr.concat([year_1999, ds_expanded_, year_2025], dim="time").sortby("time")

        # --- Check time records ---
        print("Time coverage after expansion:")
        print(len(ds_expanded_.time.values))

        # Save to NetCDF
        ds_expanded_["fch4"].fillna(0.0).to_netcdf(self.out_dir + "gcp25_c1.nc")

