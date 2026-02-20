# -*- coding: utf-8 -*-
"""
Finalized:  26/xx/xx
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
import matplotlib.pyplot as plt

import _set_case


class CnvOHNc(_set_case.SetCase):
    def __init__(self):
        super().__init__()

        print('Read Input Daily INCA OH from dir: ')
        print('\t\t', self.inc_OH_inp_dir)

        # - original period
        self.start_year = 2000
        self.end_year = 2024

        self.inca_monthly_file = self.inc_OH_dir + f"INCA_DM_OH_monthly_{self.start_year}-{self.end_year}.nc"
        self.actm_OH_IAV = self.inc_OH_dir + 'tch4_actm_gcp21oh.nc'
        self.inca_final_file = self.inc_OH_dir + f"tch4_INCA_OH_t42.nc"

        # ---
        # self.get_inca_OH()

        # ---
        self.check_inca_OH()

    # --- check inca OH
    def check_inca_OH(self):
        # Open datasets
        ds_actm = xr.open_dataset(self.actm_OH_IAV, decode_times=False)
        start_date = "1980-01-01"
        new_time = pd.date_range(start=start_date, periods=ds_actm.sizes["time"], freq="MS")
        ds_actm = ds_actm.assign_coords(time=new_time)

        ds_inca = xr.open_dataset(self.inca_final_file)

        print(ds_actm["level"])
        # print(ds_inca)

        # Pick OH variable (change name if needed)
        oh_actm = ds_actm["NDOH"]
        oh_inca = ds_inca["OH"]

        MLO_LAT, MLO_LON = 19.54, -155.58+360  # Mauna Loa
        t_lev = 0.77

        # Select South Pole and North Pole (nearest grid point)
        oh_actm_SP = oh_actm.sel(level=t_lev, y=-90, method="nearest").mean(dim=("x"))
        oh_actm_NP = oh_actm.sel(level=t_lev, y=90, method="nearest").mean(dim=("x"))
        oh_actm_MLO = oh_actm.sel(level=t_lev,y=MLO_LAT, x=MLO_LON, method="nearest")

        oh_inca_SP = oh_inca.sel(level=t_lev, y=-90, method="nearest").mean(dim=("x"))
        oh_inca_NP = oh_inca.sel(level=t_lev, y=90, method="nearest").mean(dim=("x"))
        oh_inca_MLO = oh_inca.sel(level=t_lev,y=MLO_LAT, x=MLO_LON, method="nearest")

        # Plot
        plt.figure(figsize=(8, 5))

        plt.plot(oh_actm["time"], oh_actm_SP, label="ACTM SP", color="r")
        plt.plot(oh_actm["time"], oh_actm_MLO, label="ACTM MLO", color="g")
        plt.plot(oh_actm["time"], oh_actm_NP, label="ACTM NP", color="b")

        plt.plot(oh_inca["time"], oh_inca_SP, label="INCA SP", color="r", linestyle="--")
        plt.plot(oh_inca["time"], oh_inca_MLO, label="INCA MLO", color="g", linestyle="--")
        plt.plot(oh_inca["time"], oh_inca_NP, label="INCA NP", color="b", linestyle="--")

        plt.xlabel("")
        plt.ylabel("OH concentration")
        plt.title("OH Timeseries")
        plt.xlim(pd.Timestamp("2000-01-01"), pd.Timestamp("2005-12-31"))
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        # plt.show()
        path = 'D:/OneDrive - 国立大学法人千葉大学/prj_GCP_v25/plots/'
        plt.savefig(path + 'F_OH_sigma_' + str(t_lev) + '.png', format='png', dpi=1200, bbox_inches='tight')


    # --- get inca OH
    def get_inca_OH(self):
        # -
        def merge_daily_to_monthly_oh(how="mean"):
            ds_all = []

            for year in range(self.start_year, self.end_year + 1):
                f = Path(self.inc_OH_inp_dir) / f"INCA_DM_{year}_OH_scaled.nc"
                print(f"Reading {f}")
                ds = xr.open_dataset(f)

                if not pd.api.types.is_datetime64_any_dtype(ds.time_counter):
                    ds["time_counter"] = xr.decode_cf(ds).time_counter

                if how == "mean":
                    ds_mon = ds.resample(time_counter="MS").mean()
                else:
                    raise ValueError("how must be 'mean' or 'sum'")
                ds_all.append(ds_mon)

            ds_merged = xr.concat(ds_all, dim="time_counter")
            ds_merged.to_netcdf(self.inca_monthly_file)
            print(f"Saved → {self.inca_monthly_file}")

        # -
        def get_actm_param():
            ds = xr.open_dataset(self.actm_OH_IAV, decode_times=False)
            lon = ds["x"].values
            lat = ds["y"].values
            lev = ds["level"].values
            print(lon.shape, lat.shape, lev.shape)
            return lon, lat, lev

        # -
        def inca_to_actm_latlon(t_lon, t_lat):
            ds = xr.open_dataset(self.inca_monthly_file)
            ds = ds.rename({"time_counter": "time",
                            "presnivs": "level",
                            "lat": "y",
                            "lon": "x"
                            })

            # convert lon from -180:180 to 0:360 and replace coordinate
            ds = ds.assign_coords(x=(ds["x"] + 360) % 360)
            ds = ds.sortby("x")

            ds_int_ll = ds.interp(x=t_lon, y=t_lat, method="linear")
            ds_int_ll["OH"] = ds_int_ll["OH"].where(np.isfinite(ds_int_ll["OH"]))
            return ds_int_ll

        # -
        def inca_to_actm_lev(ds, actm_sigma):
            # Surface pressure (assume lowest level ~ surface)
            ps = ds["level"].max().item()

            # Compute sigma = p / ps
            sigma = ds["level"] / ps
            ds = ds.assign_coords(sigma=("level", sigma.values))
            ds = ds.swap_dims({"level": "sigma"}).drop_vars("level")

            # Interpolate to ACTM sigma grid
            ds_int_lv = ds.interp(sigma=actm_sigma, method="linear")

            # Rename sigma -> level for ACTM compatibility
            ds_int_lv = ds_int_lv.rename({"sigma": "level"})

            # Convert OH: molecules/m3 → molecules/cm3
            ds_int_lv["OH"] = ds_int_lv["OH"] / 1e6

            # Drop vmroh
            ds_int_lv = ds_int_lv.drop_vars("vmroh")

            # Metadata
            ds_int_lv["OH"].attrs["units"] = "molecules cm-3"
            ds_int_lv["OH"].attrs["long_name"] = "Hydroxyl radical"

            return ds_int_lv

        # -
        def extend_2_years(ds):

            ds_2000 = ds.sel(time=slice("2000-01-01", "2000-12-01"))
            ds_2024 = ds.sel(time=slice("2024-01-01", "2024-12-01"))

            # build new time coords using pandas (safe way)
            time_1999 = pd.date_range("1999-01-01", "1999-12-01", freq="MS")
            time_2025 = pd.date_range("2025-01-01", "2025-12-01", freq="MS")

            ds_1999 = ds_2000.copy(deep=True)
            ds_1999 = ds_1999.assign_coords(time=time_1999)

            ds_2025 = ds_2024.copy(deep=True)
            ds_2025 = ds_2025.assign_coords(time=time_2025)

            ds_ext = xr.concat([ds_1999, ds, ds_2025], dim="time")

            return ds_ext

        # =======================================
        # - to monthly file at server
        run = 0
        if run:
            merge_daily_to_monthly_oh()

        # - get actm lat/lon/lev
        lon_a, lat_a, lev_a = get_actm_param()

        # - inca at actm lat/lon
        ds_latlon = inca_to_actm_latlon(lon_a, lat_a)
        print(ds_latlon)

        # -
        ds_lev = inca_to_actm_lev(ds_latlon, lev_a)
        ds_ = extend_2_years(ds_lev)

        ds_.encoding.pop("unlimited_dims", None)
        print(ds_)
        ds_.to_netcdf(self.inca_final_file)
