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

        # - original period
        self.start_year = 2000
        self.end_year = 2024

        self.inca_monthly_file = self.inc_OH_dir + f"INCA_DM_OH_monthly_{self.start_year}-{self.end_year}.nc"
        self.actm_OH_IAV = self.inc_OH_dir + 'tch4_actm_gcp21oh.nc'
        self.actm_OH_INCA = self.inc_OH_dir + 'tch4_INCA_OH_t42.nc'
        # self.actm_OH_CYC = self.inc_OH_dir + 'tch4_actm_oh1.nc'
        # self.actm_OH_INCA_gt3 = self.inc_OH_dir + 'gt2nc.nc'
        # self.actm_OH_INCA_bin = self.inc_OH_dir + 'tch4_INCA_OH_t42.bin'

        # --- get INCA OH
        self.get_inca_OH()

        # --- check INCA OH
        # self.check_incaOH_ts()
        # self.check_incaOH_latlev()
        # self.check_incaOH_latlon()

    def check_incaOH_latlon(self, year=2005, month=1, sigma=0.9995):
        ds_actm = xr.open_dataset(self.actm_OH_IAV, decode_times=False)
        start_date = "1980-01-01"
        new_time = pd.date_range(start=start_date, periods=ds_actm.sizes["time"], freq="MS")
        ds_actm = ds_actm.assign_coords(time=new_time)

        ds_inca = xr.open_dataset(self.actm_OH_INCA)

        # --- Select month
        sel_time = f"{year:04d}-{month:02d}-01"
        actm_mon = ds_actm["NDOH"].sel(time=sel_time)
        inca_mon = ds_inca["OH"].sel(time=sel_time)

        # --- Select nearest sigma level
        actm_map = actm_mon.sel(level=sigma, method="nearest")
        inca_map = inca_mon.sel(level=sigma, method="nearest")

        lat = actm_map["y"].values
        lon = actm_map["x"].values

        # --- Ensure ordering: (lat, lon)
        actm_plot = actm_map.transpose("y", "x").values
        inca_plot = inca_map.transpose("y", "x").values

        # --- Common color scale (reuse your style)
        levels = np.linspace(0, 1e6 / 1000, 11)

        # --- Plot
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)

        cs1 = axes[0].contourf(lon, lat, actm_plot, levels=levels, cmap="jet", extend="both")
        axes[0].set_title(f"ACTM OH ({year}-{month:02d}, σ={actm_map.level.values:.2f})")
        axes[0].set_xlabel("Longitude")
        axes[0].set_ylabel("Latitude")

        cs2 = axes[1].contourf(lon, lat, inca_plot, levels=levels, cmap="jet", extend="both")
        axes[1].set_title(f"INCA OH ({year}-{month:02d}, σ={inca_map.level.values:.2f})")
        axes[1].set_xlabel("Longitude")

        # --- Horizontal colorbar
        cbar = fig.colorbar(cs1, ax=axes.ravel().tolist(),
                            orientation="horizontal", pad=0.15, aspect=40, shrink=0.8)
        cbar.set_ticks(levels)
        cbar.set_label("OH")
        plt.show()

    # --- check inca OH latlev
    def check_incaOH_latlev(self, year=2005, month=8):
        ds_actm = xr.open_dataset(self.actm_OH_IAV, decode_times=False)
        start_date = "1980-01-01"
        new_time = pd.date_range(start=start_date, periods=ds_actm.sizes["time"], freq="MS")
        ds_actm = ds_actm.assign_coords(time=new_time)

        ds_inca = xr.open_dataset(self.actm_OH_INCA)

        # --- Select month
        sel_time = f"{year:04d}-{month:02d}-01"
        actm_mon = ds_actm["NDOH"].sel(time=sel_time)
        inca_mon = ds_inca["OH"].sel(time=sel_time)

        # --- Average over longitude
        actm_latlev = actm_mon.mean(dim="x")
        inca_latlev = inca_mon.mean(dim="x")

        # --- Select sigma range (surface to mid-troposphere)
        actm_latlev = actm_latlev.sel(level=slice(0.99, 0.75))
        inca_latlev = inca_latlev.sel(level=slice(0.99, 0.75))

        lat = actm_latlev["y"].values
        lev = actm_latlev["level"].values

        # --- Force consistent ordering: (level, lat)
        actm_plot = actm_latlev.transpose("level", "y").values
        inca_plot = inca_latlev.transpose("level", "y").values

        # --- Common color scale
        vmin = float(min(np.nanmin(actm_plot), np.nanmin(inca_plot)))
        vmax = float(max(np.nanmax(actm_plot), np.nanmax(inca_plot)))
        levels = np.linspace(0, 1e6, 11)

        # --- Plot
        fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)

        cs1 = axes[0].contourf(lat, lev, actm_plot, levels=levels, cmap="jet", extend="both")
        axes[0].set_title(f"ACTM OH ({year}-{month:02d})")
        axes[0].set_xlabel("Latitude")
        axes[0].set_ylabel("Sigma level")

        cs2 = axes[1].contourf(lat, lev, inca_plot, levels=levels, cmap="jet", extend="both")
        axes[1].set_title(f"INCA OH ({year}-{month:02d})")
        axes[1].set_xlabel("Latitude")

        # --- Invert sigma axis (1.0 at bottom)
        axes[0].set_ylim(lev.max(), lev.min())
        axes[1].set_ylim(lev.max(), lev.min())

        # --- Horizontal colorbar
        cbar = fig.colorbar(cs1, ax=axes.ravel().tolist(), orientation="horizontal",
                            pad=0.15, aspect=40, shrink=0.8)
        cbar.set_ticks(levels)
        cbar.set_label("OH")
        plt.show()

    # --- check inca OH ts
    def check_incaOH_ts(self):
        # gt2nc
        ds_actm = xr.open_dataset(self.actm_OH_INCA_gt3, decode_times=False)
        start_date = "1999-01-01"
        new_time = pd.date_range(start=start_date, periods=ds_actm.sizes["time"], freq="MS")
        ds_actm = ds_actm.assign_coords(time=new_time)

        # Open CYC datasets
        # ds_actm = xr.open_dataset(self.actm_OH_CYC, decode_times=False)
        # start_date = "1999-01-01"
        # new_time = pd.date_range(start=start_date, periods=ds_actm.sizes["time"], freq="MS")
        # ds_actm = ds_actm.assign_coords(time=new_time)
        # print(ds_actm)

        # # Open IAV datasets
        # ds_actm = xr.open_dataset(self.actm_OH_IAV, decode_times=False)
        # start_date = "1980-01-01"
        # new_time = pd.date_range(start=start_date, periods=ds_actm.sizes["time"], freq="MS")
        # ds_actm = ds_actm.assign_coords(time=new_time)

        # Open INCA datasets
        ds_inca = xr.open_dataset(self.actm_OH_INCA)
        print(ds_inca)

        # print(ds_actm["level"])
        # print(ds_inca)

        # Pick OH variable (change name if needed)
        oh_actm = ds_actm["NDOH"]
        oh_inca = ds_inca["OH"]

        l_lat, t_lon = 19.54, -155.58 + 360  # Mauna Loa
        # l_lat, t_lon = -80, 3950
        t_lev = 0.99

        # Select South Pole and North Pole (nearest grid point)
        oh_actm_SP = oh_actm.sel(level=t_lev, y=-80, method="nearest").mean(dim=("x"))
        oh_actm_NP = oh_actm.sel(level=t_lev, y=90, method="nearest").mean(dim=("x"))
        oh_actm_MLO = oh_actm.sel(level=t_lev, y=l_lat, x=t_lon, method="nearest")

        oh_inca_SP = oh_inca.sel(level=t_lev, y=-80, method="nearest").mean(dim=("x"))
        oh_inca_NP = oh_inca.sel(level=t_lev, y=90, method="nearest").mean(dim=("x"))
        oh_inca_MLO = oh_inca.sel(level=t_lev, y=l_lat, x=t_lon, method="nearest")

        # Plot
        plt.figure(figsize=(8, 5))

        plt.plot(oh_actm["time"], oh_actm_SP, label="ACTM SP", color="r")
        plt.plot(oh_actm["time"], oh_actm_NP, label="ACTM NP", color="b")
        plt.plot(oh_actm["time"], oh_actm_MLO, label="ACTM T", color="g")

        plt.plot(oh_inca["time"], oh_inca_SP, label="INCA SP", color="r", linestyle="--")
        plt.plot(oh_inca["time"], oh_inca_NP, label="INCA NP", color="b", linestyle="--")
        plt.plot(oh_inca["time"], oh_inca_MLO, label="INCA T", color="g", linestyle="--")

        plt.xlabel("")
        plt.ylabel("OH concentration")
        plt.title("OH Timeseries")
        plt.xlim(pd.Timestamp("1998-12-01"), pd.Timestamp("2000-01-31"))
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        # path = 'D:/OneDrive - 国立大学法人千葉大学/prj_GCP_v25/plots/'
        # plt.savefig(path + 'F_OH_sigma_' + str(t_lev) + '.png', format='png', dpi=1200, bbox_inches='tight')

    # --- get inca OH
    def get_inca_OH(self):
        # -
        def merge_daily_to_monthly_oh():
            print('Read Input Daily INCA OH from dir: ')
            print('\t\t', self.inc_OH_inp_dir)

            ds_all = []

            for year in range(self.start_year, self.end_year + 1):
                f = Path(self.inc_OH_inp_dir) / f"INCA_DM_{year}_OH_scaled.nc"
                print(f"Reading {f}")
                ds = xr.open_dataset(f)

                if not pd.api.types.is_datetime64_any_dtype(ds.time_counter):
                    ds["time_counter"] = xr.decode_cf(ds).time_counter

                ds_mon = ds.resample(time_counter="MS").mean()
                ds_all.append(ds_mon)

            ds_merged = xr.concat(ds_all, dim="time_counter")
            ds_merged.to_netcdf(self.inca_monthly_file)
            print(f"Saved → {self.inca_monthly_file}")

        # -
        def get_actm_ds():
            ds_a = xr.open_dataset(self.actm_OH_IAV, decode_times=False)
            start_date = "1980-01-01"
            new_time = pd.date_range(start=start_date, periods=ds_a.sizes["time"], freq="MS")
            ds_a = ds_a.assign_coords(time=new_time)

            lon = ds_a["x"].values
            lat = ds_a["y"].values
            lev = ds_a["level"].values
            print('\t ACTM grid shape: ')
            return ds_a, lon, lat, lev,

        # -
        def get_inca_ds():
            ds_i = xr.open_dataset(self.inca_monthly_file)

            # use actm names
            ds_i = ds_i.rename({"time_counter": "time",
                                "presnivs": "level",
                                "lat": "y",
                                "lon": "x"
                                })
            return ds_i

        # -
        def regrid_oh(ds_i, ds_a):

            # Fix longitude
            ds_i = ds_i.assign_coords(
                x=(ds_i.x % 360)
            ).sortby("x")

            # Fix latitude orientation
            if ds_i.y[0] > ds_i.y[-1]:
                ds_i = ds_i.sortby("y")

            # Horizontal interpolation
            ds_out = ds_i.interp(
                y=ds_a.y,
                x=ds_a.x,
                method="linear"
            )

            return ds_out

        def inca_to_actm_lev(ds_inca, actm_sigma):
            # Surface pressure (assume lowest level ~ surface)
            ps = ds_inca["level"].max().item()

            # Compute sigma = p / ps
            sigma = ds_inca["level"] / ps
            ds_inca = ds_inca.assign_coords(sigma=("level", sigma.values))
            ds_inca = ds_inca.swap_dims({"level": "sigma"}).drop_vars("level")

            # Interpolate to ACTM sigma grid
            ds_int_lv = ds_inca.interp(sigma=actm_sigma, method="linear")

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
        def inca_to_actm_lev_pressure(ds_inca, ds_actm):

            # --- ACTM sigma levels (final vertical coordinate)
            actm_sigma = ds_actm["level"]  # 1D: (level)

            # --- Surface pressure from INCA (per column)
            ps = ds_inca["level"].max(dim="level")  # (time, y, x)

            # --- Convert ACTM sigma to pressure (Pa)
            p_actm = actm_sigma * ps

            # --- Interpolate INCA to ACTM pressure levels
            ds_int_lv = ds_inca.interp(level=p_actm, method="linear")

            # --- Overwrite vertical coordinate with ACTM sigma
            ds_int_lv = ds_int_lv.assign_coords(level=actm_sigma)

            # --- Convert OH: molecules/m3 → molecules/cm3
            ds_int_lv["OH"] = ds_int_lv["OH"] / 1e6

            # --- Drop vmroh safely
            if "vmroh" in ds_int_lv:
                ds_int_lv = ds_int_lv.drop_vars("vmroh")

            # --- Metadata
            ds_int_lv["level"].attrs.update({"standard_name": "atmosphere_sigma_coordinate",
                                             "long_name": "ACTM sigma level"
                                             })
            ds_int_lv["OH"].attrs.update({"units": "molecules cm-3",
                                          "long_name": "Hydroxyl radical"
                                          })

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

        # - get actm
        ds_a, lon_a, lat_a, lev_a = get_actm_ds()

        # - get inca
        ds_i = get_inca_ds()

        # - lat/lon regrid
        ds_ia_latlon = regrid_oh(ds_i, ds_a)

        # - lev regrid
        ds_lev = inca_to_actm_lev_pressure(ds_ia_latlon, ds_a)

        # -
        ds_ = extend_2_years(ds_lev)

        # -
        ds_["OH"] = ds_["OH"].astype("float32")
        ds_.to_netcdf(self.actm_OH_INCA, unlimited_dims=[])