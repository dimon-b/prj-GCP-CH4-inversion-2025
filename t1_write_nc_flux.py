# -*- coding: utf-8 -*-
"""
Created:    22/01/26 19:17
Project:    Write total fCH4 *.nc for GCP2021 as originally created by Naveen Negi
@author:    Dmitry Belikov
"""
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import re

import _set_case


class WriteNcFlux(_set_case.SetCase):

    def __init__(self):
        super().__init__()

        self.invcases = ['inv1']
        self.flxcases = ['gcp2021_v2_soil0_inv1']
        self.dirs = ['../results2025/CYC/', '../results2025/INCA/']
        self.unp = 'p30'
        self.unx = 'ctl'

        #
        for d in self.dirs:
            self.write_submission_flux(d)

    # --- write flux categories
    def write_submission_flux(self, dir_):
        # ---
        def get_post_flx_grd(unp, unx, invc_):
            # -
            def read_flx_bin(bfile):
                size_2d = self.d1_nlon * self.d1_nlat * self.nmonth
                with open(bfile, 'rb') as f:
                    fl2d = np.fromfile(f, dtype='<f4', count=size_2d * self.nyear_nc
                                       ).reshape(self.nyear_nc, self.nmonth, self.d1_nlat, self.d1_nlon)
                print(f'\t\tpost flux shape from grd file: {fl2d.shape}')
                return fl2d

            # - flux shape correction [yr, mn, :, :] -> [time, :, :]
            def cor_shape(flx_):
                flx1_ = np.zeros([self.nyear_nc * self.nmonth, self.d1_nlat, self.d1_nlon])
                for yr in np.arange(self.nyear_nc):
                    for mn in np.arange(self.nmonth):
                        ss = yr * self.nmonth + mn
                        flx1_[ss, :, :] = flx_[yr, mn, :, :]
                print(f'\t\tpost flux shape after reshaping: {flx1_.shape}, '
                      f'time length: {self.nyear_nc * self.nmonth}, years: {self.nyear_nc}, months: {self.nmonth}')
                return flx1_

            # ===========================================================
            # - post
            i_file = dir_ + 'flux2d/' + unp + '/' + self.icase + '_' + unx + '_' + invc_ + '.grd'
            print('\tPosterior flux reading: ', invc_, i_file)
            flx = read_flx_bin(i_file)
            flx1 = cor_shape(flx)
            return flx1

        # - get_prior_full_flx_nc
        def get_prior_full_flx_nc():
            path = self.flx_apr_dir + self.flx_apr_nc_full
            ds_1 = xr.open_dataset(path, engine="h5netcdf", decode_times=False)
            return ds_1

        # - get_joined_apr2post
        def get_joined_apr2post(apr, pst):
            # - soil from apr
            apr_soil = apr["flux_ch4_soils"]

            # - extend pst
            pst_ext = np.concatenate([pst[:12], pst, pst[-12:], ], axis=0, )

            # - pst_ext into a DataArray
            pst_da = xr.DataArray(xr.where(pst_ext < 0, 0.0, pst_ext), dims=apr_soil.dims, coords=apr_soil.coords,
                                  name="pst")
            assert pst_da.shape == apr_soil.shape
            assert pst_da.time.equals(apr_soil.time)

            pst_flux = xr.Dataset({
                "fch4_soils": apr["flux_ch4_soils"],
                "fch4_total_prior": apr["flux_ch4_prior"],
                "fch4_total_prior_soil0": apr["flux_ch4_prior_soil0"],
                "fch4_total_post": pst_da,
                "fch4_total_post_soil0": pst_da - apr["flux_ch4_soils"],
                "fch4_corr": pst_da - apr["flux_ch4_prior"],
            })
            if not np.issubdtype(pst_flux.time.dtype, np.datetime64):
                pst_flux['time'] = xr.decode_cf(pst_flux[['time']]).time

            return pst_flux

        # - check_total
        def check_total(ds_i):
            # - Grid-cell area [m2]
            dlat = np.deg2rad(ds_i.lat[1] - ds_i.lat[0])
            dlon = np.deg2rad(ds_i.lon[1] - ds_i.lon[0])
            lat_rad = np.deg2rad(ds_i.lat)

            cell_area_2d = self.R ** 2 * dlat * dlon * np.cos(lat_rad)
            cell_area = cell_area_2d.broadcast_like(ds_i["fch4_total_prior"])
            # cell_area = xr.DataArray(
            #     cell_area_2d,
            #     coords=ds_i["fch4_total_prior"].isel(time=0).coords,
            #     dims=ds_i["fch4_total_prior"].isel(time=0).dims,
            # ).broadcast_like(ds_i["fch4_total_prior"])

            # - Sanity check: total Earth surface area
            earth_area = cell_area.isel(time=0).sum(dim=("lat", "lon"))
            print("")
            print(f"Check total Earth area = {earth_area.item():.2e} m2")
            print("Expected ≈ 5.10e14 m2")
            print("")

            # - Seconds per month [s]
            time = pd.to_datetime(ds_i.time.values)
            ds_i["seconds_per_month"] = (ds_i.time.dt.days_in_month * 24 * 3600).rename("seconds_per_month")
            ds_i["seconds_per_month"].attrs.update(units="s", description="Number of seconds in each month", )

            # - Time since epoch [hours since 1970-01-01 00:00:00]
            epoch = pd.Timestamp("1970-01-01 00:00:00")

            ds_i["time_hours_since_1970"] = xr.DataArray((time - epoch) / pd.Timedelta(hours=1),
                                                         coords={"time": ds_i.time},
                                                         dims=("time",),
                                                         name="time_hours_since_1970",
                                                         )

            ds_i["time_hours_since_1970"].attrs.update(units="hours since 1970-01-01 00:00:00",
                                                       description="Time coordinate expressed as hours since Unix epoch",
                                                       )

            # - Flux variables to integrate
            FLUX_VARS = [v for v in ds_i.data_vars if v.startswith("fch4_")]

            records = []

            for var in FLUX_VARS:
                # - Monthly flux (time integration only) [kg m-2]
                out_name = f"monthly_{var}"
                ds_i[out_name] = ds_i[var] * ds_i["seconds_per_month"]
                ds_i[out_name].attrs.update(units="kg m-2", description=f"Monthly integrated flux from {var}")

                # - Global monthly total [kg]
                # noinspection PyTypeChecker
                monthly_global: xr.DataArray = (ds_i[out_name] * cell_area).sum(dim=("lat", "lon"))
                # - Global annual total [Tg yr-1]
                annual_global = (monthly_global.groupby("time.year").sum(dim="time") / 1e9)

                for year, value in zip(annual_global.year.values, annual_global.values):
                    records.append({"year": int(year),
                                    "variable": var.replace("fch4_", ""),
                                    "Tg_CH4_yr": float(value),
                                    })

            # - Output table
            df = pd.DataFrame(records)

            table = df.pivot(index="year", columns="variable", values="Tg_CH4_yr", )

            # - Put totals at the end
            total_cols = [c for c in table.columns if c.lower().startswith("total")]
            other_cols = [c for c in table.columns if c not in total_cols]
            table = table[other_cols + total_cols]

            print(table.round(4))

            def plot_fluxes_and_corr(table_):
                # - Flux columns (left axis)
                flux_cols = ["post", "post_soil0", "prior", "prior_soil0"]

                # - Figure
                fig, ax1 = plt.subplots(figsize=(9, 5))

                # - Flux time series
                table_[flux_cols].plot(ax=ax1, marker="o")
                ax1.set_xlabel("Year")
                ax1.set_ylabel("Total fluxes, (Tg CH4 yr⁻¹)")
                ax1.grid(True)
                ax1.set_ylim(450, 650)

                # - Correction axis
                ax2 = ax1.twinx()
                ax2.set_ylabel("Correction, soil (Tg CH4 yr⁻¹)")
                ax2.set_ylim(-60, 60)

                # - Correction
                ax2.plot(table_.index,
                         table_["corr"],
                         linestyle="--",
                         linewidth=2,
                         label="Correction",
                         )

                # - Optional soils series (if present)
                if "soils" in table_.columns:
                    ax2.plot(table_.index,
                             table_["soils"],
                             linestyle=":",
                             linewidth=2,
                             label="soils",
                             )

                # - Combined legend
                lines1, labels1 = ax1.get_legend_handles_labels()
                lines2, labels2 = ax2.get_legend_handles_labels()
                ax1.legend(lines1 + lines2,
                           labels1 + labels2, ncol=3,
                           loc="lower right",
                           )

                plt.tight_layout()
                plt.show()

            # plot_fluxes_and_corr(table)

            return ds_i

        # - write_1nc_total
        def write_1nc_total(invc_, ds_1):

            file_out = ''
            if 'CYC' in dir_:
                file_out = dir_ + 'nc_out/' + '/MIROC4-ACTM_totflux_GMB_SURF_OH_Transcom.nc'
            elif 'INCA' in dir_:
                file_out = dir_ + 'nc_out/' + '/MIROC4-ACTM_totflux_GMB_SURF_OH_INCA.nc'

            nc = netCDF4.Dataset(file_out, 'w', format='NETCDF4')
            nc.description = 'Net prior and posteor CH4 emissions resulted from the surface based inversion for GCP-CH4, 2025. The results are produced at JAMSTEC, Japan. For details please contact at prabir@jamstec.go.jp and d.belikov@chiba-u.jp.'
            nc.Institution = 'Japan Agency for Marine-Earth Science and Technology (JAMSTEC) and Chiba University'
            nc.Contact = 'd.belikov@chiba-u.jp and prabir@jamstec.go.jp'
            nc.Disclaimer = 'This data is created for GCP-CH4 2025.'
            nc.History = 'April 2026'
            nc.Additional_Info = 'MIROC4.0-based Atmospheric Chemistry-Transport Model (MIROC4-ACTM) is used for forward simulations (Prabir et al., 2018) and Time-independent Bayesian inversion is used for optimising flux (Patra et al., 2016)'

            # --- set dimensions
            nc.createDimension('time', None)
            nc.createDimension('lat', self.d1_nlat)
            nc.createDimension('lon', self.d1_nlon)

            # --- set variables
            time = nc.createVariable('time', 'f4', ('time',))
            time.units = 'hours since 1970-01-01 00:00:00'
            time.calendar = 'Gregorian'
            lat = nc.createVariable('lat', 'f4', ('lat',))
            lat.units = 'degrees north from -90'
            lat.long_name = 'latitudes'
            lon = nc.createVariable('lon', 'f4', ('lon',))
            lon.units = 'degrees east from -180'
            lon.long_name = 'longitudes'
            cellarea = nc.createVariable('carea', 'f4', ('lat',))
            fch4_prior = nc.createVariable('fch4_prior', 'f4', ('time', 'lat', 'lon'))
            fch4_post = nc.createVariable('fch4_post', 'f4', ('time', 'lat', 'lon'))

            # --- set data
            lat[:] = self.d1_lats
            lon[:] = self.d1_lons
            cellarea[:] = self.garia_d1
            fch4_prior[:, :, :] = ds_1['monthly_fch4_total_prior_soil0'].values
            fch4_post[:, :, :] = ds_1['monthly_fch4_total_post_soil0'].values
            # print(ds_1)
            # exit()

            # --- date-time
            time[:] = ds_1["time_hours_since_1970"]

            # --- units
            cellarea.units = 'm2'
            fch4_prior.units = 'g-CH4/m2/month'
            fch4_post.units = 'g-CH4/m2/month'

            # ---
            fch4_prior.valid_range = np.array((np.min(fch4_prior), np.max(fch4_prior)))
            fch4_post.valid_range = np.array((np.min(fch4_post), np.max(fch4_post)))

            print('\n SUCCESS writing total flux to >>> ', file_out)
            print('\t Check valid range')

            def pr(name, arr):
                v = arr.valid_range
                print(f"      Check {name:<15} [{v[0]:.5f}, {v[1]:5f}]")

            pr('prior tot', fch4_prior)
            pr('post  tot', fch4_post)
            nc.close()

        # - get_post_categ
        def get_post_categ(flx_prior_, flx_join_):

            flx_prior_ = flx_prior_.copy()
            flx_join_ = flx_join_.copy()
            # print(flx_prior)

            # - Harmonize time coordinate
            if flx_join_.time.dtype != flx_prior_.time.dtype:
                flx_join_ = flx_join_.assign_coords(time=flx_prior_.time)

            # - Detect categories
            cat_pattern = re.compile(r"^flux_ch4_(.+)$")
            CATEGORIES = sorted([m.group(1)
                                 for v in flx_prior_.data_vars
                                 if (m := cat_pattern.match(v)) and m.group(1) not in ["total", "prior", "prior_soil0"]
                                 ])

            if not CATEGORIES:
                raise RuntimeError("No flux_ch4_<cat> variables detected")

            # print("Detected categories:")
            # for c in CATEGORIES:
            #     print(f"  - {c}")

            # - Totals
            f_tot_prior = flx_prior_["flux_ch4_prior_soil0"]
            f_tot_post = flx_join_["fch4_total_post_soil0"]

            # - Protect against division by zero
            # f_tot_prior_safe = xr.where(np.abs(f_tot_prior) > 1e-15, f_tot_prior, np.nan)

            # - Init output dataset
            ds_out = xr.Dataset(coords=flx_join_.coords)
            ds_out["fch4_total_prior"] = flx_prior_["flux_ch4_prior_soil0"]
            ds_out["fch4_total_post"] = flx_join_["fch4_total_post_soil0"]

            # - Copy PRIOR categories
            for cat in CATEGORIES:
                ds_out[f"fch4_{cat}_prior"] = flx_prior_[f"flux_ch4_{cat}"]

            # - Redistribute POSTERIOR categories
            # print("\nProcessing categories:")
            for cat in CATEGORIES:
                # print(f"  - {cat}")
                # ratio = flx_prior[f"flux_ch4_{cat}"] / f_tot_prior_safe
                # ds_out[f"fch4_{cat}_post"] = f_tot_post * ratio
                # ds_out[f"fch4_{cat}_post"] = xr.where(np.abs(f_tot_prior) > 1e-15, f_tot_post * ratio, 0.0)
                mask = np.abs(f_tot_prior) > 1e-15
                ratio = xr.where(mask, flx_prior_[f"flux_ch4_{cat}"] / f_tot_prior, 0.0)
                ds_out[f"fch4_{cat}_post"] = f_tot_post * ratio

            # - soils
            ds_out["fch4_soils_post"] = flx_prior_[f"flux_ch4_soils"]

            # - Mass conservation check
            post_stack = xr.concat([ds_out[f"fch4_{c}_post"] for c in CATEGORIES], dim="category")
            sum_post = post_stack.sum("category", skipna=True)
            imbalance = sum_post - ds_out["fch4_total_post"]
            mean_imb = float(imbalance.mean(skipna=True))
            max_imb = float(np.abs(imbalance).max(skipna=True))

            print("\nMass conservation check:")
            print(f"  Mean imbalance : {mean_imb:.3e}")
            print(f"  Max  imbalance : {max_imb:.3e}")

            tol = 1e-6
            if max_imb > tol:
                print("WARNING: mass conservation tolerance exceeded")

            if not np.issubdtype(ds_out.time.dtype, np.datetime64):
                ds_out['time'] = xr.decode_cf(ds_out[['time']]).time
            return ds_out

        # --- Write the netcdf file
        def write_1nc_cat(invc_, ds_1):

            file_out = ''
            if 'CYC' in dir_:
                file_out = dir_ + 'nc_out/' + '/MIROC4-ACTM_catflux_GMB_SURF_OH_Transcom.nc'
            elif 'INCA' in dir_:
                file_out = dir_ + 'nc_out/' + '/MIROC4-ACTM_catflux_GMB_SURF_OH_INCA.nc'

            nc = netCDF4.Dataset(file_out, 'w', format='NETCDF4')
            nc.description = 'Net prior and posteor CH4 emissions resulted from the surface based inversion for GCP-CH4, 2025. The results are produced at JAMSTEC, Japan. For details please contact at prabir@jamstec.go.jp and d.belikov@chiba-u.jp.'
            nc.Institution = 'Japan Agency for Marine-Earth Science and Technology (JAMSTEC) and Chiba University'
            nc.Contact = 'd.belikov@chiba-u.jp and prabir@jamstec.go.jp'
            nc.Disclaimer = 'This data is created for GCP-CH4 2025.'
            nc.History = 'April 2026'
            nc.Additional_Info = 'MIROC4.0-based Atmospheric Chemistry-Transport Model (MIROC4-ACTM) is used for forward simulations (Prabir et al., 2018) and Time-independent Bayesian inversion is used for optimising flux (Patra et al., 2016)'

            # --- set dimensions
            nc.createDimension('time', None)
            nc.createDimension('lat', self.d1_nlat)
            nc.createDimension('lon', self.d1_nlon)

            # --- set variables
            time = nc.createVariable('time', 'f4', ('time',))
            time.units = 'hours since 1970-01-01 00:00:00'
            time.calendar = 'Gregorian'
            lat = nc.createVariable('lat', 'f4', ('lat',))
            lat.units = 'degrees north from -90'
            lat.long_name = 'latitudes'
            lon = nc.createVariable('lon', 'f4', ('lon',))
            lon.units = 'degrees east from -180'
            lon.long_name = 'longitudes'
            # cellarea = nc.createVariable('carea', 'f4', ('lat', 'lon'))
            cellarea = nc.createVariable('carea', 'f4', ('lat',))

            fch4_tot_prior = nc.createVariable('fch4_tot_prior', 'f4', ('time', 'lat', 'lon'))
            fch4_tot_post = nc.createVariable('fch4_tot_post', 'f4', ('time', 'lat', 'lon'))
            fch4_wet_prior = nc.createVariable('fch4_wet_prior', 'f4', ('time', 'lat', 'lon'))
            fch4_wet_post = nc.createVariable('fch4_wet_post', 'f4', ('time', 'lat', 'lon'))
            fch4_bb_prior = nc.createVariable('fch4_bb_prior', 'f4', ('time', 'lat', 'lon'))
            fch4_bb_post = nc.createVariable('fch4_bb_post', 'f4', ('time', 'lat', 'lon'))
            fch4_biofuel_prior = nc.createVariable('fch4_biofuel_prior', 'f4', ('time', 'lat', 'lon'))
            fch4_biofuel_post = nc.createVariable('fch4_biofuel_post', 'f4', ('time', 'lat', 'lon'))
            fch4_oilgas_prior = nc.createVariable('fch4_oilgas_prior', 'f4', ('time', 'lat', 'lon'))
            fch4_oilgas_post = nc.createVariable('fch4_oilgas_post', 'f4', ('time', 'lat', 'lon'))
            fch4_coal_prior = nc.createVariable('fch4_coal_prior', 'f4', ('time', 'lat', 'lon'))
            fch4_coal_post = nc.createVariable('fch4_coal_post', 'f4', ('time', 'lat', 'lon'))
            fch4_agri_prior = nc.createVariable('fch4_agri_prior', 'f4', ('time', 'lat', 'lon'))
            fch4_agri_post = nc.createVariable('fch4_agri_post', 'f4', ('time', 'lat', 'lon'))
            fch4_waste_prior = nc.createVariable('fch4_waste_prior', 'f4', ('time', 'lat', 'lon'))
            fch4_waste_post = nc.createVariable('fch4_waste_post', 'f4', ('time', 'lat', 'lon'))
            fch4_termite_prior = nc.createVariable('fch4_termite_prior', 'f4', ('time', 'lat', 'lon'))
            fch4_termite_post = nc.createVariable('fch4_termite_post', 'f4', ('time', 'lat', 'lon'))
            fch4_geol_prior = nc.createVariable('fch4_geol_prior', 'f4', ('time', 'lat', 'lon'))
            fch4_geol_post = nc.createVariable('fch4_geol_post', 'f4', ('time', 'lat', 'lon'))
            fch4_oce_prior = nc.createVariable('fch4_oce_prior', 'f4', ('time', 'lat', 'lon'))
            fch4_oce_post = nc.createVariable('fch4_oce_post', 'f4', ('time', 'lat', 'lon'))
            fch4_soils_prior = nc.createVariable('fch4_soils_prior', 'f4', ('time', 'lat', 'lon'))
            fch4_soils_post = nc.createVariable('fch4_soils_post', 'f4', ('time', 'lat', 'lon'))

            # --- set data
            lat[:] = self.d1_lats
            lon[:] = self.d1_lons
            cellarea[:] = self.garia_d1

            fch4_tot_prior[:, :, :] = ds_1['monthly_fch4_total_prior'].values
            fch4_tot_post[:, :, :] = ds_1['monthly_fch4_total_post'].values

            fch4_wet_prior[:, :, :] = ds_1['monthly_fch4_wetlands_prior'].values
            fch4_wet_post[:, :, :] = ds_1['monthly_fch4_wetlands_post'].values
            fch4_bb_prior[:, :, :] = ds_1['monthly_fch4_biomass_prior'].values
            fch4_bb_post[:, :, :] = ds_1['monthly_fch4_biomass_post'].values
            fch4_biofuel_prior[:, :, :] = ds_1['monthly_fch4_biofuels_prior'].values
            fch4_biofuel_post[:, :, :] = ds_1['monthly_fch4_biofuels_post'].values
            fch4_oilgas_prior[:, :, :] = ds_1['monthly_fch4_oilgasind_prior'].values
            fch4_oilgas_post[:, :, :] = ds_1['monthly_fch4_oilgasind_post'].values
            fch4_coal_prior[:, :, :] = ds_1['monthly_fch4_coal_prior'].values
            fch4_coal_post[:, :, :] = ds_1['monthly_fch4_coal_post'].values
            fch4_agri_prior[:, :, :] = ds_1['monthly_fch4_livestock_prior'].values
            fch4_agri_post[:, :, :] = ds_1['monthly_fch4_livestock_post'].values
            fch4_waste_prior[:, :, :] = ds_1['monthly_fch4_waste_prior'].values
            fch4_waste_post[:, :, :] = ds_1['monthly_fch4_waste_post'].values
            fch4_geol_prior[:, :, :] = ds_1['monthly_fch4_geological_prior'].values
            fch4_geol_post[:, :, :] = ds_1['monthly_fch4_geological_post'].values
            fch4_termite_prior[:, :, :] = ds_1['monthly_fch4_termites_prior'].values
            fch4_termite_post[:, :, :] = ds_1['monthly_fch4_termites_post'].values
            fch4_oce_prior[:, :, :] = ds_1['monthly_fch4_ocean_prior'].values
            fch4_oce_post[:, :, :] = ds_1['monthly_fch4_ocean_post'].values
            fch4_soils_prior[:, :, :] = ds_1['monthly_fch4_soils_prior'].values
            fch4_soils_post[:, :, :] = ds_1['monthly_fch4_soils_post'].values

            # --- date-time
            time[:] = ds_1["time_hours_since_1970"]

            # --- units
            cellarea.units = 'm2'
            flux_vars = [fch4_tot_prior, fch4_tot_post,
                         fch4_wet_prior, fch4_wet_post,
                         fch4_bb_prior, fch4_bb_post,
                         fch4_biofuel_prior, fch4_biofuel_post,
                         fch4_oilgas_prior, fch4_oilgas_post,
                         fch4_coal_prior, fch4_coal_post,
                         fch4_agri_prior, fch4_agri_post,
                         fch4_waste_prior, fch4_waste_post,
                         fch4_geol_prior, fch4_geol_post,
                         fch4_termite_prior, fch4_termite_post,
                         fch4_oce_prior, fch4_oce_post,
                         fch4_soils_prior, fch4_soils_post]

            for v in flux_vars:
                v.units = 'g-CH4/m2/month'

            print('\n SUCCESS writing total flux to >>> ', file_out)
            print('\t Check valid range')
            flux_pairs = {
                'tot': (fch4_tot_prior, fch4_tot_post),
                'wet': (fch4_wet_prior, fch4_wet_post),
                'bb': (fch4_bb_prior, fch4_bb_post),
                'biofuel': (fch4_biofuel_prior, fch4_biofuel_post),
                'oilgas': (fch4_oilgas_prior, fch4_oilgas_post),
                'coal': (fch4_coal_prior, fch4_coal_post),
                'agri': (fch4_agri_prior, fch4_agri_post),
                'waste': (fch4_waste_prior, fch4_waste_post),
                'geol': (fch4_geol_prior, fch4_geol_post),
                'termite': (fch4_termite_prior, fch4_termite_post),
                'oce': (fch4_oce_prior, fch4_oce_post),
                'soils': (fch4_soils_prior, fch4_soils_post),
            }

            for prior, post in flux_pairs.values():
                prior.valid_range = np.array((np.min(prior), np.max(prior)))
                post.valid_range = np.array((np.min(post), np.max(post)))

            def pr(label, arr):
                v1 = arr.valid_range
                print(f"      Check {label:<15} [{v1[0]:.5f}, {v1[1]:.5f}]")

            for name, (prior, post) in flux_pairs.items():
                pr(f'prior {name}', prior)
                pr(f'post  {name}', post)

            nc.close()

        # ===========================================================
        print('\n==========================================================================')
        print(f'\tApr flux from: {self.flx_apr_dir}')
        print(f'\tPst flux from: {dir()}')
        print(f'\tOutput *.nc file length: {self.nyear_nc} years for {self.years_nc} \n')

        for j1, invc in enumerate(self.invcases):
            flx_prior = get_prior_full_flx_nc()
            flx_post = get_post_flx_grd(self.unp, self.unx, invc)
            flx_join = get_joined_apr2post(flx_prior, flx_post)

            # - to write total
            print('\n>>>>>>>>>>>>')
            print('\t>>>>>> Write Total')
            ds = flx_join.sel(time=slice(f"{self.years_nc[0]}-01-01",
                                         f"{self.years_nc[1]}-12-31"))
            ds_ = check_total(ds)
            write_1nc_total(invc, ds_)

            # - to write categ
            print('\n>>>>>>>>>>>>')
            print('\t>>>>>> Write Category')
            ds_cat = get_post_categ(flx_prior, flx_join)
            ds = ds_cat.sel(time=slice(f"{self.years_nc[0]}-01-01",
                                       f"{self.years_nc[1]}-12-31"))
            ds_ = check_total(ds)
            write_1nc_cat(invc, ds_)
            # todo
            # cellarea = nc.createVariable('carea', 'f4', ('lat', 'lon'))
