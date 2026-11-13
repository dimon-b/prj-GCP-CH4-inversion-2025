# -*- coding: utf-8 -*-
"""
Finalized:  26/xx/xx
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import os
import netCDF4
import pandas as pd
from netCDF4 import date2num
from pathlib import Path
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors

import _set_case


class CheckHippo(_set_case.SetCase):

    def __init__(self):
        super().__init__()

        # === HIPPO
        self.f_hippo_obs = ('/S/data01/G5070/y0715/ccode/prj_apack/obs/o_orig/ObsPack/'
                            + 'obspack_ch4_1_GLOBALVIEWplus_v7.0_2024-10-29/'
                            + 'data/txt/ch4_hip_aircraft-insitu_59_allvalid.txt')
        self.campaigns = [
            ('HIPPO-1 Jan 2009', '2009-01-01', '2009-02-01', 'tab:blue'),
            ('HIPPO-2 Nov 2009', '2009-11-01', '2009-12-01', 'tab:orange'),
            ('HIPPO-3 Mar–Apr 2010', '2010-03-01', '2010-05-01', 'tab:green'),
            ('HIPPO-4 Jun–Jul 2011', '2011-06-01', '2011-08-01', 'tab:red'),
            ('HIPPO-5 Aug–Sep 2011', '2011-08-01', '2011-10-01', 'tab:purple'),
        ]

        # + 'data/txt/ch4_hip_aircraft-pfp_1_allvalid.txt')
        # self.f_colloc = 'ch4_hip_aircraft-pfp_1_allvalid_collocated.csv'
        self.f_colloc = 'ch4_hip_aircraft-insitu_59_allvalid_collocated.csv'

        # === ACTM
        self.val_model_h_cases = {
            'cyc_apr': {'dir': '../trout_MIROC/CH4_t42l67_CYC_apr_H_261103/',
                        'label': 'CYC APR',
                        'color': 'b',
                        'marker': 's',
                        },
            'cyc_pst': {'dir': '../trout_MIROC/CH4_t42l67_CYC_pst_H_261112/',
                        'label': 'CYC PST',
                        'color': 'b',
                        'marker': 'o',
                        },
            'inca_apr': {'dir': '../trout_MIROC/CH4_t42l67_INCA_apr_H_261104/',
                         'label': 'INCA APR',
                         'color': 'r',
                         'marker': 's',
                         },
            'inca_pst': {'dir': '../trout_MIROC/CH4_t42l67_CYC_pst_H_261112/',
                         'label': 'INCA PST',
                         'color': 'r',
                         'marker': '*',
                         }, }

        # ===
        update = False
        if update:
            self.df_hippo_obs = self.read_hippo_obs()
            self.get_hippo_collocations()

        # ===
        ### self.fig_name = f'val3_check_hippo_{case_id}_profiles.png'
        self.val3_hippo_profiles()
        self.fig_name = 'val3_check_hippo_profiles_campaigns.png'
        self.val3_hippo_profiles_compains()
        self.fig_name = 'val3_check_hippo_latlev.png'
        self.val3_hippo_latlev()

    # ===
    def val3_hippo_latlev(self):

        case_ids = ['cyc_apr', 'inca_apr', 'cyc_pst', 'inca_pst']
        cols = ['latitude', 'altitude', 'ch4_o', 'ch4_a']

        # - read files once
        rows = []

        for case_id in case_ids:
            case = self.val_model_h_cases[case_id]
            f_inp = os.path.join(case['dir'], self.f_colloc)
            df = pd.read_csv(f_inp, usecols=['dateTime', *cols])
            df['dateTime'] = pd.to_datetime(df['dateTime'], errors='coerce', utc=True).dt.tz_localize(None)
            df[cols] = df[cols].apply(pd.to_numeric, errors='coerce')
            df[cols] = df[cols].replace([np.inf, -np.inf], np.nan)
            df = df.dropna(subset=['dateTime', *cols])

            valid = (df['latitude'].between(-90, 90) & (df['altitude'] >= 0)
                     & (df['ch4_o'] > 0) & (df['ch4_a'] > 0))
            df = df.loc[valid]

            # - retain plotted campaign periods
            mask = pd.Series(False, index=df.index)
            for _, start, stop, _ in self.campaigns:
                mask |= (df['dateTime'] >= start) & (df['dateTime'] < stop)

            df = df.loc[mask].copy()
            df['alt_km'] = df['altitude'] / 1000.0
            df['diff'] = df['ch4_a'] - df['ch4_o']

            label = case['label']
            rows.append((label, df))
            print(f'{label}: {len(df)} valid HIPPO pairs')

        if all(df.empty for _, df in rows):
            print('No valid observations in the selected HIPPO campaigns')
            return

        # - common altitude range
        zmax = max(df['alt_km'].max() for _, df in rows if not df.empty)
        zmax = max(1.0, np.ceil(zmax))

        # - common discrete colour scale
        bounds = np.arange(-60, 61, 10)
        cmap = plt.get_cmap('PuOr_r', len(bounds) - 1)
        norm = mcolors.BoundaryNorm(bounds, cmap.N)

        # - create subplots
        ncols = len(self.campaigns)
        fig, axes = plt.subplots(len(rows), ncols, figsize=(13, 13), sharex=True, sharey=True, squeeze=False)
        fig.subplots_adjust(left=0.07, right=0.94, bottom=0.15, top=0.92, wspace=0.12, hspace=0.24)

        # - loop model cases and campaigns
        for i, (label, df) in enumerate(rows):
            for j, (title, start, stop, color) in enumerate(self.campaigns):
                ax = axes[i, j]
                df_c = df.loc[(df['dateTime'] >= start) & (df['dateTime'] < stop)]

                ax.set_title(f'{chr(97 + i * ncols + j)})', fontsize=12)
                ax.set_xlim(-90, 90)
                ax.set_ylim(0, zmax)
                ax.set_xticks(np.arange(-90, 91, 30))
                ax.tick_params(axis='both', labelsize=10)
                ax.set_axisbelow(True)
                ax.grid(linestyle=':', linewidth=0.5)

                if i == 0:
                    ax.text(0.5, 1.16, title, transform=ax.transAxes, ha='center', va='bottom', fontsize=12,
                            color=color)

                if j == ncols - 1:
                    ax.text(1.08, 0.5, label, transform=ax.transAxes, rotation=270, ha='center', va='center',
                            fontsize=12)

                if df_c.empty:
                    ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center', va='center')
                    continue

                # - scatter of model minus observations
                ax.scatter(df_c['latitude'], df_c['alt_km'], c=df_c['diff'], cmap=cmap, norm=norm,
                           marker='s', s=9, edgecolors='none', rasterized=True)

                # - statistics for plotted observations
                bias = df_c['diff'].mean()
                sd = df_c['diff'].std(ddof=0)
                ax.text(0.03, 0.03, f'Bias = {bias:.2f}\nSD = {sd:.2f}', transform=ax.transAxes,
                        ha='left', va='bottom', fontsize=10,
                        bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))

                print(f'{label:12s} {title:25s} N = {len(df_c):7d}, Bias = {bias:.2f}, SD = {sd:.2f}')

        # - common axis labels
        fig.text(0.5, 0.11, 'Latitude', ha='center', fontsize=13)
        fig.text(0.025, 0.54, 'Altitude (km)', va='center', rotation=90, fontsize=13)

        # - common colour bar
        cax = fig.add_axes([0.24, 0.055, 0.52, 0.018])
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        cb = fig.colorbar(sm, cax=cax, orientation='horizontal', ticks=bounds, extend='both')
        cb.set_label('Model − observations CH₄ (ppb)', fontsize=12)
        cb.ax.tick_params(labelsize=10)

        # - save plot
        os.makedirs(self.plt_dir, exist_ok=True)
        f_out = os.path.join(self.plt_dir, self.fig_name)

        try:
            fig.savefig(f_out, format='png', dpi=600, bbox_inches='tight')
            # plt.show()
        finally:
            plt.close(fig)

        print('Plot saved to', f_out)

    # ===
    def val3_hippo_profiles_compains(self):

        cols = ['dateTime', 'altitude', 'ch4_a']
        dz = 200.0

        styles = {'apr': '-', 'pst': '--'}

        def read_profile(case):
            path = os.path.join(case['dir'], self.f_colloc)
            df = pd.read_csv(path, usecols=cols)
            df['dateTime'] = pd.to_datetime(df['dateTime'], errors='coerce', utc=True, ).dt.tz_localize(None)
            df[['altitude', 'ch4_a']] = df[['altitude', 'ch4_a']].apply(pd.to_numeric, errors='coerce', )
            return df.replace([np.inf, -np.inf], np.nan).dropna(subset=cols)

        # - read the four collocated CSV files
        profiles = {}
        for case_id in ('cyc', 'inca'):
            for stage in ('apr', 'pst'):
                key = f'{case_id}_{stage}'
                profiles[key] = read_profile(self.val_model_h_cases[key])

        # - common altitude bins for both panels
        altitude = pd.concat([df['altitude'] for df in profiles.values()])
        if altitude.empty:
            print('HIPPO: no valid profile data')
            return

        lower = np.floor(altitude.min() / dz) * dz
        upper = max(lower + dz, np.ceil(altitude.max() / dz) * dz)
        bins = np.arange(lower, upper + dz, dz)
        centers = (bins[:-1] + bins[1:]) / 2.0

        fig, axes = plt.subplots(1, 2, figsize=(8, 6), sharex=True, sharey=True, )

        # - one panel per model, one APR/PST pair per campaign
        for ax, case_id in zip(axes, ('cyc', 'inca')):
            for title, start, stop, color in self.campaigns:
                for stage, linestyle in styles.items():
                    key = f'{case_id}_{stage}'
                    case = self.val_model_h_cases[key]
                    df = profiles[key]
                    df_c = df.loc[df['dateTime'].ge(start) & df['dateTime'].lt(stop)]

                    if df_c.empty:
                        continue

                    alt_bin = pd.cut(df_c['altitude'], bins, labels=False, include_lowest=True, )
                    sd = df_c.groupby(alt_bin)['ch4_a'].std(ddof=0)
                    sd = sd.reindex(range(len(centers)))

                    ax.plot(sd.to_numpy(), centers, color=color, linestyle=linestyle, marker=case['marker'],
                            markersize=1, lw=0.7)
                    print(f"{case['label']} | {title} | N = {len(df_c)}")

            ax.set_title(case_id.upper())
            ax.set_xlabel('Within-bin CH₄ SD (ppb)')
            ax.grid(linestyle=':')

            if not ax.lines:
                ax.text(0.5, 0.5, 'No data', transform=ax.transAxes,
                        ha='center', va='center', )

        axes[0].set_ylabel('Altitude (m)')

        campaign_handles = [Line2D([0], [0], color=color, lw=2, label=title)
                            for title, _, _, color in self.campaigns]
        stage_handles = [Line2D([0], [0], color='k', linestyle=linestyle, lw=1.5, label=stage.upper())
                         for stage, linestyle in styles.items()]
        fig.legend(handles=campaign_handles + stage_handles,
                   loc='lower center', ncol=4, frameon=False, fontsize=8)
        fig.tight_layout(rect=(0, 0.13, 1, 1))

        f_out = os.path.join(self.plt_dir, self.fig_name)
        try:
            fig.savefig(f_out, dpi=600, bbox_inches='tight')
        finally:
            plt.close(fig)

        print('Plot saved to', f_out)

    # ===
    def val3_hippo_profiles(self):
        cols = ['altitude', 'ch4_a']
        dz = 100.0

        def read_profile(case):
            path = os.path.join(case['dir'], self.f_colloc)
            df = pd.read_csv(path, usecols=cols)
            df = df.apply(pd.to_numeric, errors='coerce')
            return df.replace([np.inf, -np.inf], np.nan).dropna()

        def bin_profile(df, bins):
            alt_bin = pd.cut(df['altitude'], bins, labels=False, include_lowest=True)
            sd = df.groupby(alt_bin)['ch4_a'].std(ddof=0)
            centers = (bins[:-1] + bins[1:]) / 2.0
            return centers[sd.index.to_numpy(dtype=int)], sd.to_numpy()

        # - pair each APR case with its PST case
        for apr_id, apr_case in self.val_model_h_cases.items():
            if not apr_id.endswith('_apr'):
                continue

            case_id = apr_id.removesuffix('_apr')
            pst_case = self.val_model_h_cases[f'{case_id}_pst']

            df_apr = read_profile(apr_case)
            df_pst = read_profile(pst_case)

            if df_apr.empty or df_pst.empty:
                print(f'{case_id}: missing valid APR or PST profile data')
                continue

            zmin = min(df_apr['altitude'].min(), df_pst['altitude'].min())
            zmax = max(df_apr['altitude'].max(), df_pst['altitude'].max())
            lower = np.floor(zmin / dz) * dz
            upper = max(lower + dz, np.ceil(zmax / dz) * dz)
            bins = np.arange(lower, upper + dz, dz)

            z_apr, sd_apr = bin_profile(df_apr, bins)
            z_pst, sd_pst = bin_profile(df_pst, bins)

            fig, ax = plt.subplots(figsize=(6, 6))
            ax.plot(sd_apr, z_apr, color='r', markersize=4, lw=1.0,
                    marker=apr_case['marker'], label=apr_case['label'], )
            ax.plot(sd_pst, z_pst, color='b', linestyle='--', markersize=4, lw=1.0,
                    marker=pst_case['marker'], label=pst_case['label'], )

            ax.set_xlabel('Within-bin CH₄ SD (ppb)')
            ax.set_ylabel('Altitude (m)')
            ax.set_title(f'HIPPO CH₄ profiles | {case_id.upper()}')
            ax.grid(linestyle=':')
            ax.legend(loc='best')
            fig.tight_layout()

            self.fig_name = f'val3_check_hippo_{case_id}_profiles.png'
            f_out = os.path.join(self.plt_dir, self.fig_name)
            try:
                fig.savefig(f_out, dpi=600, bbox_inches='tight')
            finally:
                plt.close(fig)

            print(f'{case_id}: APR N = {len(df_apr)}, PST N = {len(df_pst)}')
            print('Plot saved to', f_out)

    # ===
    def read_hippo_obs(self):
        cols = ['year', 'month', 'day', 'hour', 'minute',
                'latitude', 'longitude', 'altitude', 'value', ]

        df = pd.read_csv(self.f_hippo_obs, sep=r'\s+', comment='#',
                         usecols=cols, na_values=[-1e34], )
        df = df.loc[:, cols].copy()

        # - convert CH4 mole fraction from mol/mol to ppb
        df['value'] = df['value'] * 1e9

        print(f'HIPPO ObsPack: read {len(df)} observations from {self.f_hippo_obs}')
        return df

    # ===
    def get_hippo_collocations(self):
        f_save = self.f_colloc
        time_cols = ['year', 'month', 'day', 'hour', 'minute']
        cols = ['dateTime', 'latitude', 'longitude', 'altitude', 'ch4_o']

        # - sample model CH4 at the observation locations
        def get_model_df(df_o, actm_dir, z_dir):
            inq = self.invq[0]
            frames = []

            for year, df_yr in df_o.groupby(df_o['dateTime'].dt.year):
                df_yr = df_yr.copy()
                ch4_path = os.path.join(actm_dir, str(year), 'ch4_inv1.nc')
                z_path = os.path.join(z_dir, str(year), 'z.nc')

                times = pd.date_range(f'{year}-01-01', f'{year + 1}-01-01',
                                      freq='h', inclusive='left', )
                hour_idx = times.get_indexer(df_yr['dateTime'], method='nearest')
                values = np.full(len(df_yr), np.nan)

                with xr.open_dataset(ch4_path, decode_times=False) as ds, \
                        xr.open_dataset(z_path, decode_times=False) as ds_z:
                    ch4 = ds[inq]
                    z = ds_z['Z']

                    for path, field in ((ch4_path, ch4), (z_path, z)):
                        if field.sizes['time'] != len(times):
                            raise ValueError(f"{path}: expected {len(times)} hourly records, "
                                             f"found {field.sizes['time']}")

                    for dim in ('x', 'y', 'level'):
                        if not ch4.get_index(dim).equals(z.get_index(dim)):
                            raise ValueError(f'{ch4_path} and {z_path}: different {dim} coordinates')

                    x = ch4['x'].to_numpy()
                    y_idx = ch4.get_index('y').get_indexer(df_yr['latitude'], method='nearest', )
                    lon = df_yr['longitude'].to_numpy()
                    alt = df_yr['altitude'].to_numpy()

                    # - convert mol/mol to ppb when declared in the NetCDF
                    units = str(ch4.attrs.get('units', '')).lower().replace(' ', '')
                    ch4_factor = 1e9 if units in {'mol/mol', 'molmol-1', 'molmol^-1'} else 1.0

                    # - sample observations grouped by nearest model hour
                    for hour, pos in df_yr.groupby(hour_idx).indices.items():
                        dlon = (x[None, :] - lon[pos, None] + 180) % 360 - 180
                        ix = xr.DataArray(np.abs(dlon).argmin(axis=1), dims='point')
                        iy = xr.DataArray(y_idx[pos], dims='point')

                        z_profile = (
                            z.isel(time=hour, y=iy, x=ix)
                            .transpose('point', 'level')
                            .to_numpy()
                        )
                        dz = np.abs(z_profile - alt[pos, None])
                        dz = np.where(np.isfinite(dz), dz, np.inf)
                        valid_z = np.isfinite(dz).any(axis=1)
                        iz = xr.DataArray(dz.argmin(axis=1), dims='point')

                        sampled = ch4.isel(time=hour, y=iy, x=ix, level=iz).to_numpy()
                        values[pos[valid_z]] = sampled[valid_z] * ch4_factor

                df_yr['ch4_a'] = values
                frames.append(df_yr)

            if not frames:
                return df_o.assign(ch4_a=np.nan)

            return pd.concat(frames, ignore_index=True)

        # - prepare original HIPPO observations, CH4 in ppb
        obs = self.df_hippo_obs.rename(columns={'value': 'ch4_o'}).copy()
        obs['dateTime'] = pd.to_datetime(obs[time_cols], errors='coerce')
        obs = obs[cols].replace([np.inf, -np.inf], np.nan).dropna()

        valid = (obs['latitude'].between(-90, 90)
                 & obs['longitude'].between(-180, 360)
                 & obs['altitude'].ge(0)
                 & obs['ch4_o'].gt(0))
        obs = obs.loc[valid]

        # - collocate and save one CSV per case
        for case in self.val_model_h_cases.values():
            actm_dir = case['dir']
            z_dir = case.get('z_dir', actm_dir)

            df = get_model_df(obs, actm_dir, z_dir)
            df = df.loc[np.isfinite(df['ch4_a']) & df['ch4_a'].gt(0)]

            f_out = os.path.join(actm_dir, f_save)
            df.to_csv(f_out, index=False)
            print(f"{case['label']}: saved {len(df)} / {len(obs)} HIPPO pairs to {f_out}")

    #
    # --- OLD ---
    #
    # ---
    def val2_hippo_scatter_post_(self):

        paths = ['../trout_MIROC/CH4_t42l67_CYC_H_post_260323/',
                 '../trout_MIROC/CH4_t42l67_INCA_H_post_260324/'
                 ]
        labels = ['Post CYC', 'Post INCA']
        f_name = 'ch4_hip_aircraft-insitu_59_allvalid_out.csv'

        # - campaign periods, exclusive stop dates
        campaigns = [['HIPPO-1 Jan 2009', '2009-01-01', '2009-02-01'],
                     ['HIPPO-2 Nov 2009', '2009-11-01', '2009-12-01'],
                     ['HIPPO-3 Mar–Apr 2010', '2010-03-01', '2010-05-01'],
                     ['HIPPO-4 Jun–Jul 2011', '2011-06-01', '2011-08-01'],
                     ['HIPPO-5 Aug–Sep 2011', '2011-08-01', '2011-10-01']
                     ]

        # - read files once
        dfs = []

        for path, label in zip(paths, labels):
            f_inp = os.path.join(path, f_name)
            df = pd.read_csv(f_inp)
            df['dateTime'] = pd.to_datetime(df['dateTime'], errors='coerce', utc=True).dt.tz_localize(None)

            cols = ['latitude', 'altitude', 'ch4_o', 'ch4_a']
            for col in cols:
                df[col] = pd.to_numeric(df[col], errors='coerce')

            df[cols] = df[cols].replace([np.inf, -np.inf], np.nan)
            df = df.dropna(subset=['dateTime'] + cols).copy()
            df = df[
                df['latitude'].between(-90, 90) & (df['altitude'] >= 0) & (df['ch4_o'] > 0) & (df['ch4_a'] > 0)].copy()

            # - altitude in km, concentrations already in ppb
            df['alt_km'] = df['altitude'] / 1000.0
            df['diff'] = df['ch4_a'] - df['ch4_o']

            # - retain plotted campaign periods
            mask = pd.Series(False, index=df.index)
            for title, start, stop in campaigns:
                mask = mask | ((df['dateTime'] >= start) & (df['dateTime'] < stop))

            df = df[mask].copy()
            dfs.append(df)
            print(label, ' plotted observations : ', len(df))

        if all(df.empty for df in dfs):
            print('No valid observations in the selected campaign periods')
            return

        # - common altitude range
        zmax = max(df['alt_km'].max() for df in dfs if not df.empty)
        zmax = max(1.0, np.ceil(zmax))

        # - common discrete colour scale
        bounds = np.arange(-60, 61, 10)
        cmap = plt.get_cmap('PuOr_r', len(bounds) - 1)
        norm = mcolors.BoundaryNorm(bounds, cmap.N)

        # - create subplots
        fig, axes = plt.subplots(2, 5, figsize=(13, 6), sharex=True, sharey=True)
        fig.subplots_adjust(left=0.07, right=0.94, bottom=0.20, top=0.87, wspace=0.12, hspace=0.24)

        # - loop model cases and campaigns
        for i, (df, label) in enumerate(zip(dfs, labels)):
            for j, (title, start, stop) in enumerate(campaigns):
                ax = axes[i, j]
                df_y = df[(df['dateTime'] >= start) & (df['dateTime'] < stop)].copy()

                ax.set_title(f'{chr(97 + i * 5 + j)})', fontsize=12)
                ax.set_xlim([-90, 90])
                ax.set_ylim([0, zmax])
                ax.set_yticks(np.arange(0, 16, 3))
                ax.set_xticks(np.arange(-90, 91, 30))
                ax.tick_params(axis='both', labelsize=10)
                ax.set_axisbelow(True)
                ax.grid(linestyle=':', linewidth=0.5)

                if i == 0:
                    ax.text(0.5, 1.16, title, transform=ax.transAxes, ha='center', va='bottom', fontsize=12)

                if j == 4:
                    ax.text(1.08, 0.5, label, transform=ax.transAxes, rotation=270, ha='center', va='center',
                            fontsize=12)

                if df_y.empty:
                    ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center', va='center')
                    continue

                # - scatter of model minus observations
                ax.scatter(df_y['latitude'], df_y['alt_km'], c=df_y['diff'], cmap=cmap, norm=norm, marker='s', s=9,
                           edgecolors='none', rasterized=True)

                # - statistics for plotted observations
                bias = df_y['diff'].mean()
                sd = df_y['diff'].std(ddof=0)
                ax.text(0.03, 0.03, f'Bias = {bias:.2f}\nSD = {sd:.2f}', transform=ax.transAxes, ha='left', va='bottom',
                        fontsize=10, bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))
                print(f'{label:12s} {title:25s} N = {len(df_y):7d}, Bias = {bias:.2f}, SD = {sd:.2f}')

        # - common axis labels
        # fig.text(0.5, 0.13, 'Latitude', ha='center', fontsize=13)
        fig.text(0.025, 0.54, 'Altitude (km)', va='center', rotation=90, fontsize=13)

        # - common colour bar
        cax = fig.add_axes([0.24, 0.095, 0.52, 0.025])
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        cb = fig.colorbar(sm, cax=cax, orientation='horizontal', ticks=bounds, extend='both')
        cb.set_label('Model − observations CH₄ (ppb)', fontsize=12)
        cb.ax.tick_params(labelsize=10)

        # - save plot
        f_out = os.path.join(self.plt_dir, 'val2_hippo_scatter_post.png')
        fig.savefig(f_out, format='png', dpi=600, bbox_inches='tight')
        print('\nPlot saved to : ', f_out)

        plt.show()
        plt.close(fig)

    # ---
    def val2_hippo_rmse_(self):
        import os
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D

        paths = [
            '../trout_MIROC/CH4_t42l67_CYC_H_post_260323/',
            '../trout_MIROC/CH4_t42l67_INCA_H_post_260324/'
        ]
        labels = ['Post CYC', 'Post INCA']
        styles = ['-', '--']
        colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
        f_name = 'ch4_hip_aircraft-insitu_59_allvalid_out.csv'

        # - campaign periods, exclusive stop dates
        campaigns = [
            ['HIPPO-1 Jan 2009', '2009-01-01', '2009-02-01'],
            ['HIPPO-2 Nov 2009', '2009-11-01', '2009-12-01'],
            ['HIPPO-3 Mar–Apr 2010', '2010-03-01', '2010-05-01'],
            ['HIPPO-4 Jun–Jul 2011', '2011-06-01', '2011-08-01'],
            ['HIPPO-5 Aug–Sep 2011', '2011-08-01', '2011-10-01']
        ]

        # - altitude bins, 500 m
        bins = np.arange(0, 15001, 500)
        centers = (bins[:-1] + bins[1:]) / 2000.0

        # --- RMSE profile
        def bin_profile(df):
            df = df.copy()
            df['alt_bin'] = pd.cut(df['altitude'], bins, labels=False, include_lowest=True)
            df['diff2'] = df['diff'] ** 2
            mean_diff2 = df.groupby('alt_bin')['diff2'].mean()
            rmse = np.sqrt(mean_diff2.reindex(range(len(centers))))
            return centers, rmse.to_numpy()

        # - read files once
        dfs = []

        for path, label in zip(paths, labels):
            df = pd.read_csv(os.path.join(path, f_name))
            df['dateTime'] = pd.to_datetime(df['dateTime'], errors='coerce', utc=True).dt.tz_localize(None)

            cols = ['latitude', 'altitude', 'ch4_o', 'ch4_a']
            for col in cols:
                df[col] = pd.to_numeric(df[col], errors='coerce')

            df[cols] = df[cols].replace([np.inf, -np.inf], np.nan)
            df = df.dropna(subset=['dateTime'] + cols).copy()
            df = df[df['latitude'].between(-90, 90) & df['altitude'].between(0, 15000) & (df['ch4_o'] > 0) & (
                    df['ch4_a'] > 0)].copy()

            # - model minus observations, ppb
            df['diff'] = df['ch4_a'] - df['ch4_o']
            dfs.append(df)

        # - one panel
        fig, ax = plt.subplots(figsize=(8, 6))
        campaign_handles = []
        plotted_cases = set()

        # - campaign colour; model line style
        for (title, start, stop), color in zip(campaigns, colors):
            plotted = False

            for i, (df, label, style) in enumerate(zip(dfs, labels, styles)):
                df_y = df[(df['dateTime'] >= start) & (df['dateTime'] < stop)].copy()

                if df_y.empty:
                    continue

                z, rmse = bin_profile(df_y)
                ax.plot(rmse, z, color=color, linestyle=style, linewidth=2)
                plotted = True
                plotted_cases.add(i)

                rmse_all = np.sqrt(np.mean(df_y['diff'] ** 2))
                print(f'{label:12s} {title:25s} N = {len(df_y):7d}, RMSE = {rmse_all:.2f} ppb')

            if plotted:
                campaign_handles.append(Line2D([0], [0], color=color, linewidth=2, label=title))

        if not campaign_handles:
            print('No valid observations in the selected campaign periods')
            plt.close(fig)
            return

        # - axes
        ax.set_xlabel('RMSE of model − observations CH₄ (ppb)', fontsize=12)
        ax.set_ylabel('Altitude (km)', fontsize=12)
        ax.set_xlim(left=0)
        ax.set_ylim([0, 15])
        ax.set_yticks(np.arange(0, 16, 3))
        ax.tick_params(axis='both', labelsize=11)
        ax.grid(linestyle=':', linewidth=0.5)

        # - legend in lower right
        model_handles = [Line2D([0], [0], color='black', linestyle=styles[i], linewidth=2, label=labels[i]) for i in
                         sorted(plotted_cases)]
        ax.legend(handles=campaign_handles + model_handles, loc='lower right', ncol=1, fontsize=10, frameon=False)

        # - save plot
        # fig.subplots_adjust(left=0.13, right=0.97, top=0.97, bottom=0.30)
        f_out = os.path.join(self.plt_dir, 'val2_hippo_rmse_profiles_post.png')
        fig.savefig(f_out, format='png', dpi=600, bbox_inches='tight')
        print('\nPlot saved to : ', f_out)

        plt.show()
        plt.close(fig)

    # ---
    def val2_hippo_profiles_yr_(self):
        dirs = [
            ['../trout_MIROC/CH4_t42l67_CYC_H_251014/', '../trout_MIROC/CH4_t42l67_CYC_H_post_260323/'],
            ['../trout_MIROC/CH4_t42l67_INCA_H_260306/', '../trout_MIROC/CH4_t42l67_INCA_H_post_260324/']
        ]

        f_name = 'ch4_hip_aircraft-insitu_59_allvalid_out.csv'

        # --- read reference
        df_ref = pd.read_csv(os.path.join(dirs[0][0], f_name))
        print(df_ref)
        # exit()
        df_ref['dateTime'] = pd.to_datetime(df_ref['dateTime'])

        bins = np.arange(df_ref['altitude'].min(), df_ref['altitude'].max() + 500, 500)

        def bin_profile(df, var):
            df = df.copy()
            df['alt_bin'] = pd.cut(df['altitude'], bins)
            grouped = df.groupby('alt_bin', observed=True)[var]

            mean = grouped.mean()
            rmse = np.sqrt(grouped.apply(lambda x: np.mean((x - x.mean()) ** 2)))

            centers = [interval.mid for interval in mean.index]
            return centers, mean.values, rmse.values

        # --- get years
        years = sorted(df_ref['dateTime'].dt.year.unique())

        # --- create subplots
        n = len(years)
        fig, axes = plt.subplots(1, n, figsize=(5 * n, 6), sharey=True)

        if n == 1:
            axes = [axes]

        # --- loop years
        for ax, year in zip(axes, years):

            ch4_o_plotted = False

            for d_apr, d_pst in dirs:

                df_apr = pd.read_csv(os.path.join(d_apr, f_name))
                df_pst = pd.read_csv(os.path.join(d_pst, f_name))

                df_apr['dateTime'] = pd.to_datetime(df_apr['dateTime'])
                df_pst['dateTime'] = pd.to_datetime(df_pst['dateTime'])

                # --- filter by year
                df_apr_y = df_apr[df_apr['dateTime'].dt.year == year]
                df_pst_y = df_pst[df_pst['dateTime'].dt.year == year]

                if df_apr_y.empty:
                    continue

                name_apr = os.path.basename(os.path.normpath(d_apr)).split('_')[2]
                name_pst = name_apr + '_post'

                # --- OBS once per subplot
                if not ch4_o_plotted:
                    z_o, ch4_o, err_o = bin_profile(df_apr_y, 'ch4_o')
                    ax.errorbar(ch4_o, z_o, xerr=err_o,
                                color='black', linewidth=2, label='OBS', fmt='-')
                    ch4_o_plotted = True

                # --- model
                z_apr, ch4_apr, err_apr = bin_profile(df_apr_y, 'ch4_a')
                z_pst, ch4_pst, err_pst = bin_profile(df_pst_y, 'ch4_a')

                # ax.errorbar(ch4_apr, z_apr, xerr=err_apr, linestyle='-', label=name_apr, fmt='o', alpha=0.7)
                ax.errorbar(ch4_pst, z_pst, xerr=err_pst, linestyle='--', label=name_pst, fmt='s', alpha=0.7)

            ax.set_title(f'{year}')
            ax.grid()

        axes[0].set_ylabel('Altitude (m)')
        for ax in axes:
            ax.set_xlabel('CH₄ (ppb)')

        axes[-1].legend(loc='best')

        plt.suptitle('CH₄ profiles by year (binned, 500 m, RMSE)')
        plt.tight_layout()
        plt.savefig(self.plt_dir + 'val2_hippo_profiles_yr.png', format='png', dpi=600, bbox_inches='tight')
        plt.show()

    # ---
    def val2_hippo_profiles_(self):

        # ---
        def valis_comp(apr_csv_, pst_csv_, id_):

            # ---
            def read_comp_csv():
                def get_1file(ifile_: str) -> pd.DataFrame:
                    df = pd.read_csv(ifile_)  # , sep=r'\s+', skiprows=0)
                    return df

                ifile = apr_csv_
                df_apr = get_1file(ifile)
                len_apr = len(df_apr)
                # print(df_apr.head())

                ifile = pst_csv_
                df_pst = get_1file(ifile)
                len_pst = len(df_pst)
                # print(df_pst.head())

                return df_apr, df_pst

            # ---
            def check_1_comp(apr, pst):
                # use common altitude range for both datasets
                zmin = min(apr['altitude(m)'].min(), pst['altitude(m)'].min())
                zmax = max(apr['altitude(m)'].max(), pst['altitude(m)'].max())

                bins = np.arange(zmin, zmax + 100, 100)

                def compute_rmse(df):
                    df = df.copy()
                    df['alt_bin'] = pd.cut(df['altitude(m)'], bins)
                    grouped = df.groupby('alt_bin', observed=True)['ch4_a']

                    # RMSE around bin mean (i.e. variability)
                    rmse = np.sqrt(grouped.apply(lambda x: np.mean((x - x.mean()) ** 2)))

                    # bin centers
                    centers = np.array([interval.mid for interval in rmse.index])

                    return centers, rmse.values

                # compute for both
                z_apr, rmse_apr = compute_rmse(apr)
                z_pst, rmse_pst = compute_rmse(pst)

                # plot
                plt.plot(rmse_apr, z_apr, label='APR', marker='o')
                plt.plot(rmse_pst, z_pst, label='PST', marker='s')

                plt.xlabel('RMSE CH₄ (ppb)')
                plt.ylabel('Altitude (m)')
                plt.title('RMSE of CH₄ vs Altitude')
                plt.grid()
                plt.legend()
                plt.show()
                # exit()

            # ===========================================================
            print(f'\n\t\tRun read_comp_txt for: {apr_csv_, pst_csv_, id_}')
            for j1, invc in enumerate(self.invcases[:]):
                a_conc, p_conc = read_comp_csv()
                check_1_comp(a_conc, p_conc)

        # ---
        def get_csv_files(dir_):
            return [os.path.join(dir_, f) for f in os.listdir(dir_) if f.lower().endswith('.csv')]

        # ===========================================================
        self.invcases = ['inv1']
        dirs = [
            ['../trout_MIROC/CH4_t42l67_CYC_H_251014/', '../trout_MIROC/CH4_t42l67_CYC_H_post_260323/'],
            ['../trout_MIROC/CH4_t42l67_INCA_H_260306/', '../trout_MIROC/CH4_t42l67_INCA_H_post_260324/']
        ]

        for ddir in dirs:
            print('')
            print("\n List of apr csv files :")
            apr_csv_files = get_csv_files(ddir[0])
            apr_csv_files = [f for f in apr_csv_files if 'station' in f]
            print(",\n".join(f'\t\t{f}' for f in apr_csv_files))

            print("\n List of pst csv files :")
            pst_csv_files = get_csv_files(ddir[1])
            pst_csv_files = [f for f in pst_csv_files if 'station' in f]
            print(",\n".join(f'\t\t{f}' for f in pst_csv_files))

            # --- safer pairing (recommended)
            apr_csv_files.sort()
            pst_csv_files.sort()

            for apr_csv, pst_csv in zip(apr_csv_files, pst_csv_files):
                filename = os.path.basename(str(apr_csv))
                cmp_id = filename.split("_")[0]

                valis_comp(apr_csv, pst_csv, cmp_id)
