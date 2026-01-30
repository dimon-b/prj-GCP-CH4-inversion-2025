# -*- coding: utf-8 -*-
"""
Created/Corrected:    25/09/09 09:09
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset, stringtochar

import _set_case


class ObsFilesTxt(_set_case.SetCase):
    def __init__(self):
        super().__init__()

        print('Read Obs from dir: ')
        print('\t\t', self.flx_inp_dir)

        # ---
        # self.get_obs_files_obspack()
        self.get_obs_files_wdcgg()

    # --- get_obs_files_wdcgg
    def get_obs_files_wdcgg(self):
        # - get_obs_sites
        def get_obs_sites():
            df = pd.read_csv(self.sites_f, sep='\\s+', skiprows=2, header=None, usecols=range(6))
            df.columns = ['idx', 'sno', 'Sitename', 'Lat_o', 'Lon_o', 'Alt_o']
            return df

        # --- find_files
        def find_files(sitename):
            site_code = sitename.split('_')[1]
            search_order = ['hourly', 'daily', 'event']
            all_filtered_files = []

            for subdir in search_order:
                dir_path = os.path.join(self.wdcgg_dir, subdir, '**')
                full_pattern = os.path.join(dir_path, f'*ch4_{site_code}*')
                files = glob.glob(full_pattern, recursive=True)

                if not files:
                    continue  # No files found, try next directory

                # Filter out files containing '_met'
                filtered_files = [f for f in files if '_met' not in os.path.basename(f)]

                if not filtered_files:
                    print(f"For {sitename}: Found files in {subdir}/ but all contain '_met'")
                    continue

                # Add to total list
                all_filtered_files.extend(filtered_files)

                # Print found files for this directory
                clean_files = [f.replace(self.wdcgg_dir, '').lstrip('/') for f in filtered_files]
                # print(f"For {sitename} found in {subdir}/: {len(filtered_files)} files", clean_files)

            if all_filtered_files:
                print(f"For {sitename}: Total {len(all_filtered_files)} files found across all directories")
            else:
                print(f"For {sitename}: No suitable files found in any directory")

            return all_filtered_files

        # ---
        def find_file_with_largest_time_span(files):
            """Find best file with rules"""

            hourly_daily_files = []
            event_files = []

            for file in files:
                if '_met' in os.path.basename(file):
                    continue

                try:
                    # NORMALIZE PATH for consistent checking
                    normalized_path = os.path.normpath(file)

                    with open(file, 'r', encoding='utf-8', errors='ignore') as f:
                        lines = f.readlines()

                    # Find data rows
                    data_lines = []
                    for line in lines:
                        if line.strip() and not line.startswith('#'):
                            data_lines.append(line.strip())

                    if len(data_lines) < 2:
                        continue

                    first_row = data_lines[0].split()
                    last_row = data_lines[-1].split()

                    year1 = int(first_row[1])
                    year2 = int(last_row[1])

                    # Calculate time span
                    month1 = int(first_row[2]) if len(first_row) > 2 else 1
                    month2 = int(last_row[2]) if len(last_row) > 2 else 12
                    day1 = int(first_row[3]) if len(first_row) > 3 else 1
                    day2 = int(last_row[3]) if len(last_row) > 3 else 1
                    diff_days = (year2 - year1) * 365 + (month2 - month1) * 30 + (day2 - day1)

                    # Check condition
                    meets_condition = (year1 < 1999 and year2 >= 2024)

                    file_info = {
                        'file': file,
                        'year1': year1,
                        'year2': year2,
                        'diff_days': diff_days,
                        'meets_condition': meets_condition
                    }

                    # CORRECT: Check using normalized path
                    if f'{os.path.sep}hourly{os.path.sep}' in normalized_path or f'{os.path.sep}daily{os.path.sep}' in normalized_path:
                        hourly_daily_files.append(file_info)
                        # print(f"Added as hourly/daily: {os.path.basename(file)}")
                    elif f'{os.path.sep}event{os.path.sep}' in normalized_path:
                        event_files.append(file_info)
                        # print(f"Added as event: {os.path.basename(file)}")
                    else:
                        print(f"WARNING: Could not categorize {file}")

                except Exception as e:
                    print(f"Error processing {os.path.basename(file)}: {e}")
                    continue

            print(f"\tFound {len(hourly_daily_files)} hourly/daily files, {len(event_files)} event files")

            # Rest of the selection logic remains the same...
            # RULE 1: Check if any hourly/daily meets condition
            hourly_daily_with_condition = [f for f in hourly_daily_files if f['meets_condition']]

            if hourly_daily_with_condition:
                best = max(hourly_daily_with_condition, key=lambda x: x['diff_days'])
                print(f"\tSELECT: Hourly/daily file (meets condition)")
                print(f"\tFile: {os.path.basename(best['file'])}")
                print(f"\tTime span: {best['diff_days']} days ({best['year1']}-{best['year2']})")
                return best['file']

            # RULE 2: Select best event file
            if event_files:
                best = max(event_files, key=lambda x: x['diff_days'])
                print(f"\tSELECT: Event file")
                print(f"\tFile: {os.path.basename(best['file'])}")
                print(f"\tTime span: {best['diff_days']} days ({best['year1']}-{best['year2']})")
                return best['file']

            # RULE 3: Fallback
            if hourly_daily_files:
                best = max(hourly_daily_files, key=lambda x: x['diff_days'])
                print(f"\tSELECT: Hourly/daily file (fallback)")
                print(f"\tFile: {os.path.basename(best['file'])}")
                print(f"\tTime span: {best['diff_days']} days ({best['year1']}-{best['year2']})")
                return best['file']

            print("\nNo suitable files found")
            return None

        # ---
        def read_wdcgg_file(filepath):
            # Read the first few lines to find header_lines
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                first_line = f.readline().strip()

            # Check if first line contains header_lines info
            if first_line.startswith('# header_lines :'):
                header_lines = int(first_line.split(':')[1].strip())
            else:
                exit('No header_lines found in file')

            # Read the data
            try:
                # Skip header lines and read the data
                df_r = pd.read_csv(filepath, skiprows=header_lines - 1, delimiter='\\s+')
                cols = df_r.columns.tolist()
                df_r.columns = df_r.columns[1:].tolist() + ['scale']
                df_r = df_r.iloc[:, :-1]
                # df = df_r[['site_wdcgg_id', 'st_year', 'st_month', 'st_day', 'st_hour', 'st_minute', 'value']].copy()
                df = df_r[['st_year', 'st_month', 'st_day', 'st_hour', 'value']].copy()
                print(f"\tGet file date from WDCGG: {f_name}")
                print(df.head())
                print(df.tail())
                df.to_csv(self.obsout_dir + f_name + '.txt', sep='\t', index=False, header=False)

            except Exception as e:
                print(f"Error reading {filepath}: {e}")
                return None

        # =======================================
        # - obs
        df_obs = get_obs_sites()
        print(df_obs.head())
        print('\n\t Start sites iteration')

        # - loop sites
        for index, row in df_obs.iloc[:94].iterrows():
            idx, sno, f_name, lat, lon, alt = (row.iloc[0], row.iloc[1], row.iloc[2],
                                               row.iloc[3], row.iloc[4], row.iloc[5])
            files = find_files(f_name)
            if files != []:
                selected_file = find_file_with_largest_time_span(files)
                print('')
                read_wdcgg_file(selected_file)

    # --- get_obs_files_obspack
    def get_obs_files_obspack(self):
        # - get_obs_sites
        def get_obs_sites():
            df = pd.read_csv(self.sites_f, sep='\\s+', skiprows=2, header=None, usecols=range(6))
            df.columns = ['idx', 'sno', 'Sitename', 'Lat_o', 'Lon_o', 'Alt_o']
            return df

        # --- check_obspack
        def check_obspack(dir, f_name_):
            # - read file
            def get_obspack_ts():
                df_r = pd.read_csv(dir + filename, sep='\\s+', comment='#')
                df = df_r[['year', 'month', 'day', 'hour', 'value']].copy()
                df.rename(columns={'value': 'o_ch4'}, inplace=True)
                df.loc[:, 'o_ch4'] = df['o_ch4'] * 1e9
                return df

            found = False
            df_list = []

            # - check exact filename first
            for root, dirs, files in os.walk(dir):
                for filename in files:
                    if f_name_.lower() in filename.lower():
                        found = True
                        print(f"\tFound file: {filename}")
                        df = get_obspack_ts()
                        return found, df, None

            # - if no exact match, check first 7 chars
            if not found:
                for root, dirs, files in os.walk(dir):
                    for filename in files:
                        if f_name_[:11].lower() in filename.lower():
                            found = True
                            print(f"\tFound partial file: {filename}")
                            df = get_obspack_ts()
                            return found, df, filename

            if not found:
                print(f"\tNo file found for pattern: {f_name_}")
                return False, None

        # =======================================
        # - obs
        df_obs = get_obs_sites()
        print(df_obs.head())

        # - loop sites
        for index, row in df_obs.iloc[:94].iterrows():
            idx, sno, f_name, lat, lon, alt = (row.iloc[0], row.iloc[1], row.iloc[2],
                                               row.iloc[3], row.iloc[4], row.iloc[5])
            # if sno == 9999:
            if sno > -1:
                print(idx, sno, f_name)
                found = False

                # check obs
                found, df, f_name_n = check_obspack(self.obspack_dir, f_name)
                if found:
                    if f_name_n == None:
                        print(f"\tGet file date from Obspack: {f_name}")
                        df.to_csv(self.obsout_dir + f_name + '.txt', sep='\t', index=False, header=False)
                    else:
                        print(f"\t\tGet partial file date from Obspack: {f_name}")
                        df.to_csv(self.obsout_dir + f_name + '.txt', sep='\t', index=False, header=False)

                if not found:
                    print(f"\tNo file found for pattern: {f_name}")
                    exit()
