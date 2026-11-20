# -*- coding: utf-8 -*-
"""D2: ObsPack/WDCGG observation preparation."""
import importlib


def run_stage(case='INCA', dry_run=False):
    print('D2 OBS: prepare station observation files')
    if dry_run:
        print('DRY-RUN: no observation files created')
        return
    importlib.import_module('d3_obs_files').ObsFilesTxt()
