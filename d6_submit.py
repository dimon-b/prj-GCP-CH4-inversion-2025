# -*- coding: utf-8 -*-
"""D6: model experiment/submission preparation; delegates to live implementation."""
import importlib


def run_stage(case='INCA', dry_run=False):
    print('D6 SUBMITION:', str(case).upper())
    print('Procedure: existing t0_submit.Write2Submit')
    if dry_run:
        print('DRY-RUN: no submission scripts written or submitted')
        return
    m = importlib.import_module('t0_submit')
    m.Write2Submit()
