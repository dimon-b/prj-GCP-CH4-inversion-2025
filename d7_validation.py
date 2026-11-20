# -*- coding: utf-8 -*-
"""D7: independently validate forward CH4: surface, HIPPO, and growth rate."""
import importlib


VALIDATION = (
    ('v1_check_3sites', 'Check3Sites'),
    ('v2_check_all_surf', 'CheckAllSurf'),
    ('v3_check_hippo', 'CheckHippo'),
    ('v4_check_GR', 'CheckGR'),
)


def run_stage(case='INCA', dry_run=False):
    print('D7 VALIDATION:', str(case).upper())
    print('Sequence:', ', '.join(s[0] for s in VALIDATION))
    if dry_run:
        print('DRY-RUN: no validation files written')
        return
    for name, clsname in VALIDATION:
        print('\n===', name, '===')
        module = importlib.import_module(name)
        cls = getattr(module, clsname)
        cls()
