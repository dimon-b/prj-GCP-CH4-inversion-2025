# -*- coding: utf-8 -*-
"""Check GT3 framing and dates without NetCDF-Fortran compiler."""

import importlib.util
import struct
import tempfile
from pathlib import Path


def main():
    src=Path(__file__).resolve().parents[1]/'d0_inca_oh.py'
    spec=importlib.util.spec_from_file_location('d8',src)
    d8=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(d8)
    with tempfile.TemporaryDirectory() as td:
        p=Path(td)/'test.gt3'
        words=[' '*16 for _ in range(64)]
        words[2]='NDOH'.ljust(16)
        words[15]='molec/cm3'.ljust(16)
        words[26]='20000116 120000'.ljust(16)
        words[63]=str(128*64*67).rjust(16)
        data_bytes=128*64*67*4
        with p.open('wb') as f:
            header=''.join(words).encode('ascii')
            f.write(struct.pack('>I',len(header)))
            f.write(header)
            f.write(struct.pack('>I',len(header)))
            f.write(struct.pack('>I',data_bytes))
            f.write(b'\0'*data_bytes)
            f.write(struct.pack('>I',data_bytes))
        d8.verify_gt3(p,2000,1)
        with p.open('r+b') as f:
            f.seek(0)
            f.write(struct.pack('>I',666))
        try:
            d8.verify_gt3(p,2000,1)
        except ValueError:
            print('Malformed GTOOL3 frame rejected: PASS')
        else:
            raise AssertionError('Malformed GT3 was accepted')
    print('INCA OH GT3 structural test: PASS')


if __name__ == '__main__':
    main()
