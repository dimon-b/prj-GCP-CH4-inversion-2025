# -*- coding: utf-8 -*-
"""
Created:    21/11/08 17:48
Project:    GCP-21
@author:    Dmitry Belikov
"""

import os
import subprocess

import _set_case


class Write_grd(_set_case.SetCase):

    def __init__(self):
        super().__init__()

        # --- convert predicted.sources to 2D grd files by FORTRAN
        self.run_fort_grd()

    # --- run_fort_grd
    def run_fort_grd(self):

        src_dir = "../c_gcpv3_f"
        f_file = os.path.join(src_dir, "s2_ch4_grd.f90")
        exe_file = os.path.join(src_dir, "a.out")

        # --- compile
        if os.path.exists(exe_file):
            os.remove(exe_file)

        compile_cmd = ["ifort", "-O3", "s2_ch4_grd.f90"]
        process = subprocess.Popen(compile_cmd,
                                   cwd=src_dir,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   )
        output, error = process.communicate()

        if error:
            print("\trun_fort_grd compilation stderr:")
            print(error.decode())
            raise RuntimeError("Compilation failed")

        if not os.path.exists(exe_file):
            raise FileNotFoundError("s2_ch4_grd.exe not found")

        # --- run
        run_cmd = ["./a.out"]
        process = subprocess.Popen(
            run_cmd,
            cwd=src_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        output, error = process.communicate()

        if error:
            print("\trun_fort_grd execution stderr:")
            print(error.decode())
            raise RuntimeError("Execution failed")

        print("\trun_fort_grd output ::")
        for line in output.decode().splitlines():
            if line.strip():
                print("\t>>", line)