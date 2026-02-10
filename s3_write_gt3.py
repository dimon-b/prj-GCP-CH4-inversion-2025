# -*- coding: utf-8 -*-
"""
Finalized:  26/xx/xx
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""

import os
import subprocess

import _set_case


class Write_gt3(_set_case.SetCase):

    def __init__(self):
        super().__init__()

        # --- write gt3 flux by FORTRAN
        self.run_fort_gt3()

    # ---
    def run_fort_gt3(self):
        src_dir = "../c_gcpv3_f"
        f_file = os.path.join(src_dir, "s3_ch4_grd2gt3.f90")
        exe_file = os.path.join(src_dir, "a.out")

        # --- compile
        compile_cmd = ["ifort", "-O3", f_file]
        process = subprocess.Popen(
            compile_cmd,
            cwd=src_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        output, error = process.communicate()

        if error:
            print("\trun_fort_gt3 compilation stderr:")
            print(error.decode())
            raise RuntimeError("Compilation failed")

        if not os.path.exists(exe_file):
            raise FileNotFoundError("a.out not found")

        # --- run
        run_cmd = ["./a.out"]
        process = subprocess.Popen(run_cmd,
                                   cwd=src_dir,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   )
        output, error = process.communicate()

        if error:
            print("\trun_fort_gt3 execution stderr:")
            print(error.decode())
            raise RuntimeError("Execution failed")

        print("\trun_fort_gt3 output ::")
        for line in output.decode().splitlines():
            if line.strip():
                print("\t>>", line)
