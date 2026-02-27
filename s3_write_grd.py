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

        # --- Paths
        source_dir = os.path.abspath("../c_gcpv3_f")
        fortran_file = os.path.join(source_dir, "s2_ch4_grd.f90")
        exe_name = "a.out"
        exe_path = os.path.join(source_dir, exe_name)

        # --- Patch FORTRAN source (line 6 -> index 5)
        try:
            with open(fortran_file, "r", encoding="utf-8") as f:
                lines = f.readlines()

            lines[5] = f"character(len=*),parameter    :: work_dir = '{self.inv_wrk_dir}'\n"

            with open(fortran_file, "w", encoding="utf-8") as f:
                f.writelines(lines)

            print("Patched s2_ch4_grd.f90:")
            print(f"    work_dir = {self.inv_wrk_dir}")
            print('')


        except Exception as e:
            raise RuntimeError(f"ERROR patching s2_ch4_grd.f90: {e}")

        # --- Remove existing executable
        if os.path.isfile(exe_path):
            try:
                os.remove(exe_path)
                print(f"Removed existing {exe_path}")
            except Exception as e:
                print(f"Warning: could not remove {exe_path}: {e}")

        # --- Compile (IN source_dir)
        compile_cmd = ["ifort", "-O3", os.path.basename(fortran_file), "-o", exe_name]

        print(f"Compiling in {source_dir}")
        print("Command:", " ".join(compile_cmd))

        try:
            proc = subprocess.run(
                compile_cmd,
                cwd=source_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            if proc.returncode != 0:
                print(proc.stdout)
                print(proc.stderr)
                raise RuntimeError("Compilation failed")

            if proc.stderr.strip():
                print("Compiler warnings:")
                print(proc.stderr)

            print("Compilation successful")

        except FileNotFoundError:
            raise RuntimeError("ERROR: ifort not found in PATH")
        except Exception as e:
            raise RuntimeError(f"ERROR during compilation: {e}")

        # --- Verify executable
        if not os.path.isfile(exe_path):
            raise RuntimeError(f"ERROR: executable not created: {exe_path}")

        # --- Execute (IN source_dir)
        execute_cmd = [f"./{exe_name}"]

        print(f"Executing in {source_dir}: {' '.join(execute_cmd)}")

        try:
            proc = subprocess.run(
                execute_cmd,
                cwd=source_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            print(f"Return code: {proc.returncode}")

            if proc.stderr.strip():
                print("Runtime stderr:")
                print(proc.stderr)

            if proc.returncode != 0:
                raise RuntimeError(f"{exe_name} failed")

            print(f"{exe_name} output:")
            if proc.stdout.strip():
                for line in proc.stdout.splitlines():
                    print(f"\t>> {line}")
            else:
                print("\t>> No output")

        except Exception as e:
            raise RuntimeError(f"ERROR executing {exe_name}: {e}")

        print(f"{exe_name} completed successfully")