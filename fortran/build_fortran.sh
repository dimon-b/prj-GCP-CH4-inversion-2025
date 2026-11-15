#!/usr/bin/env bash
# ============================================================
# === Compile GCP/GMB 2025 Fortran programs once for CYC and INCA
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p bin

if [[ -n "${FC:-}" ]]; then
    compiler="$FC"
elif command -v ifort >/dev/null 2>&1; then
    compiler=ifort
else
    compiler=gfortran
fi

if ! command -v "$compiler" >/dev/null 2>&1; then
    echo "Compiler not found: $compiler" >&2
    exit 1
fi

echo "Compiler: $compiler"
flags=(-O2)
if [[ "$compiler" == *gfortran* ]]; then
    flags+=(-std=legacy -ffree-line-length-none)
fi

"$compiler" "${flags[@]}" f2_actm_ch4_at_obs_sites.f90 -o bin/f2_actm_ch4_at_obs_sites
"$compiler" "${flags[@]}" f3_combch4_obs_actm.f90 -o bin/f3_combch4_obs_actm
"$compiler" "${flags[@]}" f4_ch4bu_mon_reg_t42.f90 -o bin/f4_ch4bu_mon_reg_t42
"$compiler" "${flags[@]}" f6_nakazawa_filter.f90 lsfbf4_main.f lsfbf4_sub.f -o bin/f6_nakazawa_filter

if command -v nf-config >/dev/null 2>&1; then
    read -ra nf_fflags <<< "$(nf-config --fflags)"
    read -ra nf_libs <<< "$(nf-config --flibs)"
    "$compiler" "${flags[@]}" "${nf_fflags[@]}" f0_fch4_gcp25_nc2gt3_t42.f90 -o bin/f0_fch4_gcp25_nc2gt3_t42 "${nf_libs[@]}"
    "$compiler" "${flags[@]}" "${nf_fflags[@]}" f5_write_prior_ch4.f90 -o bin/f5_write_prior_ch4 "${nf_libs[@]}"
    "$compiler" "${flags[@]}" "${nf_fflags[@]}" f7_write_tdi_obsrvCH4dif.f90 -o bin/f7_write_tdi_obsrvCH4dif "${nf_libs[@]}"
else
    echo 'WARNING: nf-config is not installed/activated; F0, F5, F7 were NOT compiled.' >&2
fi

echo "Binaries: $(pwd)/bin/"
