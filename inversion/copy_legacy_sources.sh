#!/usr/bin/env bash
# ============================================================
# === Copy only the unchanged legacy inversion sources
#     into this directory, without overwriting the patched ones.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
from="${1:-$here/../../codeint54}"
from="$(cd "$from" && pwd)"
cd "$here"

if [[ "$from" == "$here" ]]; then
    echo 'ERROR: source and destination are the same directory.' >&2
    exit 1
fi

missing=0
copied=0
for base in \
    bayesinv bldfrech bldgreens datcor dumpgreen findname \
    getconsts getdat getlgreen getnum getsrc getsub \
    output_v3 srccor lnblnk fnblnk defgrps getgrps \
    dump recipe getblklst split zero getrecipe \
    getmap mapinterp influence diffdat intdat svdcalc \
    inrange getmats tunevar packmat read_pseudodata read_sourcep \
    output_tc3 read_monthlyobs read_sourcetdi tc3greens
do
    if [[ -f "$base.f90" || -f "$base.f" ]]; then
        continue
    fi
    if [[ -f "$from/$base.f90" ]]; then
        cp -p "$from/$base.f90" .
    elif [[ -f "$from/$base.f" ]]; then
        cp -p "$from/$base.f" .
    else
        echo "MISSING legacy source: $base.f90 / $base.f" >&2
        missing=$((missing + 1))
        continue
    fi
    copied=$((copied + 1))
done

for header in tdi.h common.h pointer.h; do
    if [[ ! -f "$header" ]]; then
        if [[ -f "$from/$header" ]]; then
            cp -p "$from/$header" .
            copied=$((copied + 1))
        else
            echo "MISSING header: $header" >&2
            missing=$((missing + 1))
        fi
    fi
done

# --- Master/recipe files are runtime data, not compiler sources.
if [[ -d "$from/files4run" ]]; then
    mkdir -p files4run
    cp -an "$from/files4run/." files4run/
    echo 'Copied missing files4run/ assets (existing files preserved).'
else
    echo 'WARNING: files4run/ not found under the original directory.' >&2
    echo 'Copy it manually before running the Bayesian inversion.' >&2
fi

printf 'Legacy sources/headers newly copied: %s; missing: %s\n' "$copied" "$missing"
if (( missing )); then exit 1; fi
