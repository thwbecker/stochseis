#!/bin/bash
# Copyright (C) 2026 Thorsten Becker. Written by Claude (Anthropic) with Thorsten Becker.
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Run the stochseis test suite and the experiment sweeps that need more
# than one core. Run from the root of the stochseis repository.
#
# usage: examples/single_fault/experiments/run_stochseis_tests.sh [ncpu] [outdir] [all]
#
#   ncpu    number of processes (default: $NR_CPUS)
#   outdir  output directory (default: experiments_out)
#   all     also rerun E5 to E10 (about 6 core-hours)
#
ncpu=${1-$NR_CPUS}
outdir=${2-experiments_out}
mode=${3-short}

if [ -z "$ncpu" ]
then
    echo "$0: give the number of CPUs as first argument or set NR_CPUS" >&2
    exit 1
fi
if [ ! -d stochseis ] || [ ! -d examples/single_fault/experiments ]
then
    echo "$0: run from the root of the stochseis repository" >&2
    exit 1
fi

mkdir -p "$outdir"
log=$outdir/run.log
echo "$(date): ncpu=$ncpu outdir=$outdir mode=$mode" | tee "$log"

# 1. unit and regression tests
python3 -m pytest tests -q 2>&1 | tee -a "$log"
if [ "${PIPESTATUS[0]}" -ne 0 ]
then
    echo "$0: tests failed, stopping" | tee -a "$log"
    exit 1
fi

# 2. half-space check (needs cutde; skipped if it is not installed)
if python3 -c "import cutde" 2>/dev/null
then
    python3 examples/checks/halfspace_check.py 2>&1 | tee "$outdir/halfspace_check.txt"
else
    echo "cutde not installed, skipping half-space check" | tee -a "$log"
fi

# 3. experiment sweeps
specs="examples/single_fault/experiments/e22_defaults_error_bars.json examples/single_fault/experiments/e21_stressdrop_deepslip.json"
if [ "$mode" = "all" ]
then
    specs="$specs $(ls examples/single_fault/experiments/e[5-9]_*.json examples/single_fault/experiments/e10_*.json)"
fi

for spec in $specs
do
    name=$(basename "$spec" .json)
    echo "$(date): $name" | tee -a "$log"
    python3 -m stochseis.cli sweep "$spec" "$outdir/$name" --nproc "$ncpu" --no-catalogs 2>&1 | tee -a "$log"
done

# 4. seed-averaged tables
python3 examples/single_fault/experiments/aggregate_sweep.py "$outdir/e22_defaults_error_bars/sweep_summary.json" \
    p_duration_yr > "$outdir/e22_table.txt"
python3 examples/single_fault/experiments/aggregate_sweep.py "$outdir/e21_stressdrop_deepslip/sweep_summary.json" \
    p_stress_drop_mpa p_deep_slip_km > "$outdir/e21_table.txt"
echo "$(date): done; tables in $outdir/e2[12]_table.txt, summaries in $outdir/*/sweep_summary.json" | tee -a "$log"
