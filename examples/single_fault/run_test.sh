#!/bin/bash
# Copyright (C) 2026 Thorsten Becker. Written by Claude (Anthropic) with Thorsten Becker.
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Single-fault parameter overview (examples/single_fault/README.md): the
# default 200 km strike-slip fault and two alternative values each of the
# branching ratio, the nucleation clock and the size rule, n_seeds
# realizations of duration_yr each (catalogs kept), then the figures and
# table referenced by the README. Run from the repository root.
# 7 cases x 10 x 10,000 yr is about 2 core-hours; restartable.
#
# usage: examples/single_fault/run_test.sh [ncpu] [outdir] [n_seeds] [duration_yr] [figure_dir]
# defaults:                                $NR_CPUS single_fault_out 10 10000 examples/single_fault/figures
#
ncpu=${1-$NR_CPUS}
outdir=${2-single_fault_out}
n_seeds=${3-10}
duration=${4-10000}
figdir=${5-examples/single_fault/figures}

here=examples/single_fault
if [ -z "$ncpu" ]
then
    echo "$0: give the number of CPUs as first argument or set NR_CPUS" >&2
    exit 1
fi
mkdir -p "$outdir"
python3 - "$here/specs/parameter_overview.json" "$outdir/spec.json" "$n_seeds" "$duration" <<'PY'
import json, sys
spec = json.load(open(sys.argv[1]))
spec["seeds"] = list(range(1, int(sys.argv[3]) + 1))
spec["base"]["duration_yr"] = float(sys.argv[4])
json.dump(spec, open(sys.argv[2], "w"), indent=1)
PY
python3 -m stochseis.cli sweep "$outdir/spec.json" "$outdir/runs" --nproc "$ncpu" || exit 1
python3 $here/make_figures.py "$outdir/spec.json" "$outdir/runs" "$figdir"
echo "done: figures and table in $figdir"
