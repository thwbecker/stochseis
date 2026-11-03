#!/bin/bash
#
# Alpine Fault campaign with the L1 dip correction: ABC round 5 on each
# geometry, posterior ensembles of the best vectors with catalogs, and the
# two overview figures. Run from the repository root with the STL files in
# data/nz_alpine/.
#
# usage: examples/alpine/run_l1_campaign.sh [ncpu] [outdir] [n_best] [n_seeds] [duration_yr] ["geometries"]
#
ncpu=${1-$NR_CPUS}
outdir=${2-alpine_l1}
n_best=${3-20}
n_seeds=${4-5}
duration=${5-20000}

if [ -z "$ncpu" ]
then
    echo "$0: give the number of CPUs as first argument or set NR_CPUS" >&2
    exit 1
fi
mkdir -p "$outdir"
# geometries: planar and variable dip by default; add segmented as a third argument list if wanted
geoms=${6-"planar variable_dip"}
for g in $geoms
do
    python3 -m stochseis.cli sweep "examples/alpine/abc5_${g}_l1.json" "$outdir/abc5_$g" --nproc "$ncpu" --no-catalogs
    python3 examples/alpine/abc_accept.py "$outdir/abc5_$g/sweep_summary.json" 2 2.0 --seismicity > "$outdir/abc5_${g}_accept.txt"
    python3 examples/alpine/make_posterior_spec.py "$outdir/abc5_$g/sweep_summary.json" \
        "examples/alpine/abc5_${g}_l1.json" "$outdir/posterior5_$g.json" 2 "$n_best" "$n_seeds" "$duration"
    python3 -m stochseis.cli sweep "$outdir/posterior5_$g.json" "$outdir/posterior5_$g" --nproc "$ncpu"
done
pairs=""
for g in $geoms
do
    pairs="$pairs $g $outdir/posterior5_$g"
done
python3 examples/alpine/overview_figures.py "$outdir/fig" $pairs
python3 examples/alpine/site_analysis.py "$outdir/sites" $pairs --threshold 2
echo "done: acceptance tables in $outdir/abc5_*_accept.txt, figures in $outdir/fig_*.png and $outdir/sites_*"
