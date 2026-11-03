#!/bin/bash
#
# Full Alpine Fault campaign: for every configuration in
# examples/alpine/campaign_configs.json (three geometries on 0.5 km meshes
# with 1 km patches, two on 1 km meshes with 2 km patches), an ABC round
# (hybrid sizes, deficit-capped slip, L1 dip correction; same parameter
# vectors in every configuration), acceptance tables, a posterior ensemble
# of the best vectors with catalogs, and then summary tables and figures
# for the fine set, the coarse set, and fine against coarse for each
# geometry available at both resolutions.
#
# Run from the repository root with the mesh files in data/nz_alpine/.
# Interrupted runs can be restarted with the same command: finished sweep
# members are skipped.
#
# usage: examples/alpine/run_full_campaign.sh [ncpu] [outdir] [n_samples] [abc_seeds] [abc_duration_yr]
#                                             [n_best] [post_seeds] [post_duration_yr] ["configurations"]
# defaults:                                   $NR_CPUS alpine_full 100 3 10000 20 5 20000
#                                             "planar_fine variable_dip_fine segmented_fine planar_coarse variable_dip_coarse"
#
ncpu=${1-$NR_CPUS}
outdir=${2-alpine_full}
n_samples=${3-100}
abc_seeds=${4-3}
abc_dur=${5-10000}
n_best=${6-20}
post_seeds=${7-5}
post_dur=${8-20000}
configs=${9-"planar_fine variable_dip_fine segmented_fine planar_coarse variable_dip_coarse"}

here=examples/alpine
if [ -z "$ncpu" ]
then
    echo "$0: give the number of CPUs as first argument or set NR_CPUS" >&2
    exit 1
fi
mkdir -p "$outdir"
python3 $here/make_campaign_specs.py $here/campaign_configs.json $here/abc5_planar_l1.json "$outdir/specs" \
    "$n_samples" "$abc_seeds" "$abc_dur" $configs || exit 1

for c in $configs
do
    echo "=== $c: ABC"
    python3 -m stochseis.cli sweep "$outdir/specs/abc_$c.json" "$outdir/abc_$c" --nproc "$ncpu" --no-catalogs || exit 1
    python3 $here/abc_accept.py "$outdir/abc_$c/sweep_summary.json" 2 2.0 > "$outdir/abc_${c}_accept_paleo.txt"
    python3 $here/abc_accept.py "$outdir/abc_$c/sweep_summary.json" 2 2.0 --seismicity > "$outdir/abc_${c}_accept_seismicity.txt"
    echo "=== $c: posterior"
    python3 $here/make_posterior_spec.py "$outdir/abc_$c/sweep_summary.json" "$outdir/specs/abc_$c.json" \
        "$outdir/specs/posterior_$c.json" 2 "$n_best" "$post_seeds" "$post_dur" --seismicity || exit 1
    python3 -m stochseis.cli sweep "$outdir/specs/posterior_$c.json" "$outdir/posterior_$c" --nproc "$ncpu" || exit 1
done

python3 $here/campaign_summary.py "$outdir" $configs > "$outdir/campaign_summary.md"

# figure sets: label/directory pairs for the configurations that were run
pairs_for()
{
    local out=""
    for c in "$@"
    do
        if [ -d "$outdir/posterior_$c" ]
        then
            out="$out $c $outdir/posterior_$c"
        fi
    done
    echo $out
}
figure_set()
{
    local name=$1
    shift
    local pairs
    pairs=$(pairs_for "$@")
    if [ $(echo $pairs | wc -w) -ge 2 ]
    then
        echo "=== figures: $name"
        python3 $here/overview_figures.py "$outdir/fig_$name" $pairs
        python3 $here/site_analysis.py "$outdir/sites_$name" $pairs --threshold 2
    fi
}
fine=""
coarse=""
for c in $configs
do
    case $c in
        *_fine) fine="$fine $c" ;;
        *_coarse) coarse="$coarse $c" ;;
    esac
done
figure_set fine $fine
figure_set coarse $coarse
for c in $coarse
do
    g=${c%_coarse}
    figure_set "resolution_$g" "${g}_fine" "$c"
done
echo "done: $outdir/campaign_summary.md, $outdir/abc_*_accept_*.txt, $outdir/fig_*.png, $outdir/sites_*"
