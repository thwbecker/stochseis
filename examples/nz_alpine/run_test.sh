#!/bin/bash
#
# Runs still needed for the Alpine Fault example (examples/nz_alpine/README.md),
# from the repository root, with the meshes in data/nz_alpine/:
#
#   1. realization set D (barriers at the segment boundaries), 3 cases x
#      n_seeds realizations, and its figure pages (about 3 core-hours for
#      50 x 20,000 yr at 2 km patches);
#   2. what each control does, sets A to D (control_summary.py);
#   3. constraint-weighted predictions over all cases of sets A to D
#      (predictions.py; reads each catalog several times, about an hour);
#   4. if the ABC campaign directories exist, the parameter effects pooled
#      over rounds 5 to 7 (parameter_effects.py).
#
# Results are copied to examples/nz_alpine/results/ (tables and the figures
# used in the README). Sets A to C are expected in sets_dir already (run
# run_realization_sets.sh for them if not); everything is restartable.
#
# usage: examples/nz_alpine/run_test.sh [ncpu] [sets_dir] [n_seeds] [duration_yr] [campaign_dir_r5] [campaign_dir_r6] [campaign_dir_r7]
# defaults:                             $NR_CPUS alpine_sets 50 20000 alpine_full alpine_r6 alpine_r7
#
ncpu=${1-$NR_CPUS}
sets_dir=${2-alpine_sets}
n_seeds=${3-50}
duration=${4-20000}
r5=${5-alpine_full}
r6=${6-alpine_r6}
r7=${7-alpine_r7}

here=examples/nz_alpine
res=$here/results
if [ -z "$ncpu" ]
then
    echo "$0: give the number of CPUs as first argument or set NR_CPUS" >&2
    exit 1
fi
bash $here/run_realization_sets.sh "$ncpu" "$sets_dir" "$n_seeds" "$duration" "D_boundaries" no || exit 1

python3 $here/control_summary.py "$sets_dir/controls" "$sets_dir" A_clock B_size C_barriers D_boundaries --site north=1

pairs=""
for s in A_clock B_size C_barriers D_boundaries
do
    while read label vec seed0
    do
        pairs="$pairs ${s}_${label} $sets_dir/$s:$vec"
    done < "$sets_dir/specs/realizations_${s}_cases.txt"
done
python3 $here/predictions.py "$sets_dir/predictions" $pairs --site north=1

summaries=""
for d in "$r5" "$r6" "$r7"
do
    summaries="$summaries $(ls $d/abc_*/sweep_summary.json 2>/dev/null)"
done
if [ -n "$(echo $summaries)" ]
then
    python3 $here/parameter_effects.py "$sets_dir/effects" $summaries --site north=1
    cp "$sets_dir/effects_effects.png" "$sets_dir/effects_dependence.png" "$sets_dir/effects_effects.md" "$res/sets/"
fi

mkdir -p "$res/sets/figures"
cp "$sets_dir/controls.png" "$sets_dir/controls.md" "$sets_dir/predictions.png" "$sets_dir/predictions.md" "$res/sets/"
for s in A_clock B_size C_barriers D_boundaries
do
    cp "$sets_dir/${s}_ensemble_p2.md" "$sets_dir/${s}_ensemble_p4.md" "$sets_dir/${s}_ensemble_p3_table.md" "$res/sets/" 2>/dev/null
    cp "$sets_dir/${s}_ensemble_p1_seismicity.png" "$sets_dir/${s}_ensemble_p2.png" "$sets_dir/${s}_ensemble_p4.png" \
        "$res/sets/figures/" 2>/dev/null
done
echo "done: results copied to $res/sets/"
