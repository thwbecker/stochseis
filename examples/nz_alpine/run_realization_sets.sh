#!/bin/bash
#
# Realization sets for the Alpine Fault: families of cases that vary one
# irregularity control at a time (examples/nz_alpine/realization_sets.json),
# each case run as n_seeds independent realizations with catalogs, and the
# three figure pages for one realization and for the ensemble:
#
#   page 1  <set>_<mode>_p1_seismicity.png   frequency-magnitude, Omori, Bath, moment release,
#                                            whole-fault recurrence (overview_figures.py; its
#                                            recurrence figure is written as <set>_<mode>_p1_recurrence.png)
#   page 2  <set>_<mode>_p2.png / .md        site recurrence in record-length windows, SED versus RI,
#                                            rupture patterns, magnitudes, rupture lengths (nz_page.py)
#   page 3  <set>_<mode>_p3_statistics.png   intercomparison panels (site_analysis.py; also
#                                            _spacetime.png, _extended.png, _table.md)
#   page 4  <set>_ensemble_p4.png / .md      how often synthetic records of the observed length match
#                                            the observed RI and CV, per site and jointly (record_match.py)
#   predictions.png / .md                    all cases pooled, realizations weighted by how often their
#                                            records match the observations (predictions.py)
#
# with mode "single" (first realization of each case) and "ensemble" (all
# realizations; bands are 5-95 % over realizations). Run from the repository
# root with the variable-dip mesh in data/nz_alpine/. Restartable: finished
# runs are skipped.
#
# usage: examples/nz_alpine/run_realization_sets.sh [ncpu] [outdir] [n_seeds] [duration_yr] ["sets"] [predictions]
# defaults:                                      $NR_CPUS alpine_sets 50 20000 "A_clock B_size C_barriers D_boundaries" yes
#
ncpu=${1-$NR_CPUS}
outdir=${2-alpine_sets}
n_seeds=${3-50}
duration=${4-20000}
sets=${5-"A_clock B_size C_barriers D_boundaries"}
predictions=${6-yes}

here=examples/nz_alpine
if [ -z "$ncpu" ]
then
    echo "$0: give the number of CPUs as first argument or set NR_CPUS" >&2
    exit 1
fi
mkdir -p "$outdir"
python3 $here/make_realization_specs.py $here/realization_sets.json "$outdir/specs" "$n_seeds" "$duration" $sets || exit 1

for s in $sets
do
    echo "=== set $s: runs"
    python3 -m stochseis.cli sweep "$outdir/specs/realizations_$s.json" "$outdir/$s" --nproc "$ncpu" || exit 1
    for mode in single ensemble
    do
        pairs=""
        while read label vec seed0
        do
            if [ "$mode" = single ]
            then
                pairs="$pairs $label $outdir/$s:${vec}_s${seed0}"
            else
                pairs="$pairs $label $outdir/$s:$vec"
            fi
        done < "$outdir/specs/realizations_${s}_cases.txt"
        echo "=== set $s: $mode pages"
        python3 $here/overview_figures.py "$outdir/${s}_${mode}_p1" $pairs
        python3 $here/nz_page.py "$outdir/${s}_${mode}_p2" $pairs --site north=1
        python3 $here/site_analysis.py "$outdir/${s}_${mode}_p3" $pairs --threshold 2
        if [ "$mode" = ensemble ]
        then
            python3 $here/record_match.py "$outdir/${s}_${mode}_p4" $pairs --site north=1
        fi
    done
done
# constraint-weighted predictions over all cases of all sets
if [ "$predictions" != yes ]
then
    echo "done: pages in $outdir/<set>_<single|ensemble>_p{1,2,3}*, $outdir/<set>_ensemble_p4*"
    exit 0
fi
pairs=""
for s in $sets
do
    while read label vec seed0
    do
        pairs="$pairs ${s}_${label} $outdir/$s:$vec"
    done < "$outdir/specs/realizations_${s}_cases.txt"
done
echo "=== predictions over all cases"
python3 $here/predictions.py "$outdir/predictions" $pairs --site north=1
echo "done: pages in $outdir/<set>_<single|ensemble>_p{1,2,3}*, $outdir/<set>_ensemble_p4*, $outdir/predictions.*"
