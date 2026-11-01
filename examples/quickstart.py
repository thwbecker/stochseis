"""Minimal end-to-end example: run, analyse, plot."""
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from stochseis import Params, analysis, plots, run

p = Params(duration_yr=1000, seed=1, fault_length_km=200)
cat = run(p)
cat.save("quickstart.npz")
m = analysis.summary_metrics(cat)
for k in ("b_fault_below_sat", "bath_mean", "omori_p", "cv_m6.5_declustered", "moment_closure"):
    print(k, m[k])
plots.summary_figure(cat, "quickstart_summary.png")
plots.deficit_figure(cat, "quickstart_deficit.png")
