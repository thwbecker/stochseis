"""Plot summary metrics of a sweep against arrest_scale, grouped by stress drop.

    python examples/plot_sweep.py sweep_out/sweep_summary.json sweep.png
"""
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import json  # noqa: E402

from stochseis.plots import sweep_figure

table = json.load(open(sys.argv[1]))
keys = ["b_fault_below_sat", "b_fault_above_sat", "bath_mean", "cv_m6.5_declustered",
        "omori_p", "frac_stop_deficit_m6", "max_mag", "closure_after_burnin"]
sweep_figure(table, "p_arrest_scale", keys, path=sys.argv[2], groupkey="p_stress_drop_mpa")
print("wrote", sys.argv[2])
