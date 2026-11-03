"""Write a sweep spec that reruns the best ABC vectors with catalogs kept.

    python3 examples/alpine/make_posterior_spec.py summary.(json|csv) base_spec.json out_spec.json \\
        [threshold_m=2] [n_best=20] [n_seeds=5] [duration_yr=20000] [--seismicity]
        [--site north=1 ...] [--bath-sigma 0.15] [--bath-m7]

Vectors are ranked with abc_scoring.score (with the seismicity terms if
--seismicity is given). The output spec has design "list": every selected
vector runs for every seed; run it without --no-catalogs so that event
catalogs are available for the figures. The base parameters, including
geometry and patch size, come from base_spec.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from abc_scoring import load_summary, parse_cli, score  # noqa: E402


def main():
    argv, kw, seis = parse_cli(sys.argv)
    summary, base_spec, out_spec = argv[1:4]
    th = float(argv[4]) if len(argv) > 4 else 2.0
    n_best = int(argv[5]) if len(argv) > 5 else 20
    n_seeds = int(argv[6]) if len(argv) > 6 else 5
    duration = float(argv[7]) if len(argv) > 7 else 20000.0
    with open(base_spec) as fh:
        spec = json.load(fh)
    varied = list(spec.get("ranges", {})) + list(spec.get("log_ranges", {})) + list(spec.get("choices", {}))
    scored = score(load_summary(summary), th, seis, **kw)[:n_best]
    vectors = [{p: x["rows"][0]["p_" + p] for p in varied if "p_" + p in x["rows"][0]} for x in scored]
    out = {"_comment": f"Posterior ensemble: best {len(vectors)} vectors of {summary} (threshold {th} m"
                       f"{', with seismicity' if seis else ''}), {n_seeds} seeds, catalogs kept.",
           "design": "list", "base": dict(spec["base"]), "vectors": vectors,
           "misfits": [x["misfit"] for x in scored], "seeds": list(range(101, 101 + n_seeds)),
           "mmain": spec.get("mmain", 6.0)}
    out["base"]["duration_yr"] = duration
    with open(out_spec, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"wrote {out_spec}: {len(vectors)} vectors x {n_seeds} seeds, "
          f"misfits {out['misfits'][0]:.2f} to {out['misfits'][-1]:.2f}")


if __name__ == "__main__":
    main()
