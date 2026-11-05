"""Write one ABC spec per campaign configuration.

    python3 examples/nz_alpine/make_campaign_specs.py configs.json base_spec.json outdir \\
        [n_samples=100] [n_seeds=3] [duration_yr=10000] [names...]

Every spec takes priors, design seed and fixed parameters from base_spec
(round 5: hybrid sizes, deficit-capped slip, L1 dip correction) and
geometry file, patch size and segment matching distance from the
configuration (and slip maps stored from M5.5 for 1 km patches), so all configurations draw the same parameter vectors and
can be compared design by design. Writes <outdir>/abc_<name>.json and
checks that each geometry file exists.
"""
import json
import os
import sys


def main():
    cfg_file, base_file, outdir = sys.argv[1:4]
    n_samples = int(sys.argv[4]) if len(sys.argv) > 4 else 100
    n_seeds = int(sys.argv[5]) if len(sys.argv) > 5 else 3
    duration = float(sys.argv[6]) if len(sys.argv) > 6 else 10000.0
    names = sys.argv[7:]
    with open(cfg_file) as fh:
        configs = json.load(fh)["configs"]
    with open(base_file) as fh:
        base = json.load(fh)
    os.makedirs(outdir, exist_ok=True)
    missing = []
    for c in configs:
        if names and c["name"] not in names:
            continue
        if not os.path.exists(c["geometry_file"]):
            missing.append(c["geometry_file"])
        spec = json.loads(json.dumps(base))
        spec["_comment"] = (f"ABC for configuration {c['name']} ({c['geometry']}, {c['resolution']}: "
                            f"{os.path.basename(c['geometry_file'])} on {c['patch_km']:g} km patches)")
        spec["n_samples"] = n_samples
        spec["seeds"] = list(range(1, n_seeds + 1))
        spec["base"].update(geometry_file=c["geometry_file"], patch_km=c["patch_km"],
                            segment_max_dist_km=c["segment_max_dist_km"], duration_yr=duration)
        # with 1 km patches an event of M5.5 to 6 can put a metre of slip on a site patch: store its slip
        if c["patch_km"] <= 1.0:
            spec["base"]["store_slip_min_mag"] = min(spec["base"].get("store_slip_min_mag", 6.0), 5.5)
        out = os.path.join(outdir, f"abc_{c['name']}.json")
        with open(out, "w") as fh:
            json.dump(spec, fh, indent=1)
        print(out)
    if missing:
        sys.exit("missing geometry files: " + ", ".join(missing))


if __name__ == "__main__":
    main()
