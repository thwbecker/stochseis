"""Write one sweep spec per realization set.

    python3 examples/alpine/make_realization_specs.py sets.json outdir [n_seeds=50] [duration_yr=20000] [set names...]

For every set in sets.json the spec has design "list": one vector per value
of the set's parameter (the reference vector with that parameter replaced,
plus the set's fixed parameters), each run for n_seeds seeds starting at the
set's seed0. Base parameters (slip mode, merging, sites, ...) come from the
base spec, geometry and patch size from the configuration. Writes
<outdir>/realizations_<set>.json and <outdir>/realizations_<set>_cases.txt
(one "label vNNNN" line per case, for the figure scripts).
"""
import json
import os
import sys


def main():
    sets_file, outdir = sys.argv[1:3]
    n_seeds = int(sys.argv[3]) if len(sys.argv) > 3 else 50
    duration = float(sys.argv[4]) if len(sys.argv) > 4 else 20000.0
    names = sys.argv[5:]
    with open(sets_file) as fh:
        cfg = json.load(fh)
    with open(cfg["base_spec"]) as fh:
        base = json.load(fh)["base"]
    base.update(cfg["configuration"])
    base["duration_yr"] = duration
    if not os.path.exists(base["geometry_file"]):
        sys.exit(f"missing geometry file {base['geometry_file']}")
    os.makedirs(outdir, exist_ok=True)
    for name, st in cfg["sets"].items():
        if names and name not in names:
            continue
        vectors = []
        for v in st["values"]:
            vec = dict(cfg["reference"])
            vec.update(st.get("fixed", {}))
            vec[st["parameter"]] = v
            vectors.append(vec)
        seed0 = int(st.get("seed0", 1))
        spec = {"_comment": f"Realization set {name}: {st['parameter']} = {st['values']} ({', '.join(st['labels'])}); "
                            f"{n_seeds} realizations of {duration:.0f} yr each.",
                "design": "list", "base": base, "vectors": vectors, "case_names": st["labels"],
                "seeds": list(range(seed0, seed0 + n_seeds)), "mmain": 6.0}
        out = os.path.join(outdir, f"realizations_{name}.json")
        with open(out, "w") as fh:
            json.dump(spec, fh, indent=1)
        with open(os.path.join(outdir, f"realizations_{name}_cases.txt"), "w") as fh:
            for i, lab in enumerate(st["labels"]):
                fh.write(f"{lab} v{i:04d} {seed0}\n")
        print(out)


if __name__ == "__main__":
    main()
