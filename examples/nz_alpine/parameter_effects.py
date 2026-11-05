"""How each model parameter affects performance, pooled over ABC sweeps.

    python3 examples/nz_alpine/parameter_effects.py out_prefix summary1 [summary2 ...] [--site north=1]

Each summary is a sweep_summary.json or .csv of an ABC round (for example
rounds 5, 6 and 7, which together span wide ranges of the area threshold,
nucleation exponent, branching ratio and b). Every design vector (mean
over its seeds) is one sample. Parameters that vary in the pooled sample
are used as predictors (log-transformed when they span more than a factor
of 5); performance metrics:

  CV, RI, SED at the sites (2 m; per-site thresholds with --site)
  b of grown ruptures below width saturation, excess of M >= 7 events over
  the b extrapolation from M >= 6, Bath gap (M >= 6; M >= 7 if summarized),
  Omori p, maximum magnitude, CV of whole-fault M >= 7 recurrence
  (declustered), central-south co-rupture; and whether the vector is
  accepted (paleo + seismicity, Bath not constrained), shown as acceptance
  rate in the dependence figure only

Outputs
  <out_prefix>_effects.png      standardized regression coefficients (metric sd per parameter sd) with
                                R^2 per metric: which parameter moves which metric
  <out_prefix>_dependence.png   binned medians and quartiles of key metrics against the most
                                influential parameters, and the acceptance rate per bin
  <out_prefix>_effects.md       the coefficient table
"""
import os
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
warnings.filterwarnings("ignore")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from abc_scoring import accepted, load_summary, parse_cli, score  # noqa: E402
from stochseis.params import Params  # noqa: E402

PARAM_FIELDS = set(Params.__dataclass_fields__)

KEY_METRICS = ["central CV", "south CV", "Bath M>=6", "excess M>=7", "central RI", "accepted"]


def _f(r, k):
    v = r.get(k)
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def design_table(rows, site_th):
    """One record per design vector: parameter values and mean metrics over its seeds."""
    groups = {}
    for r in rows:
        if r.get("error"):
            continue
        groups.setdefault((r["_src"], r.get("design_index", r.get("file"))), []).append(r)
    out = []
    for key, rs in groups.items():
        rec = {"_key": key, "params": {k[2:]: rs[0][k] for k in rs[0] if k.startswith("p_")}}
        m = {}
        for site in ("south", "central", "north"):
            th = site_th.get(site, 2.0)
            for q, lab in (("cv", "CV"), ("ri", "RI"), ("sed", "SED")):
                m[f"{site} {lab}"] = np.nanmean([_f(r, f"site_{site}_{q}_th{th:g}") for r in rs])
        b = np.nanmean([_f(r, "b_fault_below_sat") for r in rs])
        m["b"] = b
        n6 = np.nanmean([_f(r, "n_m6") for r in rs])
        n7 = np.nanmean([_f(r, "n_m7") for r in rs])
        m["excess M>=7"] = (n7 / n6) / 10 ** (-b) if n6 and np.isfinite(b) else np.nan
        m["Bath M>=6"] = np.nanmean([_f(r, "bath_mean") for r in rs])
        m["Bath M>=7"] = np.nanmean([_f(r, "bath_mean_m7") for r in rs])
        m["Omori p"] = np.nanmean([_f(r, "omori_p") for r in rs])
        m["max M"] = np.nanmean([_f(r, "max_mag") for r in rs])
        m["whole-fault CV M>=7"] = np.nanmean([_f(r, "cv_m7_declustered") for r in rs])
        m["co-rupture C-S"] = np.nanmean([_f(r, "site_corupture_central_south") for r in rs])
        rec["metrics"] = m
        out.append(rec)
    return out


def main():
    argv, kw, _ = parse_cli(sys.argv)
    out_prefix = argv[1]
    rows = []
    for i, path in enumerate(argv[2:]):
        rs = load_summary(path)
        for r in rs:
            r["_src"] = i
        rows += rs
    site_th = kw["site_thresholds"]
    designs = design_table(rows, site_th)
    # acceptance (paleo + seismicity, Bath not constrained), per source summary
    acc = set()
    for i in range(len(argv) - 2):
        sc = score([r for r in rows if r["_src"] == i], 2.0, True, site_thresholds=site_th, bath_sigma=1.0)
        acc |= {(i, x["key"]) for x in accepted(sc)}
    for d in designs:
        d["metrics"]["accepted"] = 1.0 if d["_key"] in acc else 0.0
    # predictors: numeric parameters that vary
    # predictors: numeric parameters that vary within at least one sweep (parameters that differ only
    # between sweeps, such as the geometry, would stand in for everything that differs between them)
    names = sorted({k for d in designs for k in d["params"]})
    P = []
    for k in names:
        if k in ("seed",) or k not in PARAM_FIELDS:
            continue
        try:
            v = np.array([float(d["params"].get(k, np.nan)) for d in designs])
        except (TypeError, ValueError):
            continue
        src = np.array([d["_key"][0] for d in designs])
        within = any(np.nanstd(v[src == i]) > 0 for i in np.unique(src) if np.isfinite(v[src == i]).sum() > 1)
        if not within or np.isfinite(v).sum() < 10:
            continue
        logt = (np.nanmin(v) > 0 and np.nanmax(v) / np.nanmin(v) > 5 and "sigma" not in k and "slip_cv" not in k)
        P.append((k, np.log(v) if logt else v, logt, v))
    metrics = [m for m in designs[0]["metrics"] if m != "accepted"
               and np.isfinite([d["metrics"][m] for d in designs]).sum() >= 10]
    print(f"{len(designs)} design vectors, {len(P)} varying parameters, {len(metrics)} metrics, {len(acc)} accepted")

    X = np.column_stack([p[1] for p in P])
    for j in range(X.shape[1]):              # parameters absent in some sweeps: fill with the median
        col = X[:, j]
        col[~np.isfinite(col)] = np.nanmedian(col)
    Z = (X - X.mean(0)) / X.std(0)
    coef = np.full((len(metrics), len(P)), np.nan)
    r2 = np.full(len(metrics), np.nan)
    for i, m in enumerate(metrics):
        y = np.array([d["metrics"][m] for d in designs])
        ok = np.isfinite(y)
        if ok.sum() < len(P) + 5 or np.std(y[ok]) == 0:
            continue
        ys = (y[ok] - y[ok].mean()) / y[ok].std()
        A = np.column_stack([np.ones(ok.sum()), Z[ok]])
        c = np.linalg.lstsq(A, ys, rcond=None)[0]
        coef[i] = c[1:]
        r2[i] = 1 - np.sum((ys - A @ c) ** 2) / np.sum((ys - ys.mean()) ** 2)

    plabels = [("log " if p[2] else "") + p[0] for p in P]
    fig, ax = plt.subplots(figsize=(1.0 * len(P) + 4, 0.5 * len(metrics) + 2))
    im = ax.imshow(coef, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    for i in range(len(metrics)):
        for j in range(len(P)):
            if np.isfinite(coef[i, j]) and abs(coef[i, j]) >= 0.1:
                ax.text(j, i, f"{coef[i, j]:+.2f}", ha="center", va="center", fontsize=7)
    ax.set_xticks(range(len(P)))
    ax.set_xticklabels(plabels, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(metrics)))
    ax.set_yticklabels([f"{m}  (R2 {r2[i]:.2f})" for i, m in enumerate(metrics)], fontsize=8)
    ax.set_title(f"standardized effects ({len(designs)} design vectors): metric change in sd per sd of the parameter")
    plt.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout()
    fig.savefig(out_prefix + "_effects.png", dpi=120)
    plt.close(fig)

    # dependence on the most influential parameters for the key metrics
    km = [m for m in KEY_METRICS if m in metrics or m == "accepted"]
    infl = np.nansum(np.abs(coef[[metrics.index(m) for m in km if m != "accepted"]]), axis=0)
    top = list(np.argsort(-infl)[:4])
    fig, axs = plt.subplots(len(km), len(top), figsize=(4.2 * len(top), 2.6 * len(km)), squeeze=False)
    for jj, j in enumerate(top):
        name, _, logt, raw = P[j]
        x = raw.copy()
        x[~np.isfinite(x)] = np.nanmedian(x)
        edges = np.unique(np.nanpercentile(x, np.linspace(0, 100, 7)))
        idx = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, len(edges) - 2)
        centres = np.array([np.median(x[idx == k]) for k in range(len(edges) - 1)])
        for ii, m in enumerate(km):
            ax = axs[ii, jj]
            y = np.array([d["metrics"][m] for d in designs])
            ax.plot(x, y, ".", color="0.7", ms=3)
            med = [np.nanmedian(y[idx == k]) for k in range(len(edges) - 1)]
            lo = [np.nanpercentile(y[idx == k], 25) for k in range(len(edges) - 1)]
            hi = [np.nanpercentile(y[idx == k], 75) for k in range(len(edges) - 1)]
            if m == "accepted":
                ax.plot(centres, [np.nanmean(y[idx == k]) for k in range(len(edges) - 1)], "o-", color="C3")
                ax.set_ylim(-0.05, 1.05)
            else:
                ax.errorbar(centres, med, yerr=[np.subtract(med, lo), np.subtract(hi, med)], fmt="o-", color="C0",
                            capsize=3)
            if logt:
                ax.set_xscale("log")
            if ii == len(km) - 1:
                ax.set_xlabel(name)
            if jj == 0:
                ax.set_ylabel(m if m != "accepted" else "acceptance rate")
    fig.suptitle("binned medians and quartiles over design vectors (acceptance: fraction per bin)")
    fig.tight_layout()
    fig.savefig(out_prefix + "_dependence.png", dpi=110)
    plt.close(fig)

    with open(out_prefix + "_effects.md", "w") as fh:
        fh.write(f"Standardized regression coefficients over {len(designs)} design vectors (metric sd per "
                 "parameter sd; |c| < 0.1 left blank), and R^2.\n\n")
        fh.write("| metric | R2 | " + " | ".join(plabels) + " |\n")
        fh.write("|" + "---|" * (len(P) + 2) + "\n")
        for i, m in enumerate(metrics):
            fh.write(f"| {m} | {r2[i]:.2f} | " + " | ".join(
                f"{coef[i, j]:+.2f}" if np.isfinite(coef[i, j]) and abs(coef[i, j]) >= 0.1 else ""
                for j in range(len(P))) + " |\n")
    print("wrote", out_prefix + "_effects.png", out_prefix + "_dependence.png", out_prefix + "_effects.md")


if __name__ == "__main__":
    main()
