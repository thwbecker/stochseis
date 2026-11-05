"""How often do synthetic records of the observed length match the observations?

    python3 examples/alpine/record_match.py out_prefix label1 dir1[:filter] [label2 dir2[:filter] ...]
        [--threshold 2] [--site north=1] [--step 100]

For each case (a directory of catalogs, one design vector dir:v0001, or one
realization dir:v0001_s5) synthetic paleoseismic records are sampled the
way the real ones were collected: for an end time t_e (every --step years
through each catalog, after burn-in), the record at each site consists of
the last n events before t_e, with n the number of events in the
observed record (n_ri: 13 south, 14 central, 4 north). All three sites end
at the same time, so correlations between sites (shared ruptures) are
kept. Site records use the detection thresholds (default 2 m; per site
with --site) and the merge interval of the catalogs (site_merge_yr).

Matching, per site:
  RI: |RI - RI_obs| <= 2 x standard error of the observed mean (sigma / sqrt(n - 1))
  CV: |CV - CV_obs| <= 2 x CV_obs / sqrt(2 (n - 1)), the approximate standard error of a CV
      estimated from n - 1 intervals
and jointly: all three sites at the same end time.

Outputs
  <out_prefix>.png  top: RI versus CV of the sampled records at each site with the observation and its
                    matching box; bottom: fractions of sampled record sets that match (RI, CV, both;
                    each site and all three sites jointly)
  <out_prefix>.md   the fractions, and the best-matching record set of each case

Record sets ending 100 yr apart overlap strongly (a record spans about
3,500 yr at the central site), so the number of independent record sets
is about (run length / record length) per realization; use tens of
realizations per case for stable fractions.
"""
import glob
import os
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
warnings.filterwarnings("ignore")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from stochseis import Catalog  # noqa: E402
from stochseis.sites import load_sites, site_records  # noqa: E402

DATA = os.path.join(HERE, "..", "..", "data", "nz_alpine")
SITES = ("south", "central", "north")
COLORS = ["C0", "C3", "C2", "C1", "C4", "C5"]


def parse(argv):
    opt = dict(threshold=2.0, site={}, step=100.0)
    rest = []
    k = 1
    while k < len(argv):
        a = argv[k]
        if a == "--threshold":
            opt["threshold"] = float(argv[k + 1])
            k += 2
        elif a == "--site":
            n, v = argv[k + 1].split("=")
            opt["site"][n] = float(v)
            k += 2
        elif a == "--step":
            opt["step"] = float(argv[k + 1])
            k += 2
        else:
            rest.append(a)
            k += 1
    return rest[0], list(zip(rest[1::2], rest[2::2])), opt


def case_files(spec):
    d, _, v = spec.partition(":")
    if "_s" in v:                                   # one realization: exact match (s1 must not match s10)
        return sorted(glob.glob(os.path.join(d, f"cat_{v}.npz")))
    pat = f"cat_{v}_s*.npz" if v else "cat_v*_s*.npz"
    return sorted(glob.glob(os.path.join(d, pat)))


def tolerances(obs):
    tol = {}
    for n in SITES:
        s = obs[n]
        nri = int(s["n_ri"])
        cv = s["ri_sigma_yr"] / s["ri_yr"]
        tol[n] = dict(n=nri, ri=s["ri_yr"], cv=cv, d_ri=2 * s["ri_sigma_yr"] / np.sqrt(nri - 1),
                      d_cv=2 * cv / np.sqrt(2 * (nri - 1)))
    return tol


def sample_records(path, sites, tol, opt):
    """RI and CV per site for record sets ending every opt['step'] years; returns array (n_end, 3 sites, 2)."""
    cat = Catalog.load(path)
    t0 = 0.1 * cat.duration
    merge = getattr(cat.params, "site_merge_yr", 0.0)
    times = {}
    for s in sites:
        th = opt["site"].get(s["name"], opt["threshold"])
        times[s["name"]] = site_records(cat, [s], threshold=th, t_start=t0, merge_yr=merge)[s["name"]]["times"]
    out = []
    ends = []
    for te in np.arange(t0, cat.duration + 1e-9, opt["step"]):
        rec = []
        for n in SITES:
            t = times[n][times[n] < te]
            if len(t) < tol[n]["n"]:
                rec = None
                break
            dt = np.diff(t[-tol[n]["n"]:])
            rec.append((dt.mean(), dt.std(ddof=1) / dt.mean()))
        if rec is not None:
            out.append(rec)
            ends.append(te)
    return np.array(out).reshape(-1, 3, 2), np.array(ends)


def main():
    out_prefix, pairs, opt = parse(sys.argv)
    sites = [s for s in load_sites(os.path.join(DATA, "paleo_sites.csv")) if s["name"] in SITES]
    obs = {s["name"]: s for s in sites}
    tol = tolerances(obs)
    cases = []
    for label, spec in pairs:
        files = case_files(spec)
        if not files:
            sys.exit(f"no catalogs for {spec}")
        recs, where = [], []
        for p in files:
            r, e = sample_records(p, sites, tol, opt)
            recs.append(r)
            where += [(os.path.basename(p), te) for te in e]
        R = np.concatenate(recs) if recs else np.zeros((0, 3, 2))
        cases.append((label, R, where, len(files)))
        print(f"{label}: {len(files)} catalogs, {len(R)} record sets")

    fig = plt.figure(figsize=(17, 10))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.1, 1])
    results = []
    for k, (label, R, where, nfiles) in enumerate(cases):
        c = COLORS[k % len(COLORS)]
        ri_ok = np.zeros((len(R), 3), bool)
        cv_ok = np.zeros((len(R), 3), bool)
        for j, n in enumerate(SITES):
            ri_ok[:, j] = np.abs(R[:, j, 0] - tol[n]["ri"]) <= tol[n]["d_ri"]
            cv_ok[:, j] = np.abs(R[:, j, 1] - tol[n]["cv"]) <= tol[n]["d_cv"]
            ax = fig.add_subplot(gs[0, j]) if k == 0 else fig.axes[j]
            sub = R[:: max(1, len(R) // 1500)]          # thin for plotting
            ax.plot(sub[:, j, 0], sub[:, j, 1], ".", color=c, ms=2, alpha=0.3)
            ax.plot([], [], "o", color=c, ms=5, label=label)
        both = ri_ok & cv_ok
        # misfit for the best-matching record set: sum of squared normalized residuals
        mis = sum(((R[:, j, 0] - tol[n]["ri"]) / (tol[n]["d_ri"] / 2)) ** 2
                  + ((R[:, j, 1] - tol[n]["cv"]) / (tol[n]["d_cv"] / 2)) ** 2 for j, n in enumerate(SITES))
        b = int(np.argmin(mis)) if len(R) else None
        results.append(dict(label=label, n=len(R), files=nfiles,
                            p_ri=ri_ok.mean(0), p_cv=cv_ok.mean(0), p_both=both.mean(0),
                            j_ri=ri_ok.all(1).mean(), j_cv=cv_ok.all(1).mean(), j_both=both.all(1).mean(),
                            best=(R[b], where[b], mis[b] / 6) if b is not None else None))
    for j, n in enumerate(SITES):
        ax = fig.axes[j]
        t = tol[n]
        ax.add_patch(plt.Rectangle((t["ri"] - t["d_ri"], t["cv"] - t["d_cv"]), 2 * t["d_ri"], 2 * t["d_cv"],
                                   fill=False, ec="k", lw=1.5))
        ax.plot(t["ri"], t["cv"], "ks", ms=7, label="observed (box: match)")
        ax.set(xlabel=f"mean interval of the last {t['n']} events [yr]", ylabel="CV",
               title=f"{n} site: synthetic records of the observed length")
        ax.legend(fontsize=7)
    ax = fig.add_subplot(gs[1, :])
    cats = [f"{n} RI" for n in SITES] + [f"{n} CV" for n in SITES] + [f"{n} both" for n in SITES] \
        + ["all RI", "all CV", "all both"]
    width = 0.8 / len(results)
    for k, r in enumerate(results):
        vals = list(r["p_ri"]) + list(r["p_cv"]) + list(r["p_both"]) + [r["j_ri"], r["j_cv"], r["j_both"]]
        x = np.arange(len(cats)) + (k - (len(results) - 1) / 2) * width
        ax.bar(x, vals, width, color=COLORS[k % len(COLORS)], edgecolor="k", label=r["label"])
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels(cats, rotation=30, ha="right")
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 1.05)
    ax.set(ylabel="fraction of record sets that match",
           title="matching the observations: each site, and all three sites at the same time (all)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_prefix + ".png", dpi=110)
    plt.close(fig)

    with open(out_prefix + ".md", "w") as fh:
        fh.write(f"Synthetic records of the observed length (last n events before an end time, every "
                 f"{opt['step']:g} yr; n = " + ", ".join(f"{n} {tol[n]['n']}" for n in SITES) + "). Match: RI within "
                 + ", ".join(f"{n} +-{tol[n]['d_ri']:.0f} yr" for n in SITES) + "; CV within "
                 + ", ".join(f"{n} +-{tol[n]['d_cv']:.2f}" for n in SITES) + ".\n\n")
        fh.write("| case | record sets | RI match S / C / N | CV match S / C / N | both S / C / N | "
                 "all sites RI | all sites CV | all sites both |\n|---|---|---|---|---|---|---|---|\n")
        for r in results:
            f3 = lambda v: " / ".join(f"{100 * x:.0f}%" for x in v)  # noqa: E731
            fh.write(f"| {r['label']} | {r['n']} | {f3(r['p_ri'])} | {f3(r['p_cv'])} | {f3(r['p_both'])} | "
                     f"{100 * r['j_ri']:.1f}% | {100 * r['j_cv']:.1f}% | {100 * r['j_both']:.2f}% |\n")
        fh.write("\nBest-matching record set per case (RI [yr] / CV at south, central, north; mean squared "
                 "normalized residual):\n\n")
        for r in results:
            if r["best"] is None:
                continue
            R, (f, te), m = r["best"]
            fh.write(f"* {r['label']}: " + ", ".join(f"{n} {R[j, 0]:.0f} / {R[j, 1]:.2f}" for j, n in enumerate(SITES))
                     + f" (misfit {m:.2f}; {f}, ending at {te:.0f} yr)\n")
        fh.write("* observed: " + ", ".join(f"{n} {tol[n]['ri']:.0f} / {tol[n]['cv']:.2f}" for n in SITES) + "\n")
    print("wrote", out_prefix + ".png", out_prefix + ".md")


if __name__ == "__main__":
    main()
