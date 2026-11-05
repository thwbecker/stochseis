"""Tests for the Alpine campaign scripts (scoring and spec generation)."""
import json
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "examples", "nz_alpine"))

from abc_scoring import accepted, score  # noqa: E402


def _row(i, seed, ri, cv, sed, b=1.0, bath=1.2):
    r = {"design_index": i, "p_seed": seed, "p_stress_drop_mpa": 5.0 + i, "b_fault_below_sat": b, "bath_mean": bath}
    for n in ("south", "central", "north"):
        r[f"site_{n}_ri_th2"] = ri[n]
        r[f"site_{n}_cv_th2"] = cv
        r[f"site_{n}_sed_th2"] = sed[n]
    return r


def test_score_ranks_the_data_first():
    good_ri = {"south": 268.0, "central": 249.0, "north": 413.0}
    good_sed = {"south": 7.5, "central": 7.1, "north": 4.0}
    bad_ri = {k: 2 * v for k, v in good_ri.items()}
    rows = [_row(0, 1, bad_ri, 0.8, good_sed), _row(1, 1, good_ri, 0.3, good_sed),
            _row(1, 2, good_ri, 0.3, good_sed), _row(2, 1, good_ri, 0.3, good_sed, b=0.05, bath=0.7)]
    s = score(rows, 2.0, False)
    assert s[0]["key"] in (1, 2) and s[-1]["key"] == 0
    assert len(accepted(s)) == 2
    s2 = score(rows, 2.0, True)
    assert s2[0]["key"] == 1 and len(accepted(s2)) == 1      # unphysical seismicity rejected
    assert np.isclose(s2[0]["detail"]["seismicity"][0], 1.0)


def test_campaign_specs(tmp_path):
    cfg = {"configs": [{"name": "a_fine", "geometry": "a", "resolution": "fine", "geometry_file": __file__,
                        "patch_km": 1.0, "segment_max_dist_km": 10.0},
                       {"name": "a_coarse", "geometry": "a", "resolution": "coarse", "geometry_file": __file__,
                        "patch_km": 2.0, "segment_max_dist_km": 40.0}]}
    (tmp_path / "cfg.json").write_text(json.dumps(cfg))
    base = os.path.join(ROOT, "examples", "nz_alpine", "abc5_planar_l1.json")
    subprocess.run([sys.executable, os.path.join(ROOT, "examples", "nz_alpine", "make_campaign_specs.py"),
                    str(tmp_path / "cfg.json"), base, str(tmp_path / "specs"), "7", "2", "500"], check=True,
                   capture_output=True)
    a = json.loads((tmp_path / "specs" / "abc_a_fine.json").read_text())
    b = json.loads((tmp_path / "specs" / "abc_a_coarse.json").read_text())
    assert a["base"]["patch_km"] == 1.0 and b["base"]["patch_km"] == 2.0
    assert a["n_samples"] == 7 and a["seeds"] == [1, 2] and a["base"]["duration_yr"] == 500
    assert a["design_seed"] == b["design_seed"] and a["ranges"] == b["ranges"]       # paired designs
    assert a["base"]["dip_correction"] is True
    assert a["base"]["store_slip_min_mag"] == 5.5 and b["base"]["store_slip_min_mag"] == 6.0


def test_site_specific_threshold():
    good_ri = {"south": 268.0, "central": 249.0, "north": 413.0}
    good_sed = {"south": 7.5, "central": 7.1, "north": 4.0}
    r = _row(0, 1, good_ri, 0.3, good_sed)
    # at 2 m the north site is too slow, at 1 m it fits
    r["site_north_ri_th2"] = 800.0
    for n in ("south", "central", "north"):
        for q in ("ri", "cv", "sed"):
            r[f"site_{n}_{q}_th1"] = r[f"site_{n}_{q}_th2"]
    r["site_north_ri_th1"] = 413.0
    m2 = score([r], 2.0, False)[0]["misfit"]
    m1 = score([r], 2.0, False, site_thresholds={"north": 1.0})[0]["misfit"]
    assert m1 < 0.5 < m2


def test_bath_options():
    good_ri = {"south": 268.0, "central": 249.0, "north": 413.0}
    good_sed = {"south": 7.5, "central": 7.1, "north": 4.0}
    r = _row(0, 1, good_ri, 0.3, good_sed, bath=1.5)
    r["bath_mean_m7"] = 1.2
    loose = score([r], 2.0, True)[0]["misfit"]
    tight = score([r], 2.0, True, bath_sigma=0.05)[0]["misfit"]
    m7 = score([r], 2.0, True, bath_sigma=0.05, bath_key="bath_mean_m7")[0]["misfit"]
    assert m7 < loose < tight
    from abc_scoring import parse_cli
    rest, kw, seis = parse_cli(["x", "a.json", "--seismicity", "--site", "north=1", "--bath-sigma", "0.1", "--bath-m7"])
    assert rest == ["x", "a.json"] and seis and kw["site_thresholds"] == {"north": 1.0}
    assert kw["bath_sigma"] == 0.1 and kw["bath_key"] == "bath_mean_m7"


def test_realization_specs_and_windows(tmp_path):
    cfg = json.load(open(os.path.join(ROOT, "examples", "nz_alpine", "realization_sets.json")))
    cfg["configuration"]["geometry_file"] = __file__          # any existing file
    cfg["base_spec"] = os.path.join(ROOT, "examples", "nz_alpine", "abc7_regularity.json")
    (tmp_path / "sets.json").write_text(json.dumps(cfg))
    subprocess.run([sys.executable, os.path.join(ROOT, "examples", "nz_alpine", "make_realization_specs.py"),
                    str(tmp_path / "sets.json"), str(tmp_path / "specs"), "4", "1000"], check=True, capture_output=True)
    for name, st in cfg["sets"].items():
        spec = json.loads((tmp_path / "specs" / f"realizations_{name}.json").read_text())
        assert [v[st["parameter"]] for v in spec["vectors"]] == st["values"]
        assert spec["seeds"] == list(range(st["seed0"], st["seed0"] + 4))
        for v in spec["vectors"]:
            for k, x in cfg["reference"].items():
                if k != st["parameter"]:
                    assert v[k] == x
            for k, x in st["fixed"].items():
                assert v[k] == x
        cases = (tmp_path / "specs" / f"realizations_{name}_cases.txt").read_text().split("\n")
        assert cases[0].split() == [st["labels"][0], "v0000", str(st["seed0"])]
    from nz_page import windows
    w = windows(np.array([100.0, 200, 300, 100, 200, 300, 50]), 3)
    assert len(w) == 2 and np.isclose(w[0][0], 200.0) and np.isclose(w[0][1], 0.5)


def test_parameter_effects_runs(tmp_path):
    rng = np.random.default_rng(0)
    rows = []
    for i in range(30):
        area = float(np.exp(rng.uniform(np.log(100), np.log(3000))))
        gamma = float(rng.uniform(2, 12))
        for seed in (1, 2):
            r = _row(i, seed, {"south": 260.0, "central": 250.0, "north": 420.0}, 0.1 + 0.1 * np.log(area) + 0.01 * rng.normal(),
                     {"south": 7.0, "central": 7.0, "north": 4.0})
            r.update(p_hazard_max_area_km2=area, p_nucleation_exponent=gamma, n_m6=100.0, n_m7=12.0, omori_p=1.0,
                     max_mag=8.2, cv_m7_declustered=0.5, site_corupture_central_south=0.7)
            rows.append(r)
    f = tmp_path / "s.json"
    f.write_text(json.dumps(rows))
    subprocess.run([sys.executable, os.path.join(ROOT, "examples", "nz_alpine", "parameter_effects.py"), str(tmp_path / "pe"),
                    str(f)], check=True, capture_output=True)
    table = (tmp_path / "pe_effects.md").read_text()
    line = [l for l in table.split("\n") if l.startswith("| central CV")][0]
    cells = [c.strip() for c in line.split("|")]
    header = [c.strip() for c in table.split("\n")[2].split("|")]
    assert float(cells[header.index("log hazard_max_area_km2")]) > 0.9      # CV built to depend on log area only
    assert (tmp_path / "pe_effects.png").exists() and (tmp_path / "pe_dependence.png").exists()


def test_record_match_tolerances():
    from record_match import tolerances
    from abc_scoring import load_sites as load_site_rows

    obs = {r["name"]: {k: (float(v) if k not in ("name", "sense") else v) for k, v in r.items()} for r in load_site_rows()}
    tol = tolerances(obs)
    assert tol["central"]["n"] == 14 and tol["north"]["n"] == 4
    assert abs(tol["central"]["d_ri"] - 2 * 57 / np.sqrt(13)) < 1e-9
    assert abs(tol["south"]["d_cv"] - 2 * (108 / 268) / np.sqrt(24)) < 1e-9


def test_case_files_exact_realization(tmp_path):
    for s in (1, 2, 10, 11):
        (tmp_path / f"cat_v0000_s{s}.npz").write_bytes(b"")
    (tmp_path / "cat_v0001_s1.npz").write_bytes(b"")
    import overview_figures
    import record_match
    import site_analysis
    import spectrum_figures
    for mod in (overview_figures, site_analysis, spectrum_figures, record_match):
        assert [os.path.basename(f) for f in mod.case_files(f"{tmp_path}:v0000_s1")] == ["cat_v0000_s1.npz"]
        assert len(mod.case_files(f"{tmp_path}:v0000")) == 4
        assert len(mod.case_files(str(tmp_path))) == 5


def test_weighted_quantile():
    from predictions import wquantile
    x = np.arange(10.0)
    assert abs(wquantile(x, np.ones(10), 0.5) - 4.5) < 0.6
    w = np.zeros(10)
    w[7] = 1.0
    assert wquantile(x, w, 0.05) == 7.0 and wquantile(x, w, 0.95) == 7.0
