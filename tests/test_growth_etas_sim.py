import numpy as np
from scipy import stats

from stochseis import Params, Simulation
from stochseis import etas
from stochseis.geometry import planar_fault
from stochseis.rupture import grow_rupture


def test_growth_size_distribution_is_power_law_with_b():
    """Without state arrest, N(>A) ~ A^-b on a fault much larger than the ruptures."""
    p = Params(fault_length_km=120, fault_width_km=40, patch_km=1.0, top_depth_km=3.0,
               arrest_scale=0.0, b_value=1.0)
    fault = planar_fault(p)
    graph = fault.graph()
    rng = np.random.default_rng(0)
    i0 = int(np.argmin((fault.col - 60) ** 2 + (fault.row - 20) ** 2))
    n = []
    for _ in range(4000):
        idx, dist, code, rho, asc = grow_rupture(fault, graph, i0, lambda ii: np.ones(len(ii)),
                                                 3e6, rng, p, 40e3)
        assert 1.0 <= asc < (1.0 + 1.0 / len(idx)) ** 1.5 + 1e-6
        n.append(len(idx))
    n = np.array(n)
    # P(N >= k) = 1/k exactly for b = 1
    for k in (2, 4, 8, 16):
        frac = np.mean(n >= k)
        assert abs(frac - 1.0 / k) < 3.0 * np.sqrt((1.0 / k) * (1 - 1.0 / k) / len(n)) + 0.005
    # distances are non-decreasing along the growth order
    assert np.all(np.diff(dist) >= -1e-6)


def test_state_arrest_blocks_depleted_region():
    """A zero-deficit region stops ruptures that would otherwise cross it (b=0.5 for statistics)."""
    p = Params(fault_length_km=120, fault_width_km=15, patch_km=1.0, arrest_scale=0.5, b_value=0.5)
    fault = planar_fault(p)
    graph = fault.graph()
    rng = np.random.default_rng(1)
    i0 = int(np.argmin((fault.col - 60) ** 2 + (fault.row - 7) ** 2))
    barrier = np.where(fault.col >= 70, 0.0, 5.0)
    uniform = np.full(fault.n, 5.0)
    reach = {}
    for name, deficit in (("barrier", barrier), ("uniform", uniform)):
        cnt = 0
        for _ in range(1500):
            idx, dist, code, rho, asc = grow_rupture(fault, graph, i0, lambda ii: deficit[ii], 3e6, rng, p, 30e3)
            if np.any(fault.col[idx] >= 74):
                cnt += 1
        reach[name] = cnt
    assert reach["uniform"] > 40
    assert reach["barrier"] < 0.1 * reach["uniform"]


def test_omori_times_follow_cdf():
    rng = np.random.default_rng(2)
    c, pexp = 0.01, 1.2
    t = etas.omori_times(20000, c, pexp, rng)
    d, pval = stats.kstest(t, lambda x: etas.omori_cdf(x, c, pexp))
    assert pval > 0.01


def test_gr_magnitudes_slope():
    rng = np.random.default_rng(3)
    m = etas.gr_magnitudes(50000, 1.0, 3.0, 6.0, rng)
    b = np.log10(np.e) / (m.mean() - 3.0)
    assert abs(b - 1.0) < 0.1     # truncation at 6 biases slightly high


def test_productivity_constant_alpha_equal_b():
    e = etas.gr_expected_productivity(1.0, 1.0, 3.0, 7.0)
    assert abs(e - 1.0 * np.log(10) * 4.0 / (1 - 1e-4)) < 1e-2


def test_short_simulation_reproducible_and_consistent():
    p = Params(duration_yr=150, mag_min=3.5, pilot_events=200, seed=7, fault_length_km=100)
    c1 = Simulation(p).run()
    c2 = Simulation(p).run()
    assert len(c1) == len(c2)
    assert np.allclose(c1.events["mag"], c2.events["mag"])
    ev = c1.events
    assert np.all(np.diff(ev["t"]) >= 0)
    # deficit-capped events (nucleations in relaxed zones) fall below mag_min; see HANDOFF.md item 2
    assert np.mean(ev["mag"] >= p.mag_min - 1e-9) > 0.7
    fe = c1.fault_events()
    # moment bookkeeping: released slip equals the sum of all event moments / (mu A)
    tot = np.sum(ev["m0"]) / (c1.fault.mu[0] * c1.fault.area[0])
    assert abs(c1.released.sum() - tot) / tot < 1e-6
    # onset times: duration = max distance / rupture speed, positive for multi-patch events
    multi = fe[fe["npatch"] > 1]
    assert np.all(multi["duration"] > 0)
    assert np.all(multi["duration"] * multi["vr"] <= multi["length"] + multi["width"] + 1.0)
    # children never precede parents
    par = ev["parent"]
    ok = par >= 0
    assert np.all(ev["t"][ok] >= ev["t"][par[ok]])
    # fraction of grown ruptures among all nucleations consistent with p_fault
    frac = np.mean(ev["kind"] == 1)
    assert abs(frac - c1.meta["p_fault"]) < 0.05
    # magnitudes of sub-patch events stay below those of large ruptures and above mag_min
    assert ev["mag"][ev["kind"] == 0].max() < c1.meta["m_patch"] + 0.5
    # save / load round trip
    c1.save("/tmp/_stochseis_test.npz")
    from stochseis import Catalog
    c3 = Catalog.load("/tmp/_stochseis_test.npz")
    assert len(c3) == len(c1)
    assert c3.fault.n == c1.fault.n
    s1 = c1.slip_of(int(c1.slip_event[0])) if len(c1.slip_event) else None
    if s1 is not None:
        s3 = c3.slip_of(int(c3.slip_event[0]))
        assert np.allclose(s1[1], s3[1])


def test_segmented_fault_connections():
    p = Params(fault_length_km=60, fault_width_km=10, segment_gaps_km=[[30, 2.0]], jump_max_km=5.0)
    fault = planar_fault(p)
    assert fault.n_conn == 1
    assert len(np.unique(fault.section)) == 2
    g_open = fault.graph(np.array([True]))
    g_closed = fault.graph(np.array([False]))
    from scipy.sparse.csgraph import connected_components
    assert connected_components(g_open, directed=False)[0] == 1
    assert connected_components(g_closed, directed=False)[0] == 2
    # nearest patch snaps gap positions onto a real patch
    idx = fault.nearest_patch(np.array([30.5e3, 31.0e3]), np.array([5e3, 5e3]))
    assert np.all(idx >= 0)


def test_analytic_calibration_matches_monte_carlo():
    """Analytic expected moment, productivity and grown fraction agree with a large pilot."""
    from stochseis.catalog import moment_magnitude

    p = Params(fault_length_km=60, fault_width_km=15, pilot_events=100, seed=11, duration_yr=1)
    sim = Simulation(p)
    m = sim.meta
    m0, kind, npatch, ds = sim._pilot(n=20000)
    se = m0.std() / np.sqrt(len(m0))
    assert abs(m0.mean() - m["expected_m0_per_nucleation"]) < 4.0 * se
    assert abs(kind.mean() - m["expected_grown_fraction"]) < 0.01
    prod = 10.0 ** (sim.alpha * (moment_magnitude(m0) - p.mag_min))
    assert abs(prod.mean() - m["expected_productivity"]) < 4.0 * prod.std() / np.sqrt(len(prod))


def test_deficit_weighted_nucleation_run():
    p = Params(duration_yr=200, mag_min=3.5, seed=5, fault_length_km=80,
               nucleation_weighting="deficit", cap_by_deficit=False)
    cat = Simulation(p).run()
    ev = cat.events
    assert len(ev) > 100
    # nucleations happen only on loaded patches, so no event is cut below mag_min
    assert np.all(ev["mag"] >= p.mag_min - 1e-9)
    assert cat.meta["background_bound_violations"] < 0.01
    # triggered events keep their number: realised branching ratio near the target
    assert abs(np.mean(ev["parent"] >= 0) - p.branching_ratio) < 0.2


def test_state_only_termination_and_nucleation_exponent():
    """size_hazard=False stops ruptures only by state; gamma > 1 runs and closes the budget roughly."""
    p = Params(duration_yr=300, mag_min=3.5, seed=4, fault_length_km=60, etas=False,
               nucleation_weighting="uniform", rate_feedback=False, size_hazard=False,
               budget_scale=0.1)
    cat = Simulation(p).run()
    fe = cat.fault_events()
    assert len(fe) > 10
    assert np.all(fe["stop"] != 0)        # no size-hazard stops
    p2 = Params(duration_yr=300, mag_min=3.5, seed=4, fault_length_km=60, nucleation_exponent=3.0)
    c2 = Simulation(p2).run()
    assert c2.meta["background_bound_violations"] < 0.01
    assert np.all(c2.events["mag"] >= p2.mag_min - 1e-9)


def test_defaults_produce_gutenberg_richter():
    """Regression: with default parameters the FMD follows GR across the patch magnitude."""
    from stochseis.analysis import b_value_mle

    p = Params(duration_yr=800, mag_min=3.5, seed=2, fault_length_km=100, pilot_events=200)
    cat = Simulation(p).run()
    m = cat.events["mag"]
    m = m[m >= p.mag_min]
    for M in (4.0, 4.5, 5.0, 5.5):
        ratio = np.mean(m >= M) / 10.0 ** (-p.b_value * (M - p.mag_min))
        assert 0.75 < ratio < 1.3, (M, ratio)
    b, sb, n = b_value_mle(m, cat.meta["m_patch"], 6.0)
    assert abs(b - p.b_value) < 3.0 * sb + 0.05
    assert abs(cat.meta["moment_closure"] - 1.0) < 0.25
