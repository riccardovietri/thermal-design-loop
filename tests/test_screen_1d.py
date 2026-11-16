import math
from dataclasses import replace

import numpy as np
import pytest

from tools.pareto import (capacity_ratios, knee_index, marginal_return, pareto_indices, recommend, recommend_marginal_return,
                          recommend_marginal_return_objective, cutoff_sensitivity, stable_cutoff_range, local_returns,
                          CUTOFF_SWEEP,
                          spread_along_front)
from tools.screen_1d import analyze_1d, f_re_fanning, nu_h1, screen_1d
from tools.types import Design, DesignSpace, Fluid, Requirements

REQ = Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, 40.0, 30e3, "PG25")
FLUID = Fluid("test", 40.0, 1000.0, 4000.0, 0.5, 1.0e-3)
SPACE = DesignSpace((0.15e-3, 0.5e-3), (0.15e-3, 0.5e-3), (1e-3, 4e-3), (1.0, 2.0), 0.15e-3, 10.0)

# Shah & London (1978) tabulated values: alpha -> (Fanning f*Re, Nu_H1)
TABLE = {1.0: (14.227, 3.608), 0.5: (15.548, 4.123), 0.25: (18.233, 5.331), 0.125: (20.585, 6.490)}


@pytest.mark.parametrize("alpha", TABLE)
def test_polynomials_match_shah_london_table(alpha):
    assert f_re_fanning(alpha) == pytest.approx(TABLE[alpha][0], rel=0.005)
    assert nu_h1(alpha) == pytest.approx(TABLE[alpha][1], rel=0.005)


def test_hand_calculation():
    """0.5 mm channels and fins, 2 mm tall, 1.5 LPM: every term worked by hand from the definitions."""
    r = analyze_1d(Design(0.5e-3, 0.5e-3, 2e-3, 1.5), REQ, FLUID)
    assert r.n_channels == 60
    v = 2.5e-5 / (60 * 0.5e-3 * 2e-3)  # 0.41667 m/s
    dh = 0.8e-3
    assert r.reynolds == pytest.approx(1000 * v * dh / 1e-3)
    assert r.dp_pa == pytest.approx(2 * f_re_fanning(0.25) * 1e-3 * v * 0.06 / dh**2)
    htc = nu_h1(0.25) * 0.5 / dh
    m = math.sqrt(2 * htc / (390 * 0.5e-3))
    eta = math.tanh(m * 2e-3) / (m * 2e-3)
    assert r.fin_efficiency == pytest.approx(eta)
    assert r.r_conv == pytest.approx(1 / (60 * 0.06 * htc * (0.5e-3 + 2 * eta * 2e-3)))
    assert r.r_base == pytest.approx(0.002 / (390 * 0.06 * 0.06))
    assert r.r_caloric == pytest.approx(1 / (2 * 1000 * 4000 * 2.5e-5))
    assert r.r_th_k_per_w == pytest.approx(r.r_base + r.r_conv + r.r_caloric)


def test_trends():
    base = analyze_1d(Design(0.3e-3, 0.3e-3, 2e-3, 1.0), REQ, FLUID)
    more_flow = analyze_1d(Design(0.3e-3, 0.3e-3, 2e-3, 2.0), REQ, FLUID)
    assert more_flow.dp_pa == pytest.approx(2 * base.dp_pa)
    assert more_flow.r_th_k_per_w < base.r_th_k_per_w
    assert more_flow.r_conv == pytest.approx(base.r_conv)  # fully developed laminar: h independent of flow


def test_constraints_flagged():
    tight = replace(REQ, dp_budget_pa=100.0)
    r = analyze_1d(Design(0.15e-3, 0.15e-3, 4e-3, 2.0), tight, FLUID, SPACE)
    assert not r.feasible
    assert {"dp_over_budget", "fin_aspect_ratio", "channel_aspect_ratio"} <= set(r.violations)


def test_sweep_is_full_factorial():
    assert len(screen_1d(REQ, FLUID, SPACE)) == 16


def test_pareto_front_is_non_dominated():
    x = [1, 2, 3, 4, 2.5]
    y = [5, 3, 4, 1, 2.9]
    front = pareto_indices(x, y)
    assert front == [0, 1, 4, 3]
    for i in front:
        assert not any(x[j] <= x[i] and y[j] <= y[i] and (x[j] < x[i] or y[j] < y[i]) for j in range(len(x)))


def test_spread_and_knee():
    xs = [10**e for e in (0, 0.1, 0.2, 1.0, 2.0, 2.1, 3.0)]
    picks = spread_along_front(xs, 4)
    assert picks[0] == 0 and picks[-1] == 6 and len(picks) == 4
    assert knee_index([1, 10, 100], [100, 2, 1]) == 1


# (dp_pa, r_th_k_per_w) candidates for the recommendation rule; budget 6 kPa, target 15 mK/W
CANDIDATES = [(500.0, 30e-3), (1500.0, 14.8e-3), (3000.0, 12e-3), (5500.0, 10e-3), (9000.0, 8e-3)]


def test_recommend_target_met_picks_lowest_dp_among_qualifying():
    rec = recommend(CANDIDATES, 6e3, 15e-3)
    assert rec == {"index": 1, "rule": "meets_target_lowest_dp"}


def test_recommend_falls_back_to_lowest_r_th_when_target_unmet():
    rec = recommend(CANDIDATES, 6e3, 5e-3)
    assert rec == {"index": 3, "rule": "fallback_lowest_r_th"}


def test_recommend_budget_excludes_lower_r_th_design():
    # the 8 mK/W design is the only one meeting a 9 mK/W target, but it is over budget
    rec = recommend(CANDIDATES, 6e3, 9e-3)
    assert rec == {"index": 3, "rule": "fallback_lowest_r_th"}
    assert recommend(CANDIDATES, 10e3, 9e-3) == {"index": 4, "rule": "meets_target_lowest_dp"}
    with pytest.raises(ValueError):
        recommend(CANDIDATES, 100.0, 15e-3)


# Marginal-return rule (SPEC Amendment 3): budget 5 kPa, threshold 20 mK/W, objective 15 mK/W, cut-off 0.25


def test_marginal_return_definition():
    # doubling ΔP for 16% lower R_th is a return of about 0.25
    assert marginal_return(1000.0, 20e-3, 2000.0, 20e-3 * 2 ** -0.25) == pytest.approx(0.25)
    assert marginal_return(1000.0, 20e-3, 4000.0, 10e-3) == pytest.approx(0.5)


def test_marginal_return_stops_where_return_drops_below_cutoff():
    # returns 0->1: ln(18/14)/ln(2) = 0.36; 1->2: ln(14/13)/ln(2) = 0.11, 1->3: ln(14/12)/ln(4) = 0.11
    cands = [(500.0, 18e-3), (1000.0, 14e-3), (2000.0, 13e-3), (4000.0, 12e-3)]
    rec = recommend_marginal_return(cands, 5e3, 20e-3, 15e-3, 0.25)
    assert rec["index"] == 1 and rec["rule"] == "marginal_return"
    assert rec["meets_threshold"] and rec["meets_objective"]
    assert [s["from"] for s in rec["steps"]] == [0, 1]
    assert rec["steps"][0]["to"] == 1 and rec["steps"][1]["to"] is None
    assert rec["steps"][0]["returns"][0]["return"] == pytest.approx(np.log(18 / 14) / np.log(2))
    assert all(e["return"] < 0.25 for e in rec["steps"][1]["returns"])


def test_marginal_return_skips_flat_step_to_large_gain():
    # 0->1 is nearly flat (0.03), but 0->2 doubles ΔP twice for half the R_th (0.5): move to 2, skipping 1
    cands = [(1000.0, 19e-3), (2000.0, 18.6e-3), (4000.0, 9.5e-3)]
    rec = recommend_marginal_return(cands, 5e3, 20e-3, 15e-3, 0.25)
    assert rec["index"] == 2 and rec["rule"] == "marginal_return"
    assert rec["steps"][0]["returns"][0]["return"] < 0.25
    assert rec["steps"][0]["to"] == 2
    assert rec["meets_objective"]


def test_marginal_return_moves_to_nearest_qualifying_point():
    # both 1 and 2 qualify from 0; the rule moves to the nearer one (1) first
    cands = [(1000.0, 19e-3), (2000.0, 15e-3), (4000.0, 10e-3)]
    rec = recommend_marginal_return(cands, 5e3, 20e-3, 15e-3, 0.25)
    assert rec["steps"][0]["to"] == 1
    assert rec["index"] == 2  # 1->2: ln(1.5)/ln(2) = 0.58


def test_marginal_return_fallback_when_threshold_unmet():
    cands = [(500.0, 30e-3), (2000.0, 24e-3), (4500.0, 21e-3), (9000.0, 15e-3)]
    rec = recommend_marginal_return(cands, 5e3, 20e-3, 15e-3, 0.25)
    assert rec["index"] == 2 and rec["rule"] == "fallback_lowest_r_th"
    assert not rec["meets_threshold"] and not rec["meets_objective"]
    assert rec["steps"] == [] and rec["front"] == []
    with pytest.raises(ValueError):
        recommend_marginal_return(cands, 100.0, 20e-3, 15e-3, 0.25)


def test_marginal_return_budget_excludes_better_design():
    # the 9 kPa design would be a large gain from 2 kPa, but it is over budget
    cands = [(1000.0, 18e-3), (2000.0, 17.5e-3), (9000.0, 8e-3)]
    rec = recommend_marginal_return(cands, 5e3, 20e-3, 15e-3, 0.25)
    assert rec["front"] == [0, 1]
    assert rec["index"] == 0 and not rec["meets_objective"]
    assert recommend_marginal_return(cands, 10e3, 20e-3, 15e-3, 0.25)["index"] == 2


def test_marginal_return_single_candidate_and_dominated_points():
    rec = recommend_marginal_return([(3000.0, 16e-3), (4000.0, 17e-3), (6000.0, 10e-3)], 5e3, 20e-3, 15e-3, 0.25)
    assert rec["front"] == [0] and rec["index"] == 0
    assert rec["steps"] == [{"from": 0, "returns": [], "to": None}]
    assert rec["meets_threshold"] and not rec["meets_objective"]


# Amendment 4: objective clause and minimum-ΔP-ratio guard (recommend_marginal_return_objective)
OFFSET_STRIP_LIKE = [(500.0, 18.36e-3), (790.0, 16.55e-3), (990.0, 15.86e-3), (1050.0, 15.74e-3), (1570.0, 14.58e-3),
                     (2290.0, 13.53e-3)]


def test_objective_clause_fires_when_walk_stops_short():
    # every return from the first point is 0.19-0.23 < 0.25, so the walk stops at 18.36 mK/W; the lowest-ΔP point meeting
    # 15 mK/W further along is 14.58 mK/W at 1.57 kPa
    rec = recommend_marginal_return_objective(OFFSET_STRIP_LIKE, 5e3, 20e-3, 15e-3, 0.25)
    assert all(e["return"] < 0.25 for e in rec["steps"][0]["returns"])
    assert rec["walk_stop"] == 0 and rec["index"] == 4
    assert rec["rule"] == "marginal_return_then_objective"
    assert rec["objective_step"] == {"from": 0, "to": 4}
    assert rec["meets_objective"] and not rec["objective_unreachable"]
    assert recommend_marginal_return(OFFSET_STRIP_LIKE, 5e3, 20e-3, 15e-3, 0.25)["index"] == 0


def test_objective_clause_does_not_fire_when_stop_meets_objective():
    cands = [(500.0, 18e-3), (1000.0, 14e-3), (2000.0, 13e-3), (4000.0, 12e-3)]
    rec = recommend_marginal_return_objective(cands, 5e3, 20e-3, 15e-3, 0.25)
    assert rec["index"] == 1 and rec["walk_stop"] == 1 and rec["rule"] == "marginal_return"
    assert rec["objective_step"] is None and rec["meets_objective"]


def test_objective_clause_no_later_design_meets_objective():
    # neither point meets 15 mK/W: stay at the stop point and flag it
    cands = [(1000.0, 18e-3), (2000.0, 17.5e-3), (9000.0, 8e-3)]
    rec = recommend_marginal_return_objective(cands, 5e3, 20e-3, 15e-3, 0.25)
    assert rec["index"] == 0 and rec["rule"] == "marginal_return"
    assert rec["objective_unreachable"] and not rec["meets_objective"] and rec["objective_step"] is None


def test_guard_skips_nearly_coincident_point():
    # 0 -> 1 is a 5% ΔP step for 11% lower R_th, a return of 2.3; the guard (ΔP ratio < 1.10) ignores it, and 0 -> 2 returns
    # ln(19/15.5)/ln(3) = 0.185, so the walk stops at 0. Without the guard it would move to 1.
    cands = [(1000.0, 19e-3), (1050.0, 17e-3), (3000.0, 15.5e-3)]
    rec = recommend_marginal_return_objective(cands, 5e3, 20e-3, 15e-3, 0.25, min_dp_ratio=1.10)
    first = rec["steps"][0]["returns"]
    assert first[0]["skipped"] and first[0]["return"] == pytest.approx(np.log(19 / 17) / np.log(1.05))
    assert not first[1]["skipped"] and first[1]["return"] == pytest.approx(np.log(19 / 15.5) / np.log(3))
    assert rec["walk_stop"] == 0 and rec["index"] == 0 and rec["objective_unreachable"]
    assert recommend_marginal_return_objective(cands, 5e3, 20e-3, 15e-3, 0.25, min_dp_ratio=None)["walk_stop"] == 1


@pytest.mark.parametrize("cands", [
    [(500.0, 18e-3), (1000.0, 14e-3), (2000.0, 13e-3), (4000.0, 12e-3)],
    [(1000.0, 19e-3), (2000.0, 18.6e-3), (4000.0, 9.5e-3)],
    [(1000.0, 19e-3), (2000.0, 15e-3), (4000.0, 10e-3)],
    [(500.0, 30e-3), (2000.0, 24e-3), (4500.0, 21e-3), (9000.0, 15e-3)],
    [(1000.0, 18e-3), (2000.0, 17.5e-3), (9000.0, 8e-3)],
    [(3000.0, 16e-3), (4000.0, 17e-3), (6000.0, 10e-3)],
    [(1000.0, 19e-3), (1050.0, 17e-3), (3000.0, 15.5e-3)],
    OFFSET_STRIP_LIKE,
])
def test_guard_and_clause_off_reproduce_amendment_3_rule(cands):
    old = recommend_marginal_return(cands, 5e3, 20e-3, 15e-3, 0.25)
    new = recommend_marginal_return_objective(cands, 5e3, 20e-3, 15e-3, 0.25, min_dp_ratio=None, objective_clause=False)
    strip = lambda steps: [{**st, "returns": [{"to": e["to"], "return": e["return"]} for e in st["returns"]]} for st in steps]
    assert {k: new[k] for k in old if k != "steps"} == {k: v for k, v in old.items() if k != "steps"}
    assert strip(new["steps"]) == old["steps"]


def test_cutoff_sensitivity_and_stable_range():
    # 0->1 returns ln(18/14)/ln(2) = 0.363; 1->2 0.107, 1->3 0.111, 2->3 0.115; at 0.40 the clause continues to 1
    cands = [(500.0, 18e-3), (1000.0, 14e-3), (2000.0, 13e-3), (4000.0, 12e-3)]
    rows = cutoff_sensitivity(cands, 5e3, 20e-3, 15e-3, (0.05, 0.10, 0.15, 0.25, 0.35, 0.40))
    assert [r["index"] for r in rows] == [3, 3, 1, 1, 1, 1]
    assert rows[-1]["rule"] == "marginal_return_then_objective" and rows[-1]["walk_stop"] == 0
    assert stable_cutoff_range(rows, 0.25) == (0.15, 0.40)
    assert stable_cutoff_range(rows, 0.30) is None
    assert CUTOFF_SWEEP[0] == 0.05 and CUTOFF_SWEEP[-1] == 0.6 and len(CUTOFF_SWEEP) == 12


def test_local_returns():
    lr = local_returns([1000.0, 2000.0, 4000.0], [20e-3, 15e-3, 14e-3])
    assert [e["dp_ratio"] for e in lr] == [2.0, 2.0]
    assert lr[0]["return"] == pytest.approx(np.log(20 / 15) / np.log(2))


def test_capacity_ratios_hand_calculation():
    # 1 kW, R_th 12.5 mK/W, threshold 20 mK/W -> 20 K / 12.5 mK/W = 1.6 kW; ΔP 2 kPa -> 0.8 kW/kPa;
    # 1.5 LPM = 2.5e-5 m3/s, pumping 2000 × 2.5e-5 = 0.05 W -> 1 kW / 0.05 W = 20 kW/W
    r = capacity_ratios(2000.0, 12.5e-3, 1000.0, 1.5e-3 / 60, 20e-3)
    assert r["capacity_kw"] == pytest.approx(1.6)
    assert r["capacity_kw_per_kpa"] == pytest.approx(0.8)
    assert r["kw_per_pumping_w"] == pytest.approx(20.0)
