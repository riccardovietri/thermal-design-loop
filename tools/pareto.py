"""pareto: non-dominated filtering and shortlist selection for two minimised objectives."""
import numpy as np


def pareto_indices(x, y) -> list[int]:
    """Indices of non-dominated points when minimising both x and y, sorted by x."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    order = np.lexsort((y, x))
    front, best_y = [], np.inf
    for i in order:
        if y[i] < best_y:
            front.append(int(i))
            best_y = y[i]
    return front


def spread_along_front(x_front, n: int) -> list[int]:
    """Pick up to n positions along a front (sorted by x), evenly spaced in log10(x), endpoints included."""
    logx = np.log10(np.asarray(x_front, float))
    if len(logx) <= n:
        return list(range(len(logx)))
    targets = np.linspace(logx[0], logx[-1], n)
    picks = []
    for target in targets:
        for j in np.argsort(np.abs(logx - target)):
            if int(j) not in picks:
                picks.append(int(j))
                break
    return sorted(picks)


def knee_index(x_front, y_front) -> int:
    """Knee of a front: the point farthest from the chord between its endpoints, in log-log space."""
    lx, ly = np.log10(np.asarray(x_front, float)), np.log10(np.asarray(y_front, float))
    if len(lx) < 3:
        return 0
    dx, dy = lx[-1] - lx[0], ly[-1] - ly[0]
    dist = np.abs(dx * (ly - ly[0]) - dy * (lx - lx[0])) / np.hypot(dx, dy)
    return int(np.argmax(dist))


def recommend(candidates, dp_budget_pa: float, r_th_target: float | None) -> dict:
    """Recommendation rule (SPEC Amendment 1) over (dp_pa, r_th_k_per_w) pairs, e.g. CFD-verified designs.

    Among candidates with dp <= budget and r_th <= target, pick the lowest dp ("meets_target_lowest_dp").
    If none meets the target (or there is no target), pick the lowest r_th within budget ("fallback_lowest_r_th").
    Ties break on the other objective. Returns {"index": position in candidates, "rule": name}.
    """
    pts = [(float(dp), float(r)) for dp, r in candidates]
    in_budget = [i for i, (dp, _) in enumerate(pts) if dp <= dp_budget_pa]
    if not in_budget:
        raise ValueError("No candidate meets the pressure-drop budget.")
    if r_th_target is not None:
        meets = [i for i in in_budget if pts[i][1] <= r_th_target]
        if meets:
            return {"index": min(meets, key=lambda i: (pts[i][0], pts[i][1])), "rule": "meets_target_lowest_dp"}
    return {"index": min(in_budget, key=lambda i: (pts[i][1], pts[i][0])), "rule": "fallback_lowest_r_th"}


def marginal_return(dp_i: float, r_i: float, dp_j: float, r_j: float) -> float:
    """Return of moving from i to a higher-ΔP point j: e = -ln(R_j / R_i) / ln(ΔP_j / ΔP_i) (log-log slope, sign flipped)."""
    return float(-np.log(r_j / r_i) / np.log(dp_j / dp_i))


def recommend_marginal_return(candidates, dp_budget_pa: float, r_th_threshold: float | None,
                              r_th_objective: float | None = None, cutoff: float = 0.25) -> dict:
    """Recommendation rule (SPEC Amendment 3) over (dp_pa, r_th_k_per_w) pairs, e.g. CFD-verified designs.

    Eligible: dp <= budget and r_th <= threshold. Take their non-dominated front sorted by dp ascending and start at its
    lowest-dp point. From the current point i, compute e(i -> j) for every front point j further along; if any j has
    e >= cutoff, move to the nearest such j (lowest dp) and repeat; otherwise stop and recommend i ("marginal_return").
    If no candidate within budget meets the threshold (or there is no threshold), recommend the lowest r_th within budget
    ("fallback_lowest_r_th", meets_threshold False). Raises ValueError if nothing is within budget.

    Returns {"index", "rule", "meets_threshold", "meets_objective", "front": candidate indices of the eligible front,
    "steps": [{"from": i, "returns": [{"to": j, "return": e}, ...], "to": j or None}, ...]}; all indices are positions in
    candidates. The last step has "to": None (no further point reaches the cut-off).
    """
    pts = [(float(dp), float(r)) for dp, r in candidates]
    in_budget = [i for i, (dp, _) in enumerate(pts) if dp <= dp_budget_pa]
    if not in_budget:
        raise ValueError("No candidate meets the pressure-drop budget.")
    eligible = [i for i in in_budget if r_th_threshold is not None and pts[i][1] <= r_th_threshold]

    def result(index, rule, front, steps):
        r = pts[index][1]
        return {"index": index, "rule": rule, "meets_threshold": r_th_threshold is not None and r <= r_th_threshold,
                "meets_objective": r_th_objective is not None and r <= r_th_objective, "front": front, "steps": steps}

    if not eligible:
        best = min(in_budget, key=lambda i: (pts[i][1], pts[i][0]))
        return result(best, "fallback_lowest_r_th", [], [])
    front = [eligible[k] for k in pareto_indices([pts[i][0] for i in eligible], [pts[i][1] for i in eligible])]
    pos, steps = 0, []
    while True:
        i = front[pos]
        returns = [{"to": j, "return": marginal_return(*pts[i], *pts[j])} for j in front[pos + 1:]]
        nxt = next((k for k, e in enumerate(returns) if e["return"] >= cutoff), None)
        steps.append({"from": i, "returns": returns, "to": None if nxt is None else returns[nxt]["to"]})
        if nxt is None:
            return result(i, "marginal_return", front, steps)
        pos += nxt + 1


def recommend_marginal_return_objective(candidates, dp_budget_pa: float, r_th_threshold: float | None,
                                        r_th_objective: float | None = None, cutoff: float = 0.25,
                                        min_dp_ratio: float | None = 1.10, objective_clause: bool = True) -> dict:
    """Recommendation rule of SPEC Amendment 4: the Amendment 3 marginal-return walk with a minimum-ΔP-ratio guard and an
    objective clause, over (dp_pa, r_th_k_per_w) pairs.

    Guard: from the current point i, a further front point j with ΔP_j / ΔP_i < min_dp_ratio is ignored (its return is still
    computed and listed, flagged "skipped"), so two nearly coincident designs cannot produce a large return. None or 1.0 turns
    the guard off.
    Clause: if the walk stops at a point that misses the objective and a further point on the eligible front meets it, continue
    to the lowest-ΔP front point that meets it ("marginal_return_then_objective"). If none does, stay at the stop point
    ("marginal_return", objective_unreachable True).
    With min_dp_ratio None and objective_clause False this is recommend_marginal_return.

    Returns the keys of recommend_marginal_return (each listed return also carries "skipped") plus "walk_stop" (where the
    walk stopped), "objective_step" ({"from", "to"} or None), "objective_unreachable", "cutoff" and "min_dp_ratio".
    """
    pts = [(float(dp), float(r)) for dp, r in candidates]
    in_budget = [i for i, (dp, _) in enumerate(pts) if dp <= dp_budget_pa]
    if not in_budget:
        raise ValueError("No candidate meets the pressure-drop budget.")
    eligible = [i for i in in_budget if r_th_threshold is not None and pts[i][1] <= r_th_threshold]
    meets_obj = lambda i: r_th_objective is not None and pts[i][1] <= r_th_objective

    def result(index, rule, front, steps, walk_stop=None, objective_step=None, unreachable=False):
        r = pts[index][1]
        return {"index": index, "rule": rule, "meets_threshold": r_th_threshold is not None and r <= r_th_threshold,
                "meets_objective": meets_obj(index), "front": front, "steps": steps, "walk_stop": walk_stop,
                "objective_step": objective_step, "objective_unreachable": unreachable, "cutoff": cutoff,
                "min_dp_ratio": min_dp_ratio}

    if not eligible:
        best = min(in_budget, key=lambda i: (pts[i][1], pts[i][0]))
        return result(best, "fallback_lowest_r_th", [], [])
    front = [eligible[k] for k in pareto_indices([pts[i][0] for i in eligible], [pts[i][1] for i in eligible])]
    pos, steps = 0, []
    while True:
        i = front[pos]
        returns = [{"to": j, "return": marginal_return(*pts[i], *pts[j]),
                    "skipped": min_dp_ratio is not None and pts[j][0] / pts[i][0] < min_dp_ratio} for j in front[pos + 1:]]
        nxt = next((k for k, e in enumerate(returns) if not e["skipped"] and e["return"] >= cutoff), None)
        steps.append({"from": i, "returns": returns, "to": None if nxt is None else returns[nxt]["to"]})
        if nxt is None:
            break
        pos += nxt + 1
    stop = front[pos]
    if r_th_objective is None or meets_obj(stop) or not objective_clause:
        return result(stop, "marginal_return", front, steps, stop)
    target = next((j for j in front[pos + 1:] if meets_obj(j)), None)
    if target is None:
        return result(stop, "marginal_return", front, steps, stop, unreachable=True)
    return result(target, "marginal_return_then_objective", front, steps, stop, {"from": stop, "to": target})


CUTOFF_SWEEP = tuple(round(0.05 * k, 2) for k in range(1, 13))  # 0.05 to 0.60 in steps of 0.05


def cutoff_sensitivity(candidates, dp_budget_pa: float, r_th_threshold: float | None, r_th_objective: float | None = None,
                       cutoffs=CUTOFF_SWEEP, min_dp_ratio: float | None = 1.10, objective_clause: bool = True) -> list[dict]:
    """The Amendment 4 rule re-run at each cut-off: [{"cutoff", "index", "rule", "walk_stop", "meets_objective"}, ...]."""
    out = []
    for c in cutoffs:
        sel = recommend_marginal_return_objective(candidates, dp_budget_pa, r_th_threshold, r_th_objective, float(c),
                                                  min_dp_ratio, objective_clause)
        out.append({"cutoff": float(c), "index": sel["index"], "rule": sel["rule"], "walk_stop": sel["walk_stop"],
                    "meets_objective": sel["meets_objective"]})
    return out


def stable_cutoff_range(rows: list[dict], cutoff: float, key: str = "index") -> tuple[float, float] | None:
    """Lowest and highest swept cut-off of the contiguous run, containing `cutoff`, over which rows[...][key] is unchanged.
    None if `cutoff` is not one of the swept values."""
    cs = [r["cutoff"] for r in rows]
    k = next((n for n, c in enumerate(cs) if abs(c - cutoff) < 1e-9), None)
    if k is None:
        return None
    lo = hi = k
    while lo > 0 and rows[lo - 1][key] == rows[k][key]:
        lo -= 1
    while hi < len(rows) - 1 and rows[hi + 1][key] == rows[k][key]:
        hi += 1
    return cs[lo], cs[hi]


def local_returns(x_front, y_front) -> list[dict]:
    """Return e between consecutive points of a front sorted by x: [{"i", "j", "x_i", "x_j", "dp_ratio", "return"}, ...]."""
    return [{"i": k, "j": k + 1, "x_i": float(x_front[k]), "x_j": float(x_front[k + 1]),
             "dp_ratio": float(x_front[k + 1] / x_front[k]),
             "return": marginal_return(x_front[k], y_front[k], x_front[k + 1], y_front[k + 1])} for k in range(len(x_front) - 1)]


def capacity_ratios(dp_pa: float, r_th_k_per_w: float, heat_load_w: float, flow_m3_s: float,
                    r_th_threshold: float) -> dict:
    """Ratios reported for each CFD-verified design (SPEC Amendment 3).

    capacity_kw: heat load the design could carry at the threshold rise, (threshold × heat load) / R_th / 1000, i.e. the rise
    threshold × heat load (20 K for 20 mK/W at 1 kW) divided by R_th; capacity_kw_per_kpa = capacity_kw / (ΔP / 1000);
    kw_per_pumping_w = heat_load / (ΔP × volumetric flow) / 1000, heat load in kW per watt of hydraulic pumping power.
    """
    capacity_kw = r_th_threshold * heat_load_w / r_th_k_per_w / 1e3
    return {"capacity_kw": capacity_kw, "capacity_kw_per_kpa": capacity_kw / (dp_pa / 1e3),
            "kw_per_pumping_w": heat_load_w / (dp_pa * flow_m3_s) / 1e3}
