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
