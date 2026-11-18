"""Published text (README.md, docs/index.html, the heated-width row of ASSUMPTIONS.md) against the committed results.

Every expected number is derived from the committed results files (or, for the screen counts, from re-running the analytic
1D screen on the committed requirements file, as `run_study.py --plot-only` does). Nothing the study produces is hard-coded
here except the short list of superseded values that must not reappear.

How a check works: a claim is a regex for one published sentence, table row or cell, with every number in it captured as a
complete numeric token (no digit or decimal point directly before or after). Each claim must occur at least a set number of
times, and in every occurrence each captured number must equal the value from the results rounded to the precision the text
prints (rounded up, for error bounds). Two generic checks cover every "R mK/W at P kPa" pair and every table row that names a
straight-fin geometry. `failures(texts)` takes the text as an argument, so the mutation tests at the end corrupt an in-memory
copy and show that each corruption is caught.
"""
import json
import math
import re
from functools import lru_cache
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "examples" / "output"
TEXT = {"README.md": ROOT / "README.md", "docs/index.html": ROOT / "docs" / "index.html",
        "ASSUMPTIONS.md": ROOT / "ASSUMPTIONS.md"}
PUBLISHED = ("README.md", "docs/index.html")
SUPERSEDED = ("14.93", "672", "258")  # values of earlier studies that must not reappear in the published text

# A complete numeric token (4,485 or 0.39 or 18), never part of a longer number.
N = r"(?<![\d.,])(\d{1,3}(?:,\d{3})+|\d+(?:\.\d+)?)(?![\d])(?!\.\d)"


def _json(path: Path):
    return json.loads(path.read_text())


@lru_cache(maxsize=None)
def ctx() -> dict:
    from run_study import load
    from tools.fluids import fluid_props
    from tools.screen_1d import screen_1d

    l1 = _json(OUT / "gpu_coldplate" / "results.json")
    pv = _json(OUT / "level2_planview" / "results.json")
    _, req, space, _ = load(ROOT / "examples" / "gpu_coldplate.yaml")
    sweep = screen_1d(req, fluid_props(req.coolant, req.inlet_temp_c), space)
    return {"l1": l1, "pv": pv, "en": _json(OUT / "entrance_effects" / "results.json"),
            "jets": _json(OUT / "level2_slot_jet" / "results.json"),
            "g3": _json(ROOT / "validation" / "level2_planview_straight.json"),
            "duct": _json(ROOT / "validation" / "duct_fully_developed.json"),
            "l2_entrance": _json(ROOT / "validation" / "level2_thermal_entrance.json"),
            "req": req, "n_sweep": len(sweep), "n_feasible": sum(s.feasible for s in sweep),
            "by_id": {d["cfd"]["design_id"]: d for d in l1["designs"]},
            "pv_by_id": {d["design_id"]: d for d in pv["designs"] if d.get("status") == "converged"}}


# ---------------------------------------------------------------- derived quantities
def rec(c):
    return c["by_id"][c["l1"]["recommended"]]


def geom(d):
    g = d["one_d"]["design"]
    return [g["w_ch_m"] * 1e3, g["t_fin_m"] * 1e3, g["h_fin_m"] * 1e3]


def rp(d):
    return [d["cfd"]["r_th_k_per_w"] * 1e3, d["cfd"]["dp_pa"] / 1e3]


def margin(c, r_k_per_w):
    return 100 * (1 - r_k_per_w / c["l1"]["r_th_objective_k_per_w"])


def objective_mk(c):
    return c["l1"]["r_th_objective_k_per_w"] * 1e3


def flow(c, q):
    return next(x for x in c["l1"]["flow_sensitivity"] if abs(x["flow_lpm"] - q) < 1e-9)


def n_cfd(c):
    return len(c["l1"]["designs"])


def wetted_ratios(c):
    """Heated wetted area (channel floor plus both fin sides) over the modelled heated footprint (n x pitch x length)."""
    out = {}
    for i, d in c["by_id"].items():
        g = d["one_d"]["design"]
        out[i] = (g["w_ch_m"] + 2 * g["h_fin_m"]) / (g["w_ch_m"] + g["t_fin_m"])
    return out


def heated_width_mm(c):
    g, n = rec(c)["one_d"]["design"], rec(c)["one_d"]["n_channels"]
    return n * (g["w_ch_m"] + g["t_fin_m"]) * 1e3


def entrance_rec(c):
    en = c["en"]
    return next(r for r in en["one_d"] if r["design_id"] == en["level1_recommended"])["cfd_anchored"]


def entrance_resolved(c, k):
    """k = 0: the widest channel (first resolved design); 'rec': the recommended design."""
    res = c["en"]["resolved"]
    return next(r for r in res if r["design_id"] == c["en"]["level1_recommended"]) if k == "rec" else res[k]


def step_return(steps, frm, to):
    st = next(s for s in steps if s["from"] == frm)
    return next(e["return"] for e in st["returns"] if e["to"] == to)


def walk_return(c, design_id):
    """Return shown for a walk-table row: to the design moved to, or the best return to any further design at the stop."""
    st = next(s for s in c["l1"]["selection_steps"] if s["from"] == design_id)
    if st["to"] is None:
        return max(e["return"] for e in st["returns"])
    return next(e["return"] for e in st["returns"] if e["to"] == st["to"])


def last_step_returns(c):
    """Return of the last step into the pick, fully developed and on entrance-corrected values."""
    pick = c["l1"]["recommended"]
    frm = next(st["from"] for st in c["l1"]["selection_steps"] if st["to"] == pick)
    return [step_return(c["l1"]["selection_steps"], frm, pick), step_return(c["en"]["selection_core"]["steps"], frm, pick)]


def cutoff_groups(c):
    rows, (lo, hi) = c["l1"]["cutoff_sensitivity"], c["l1"]["cutoff_stable_range"]
    above = [r for r in rows if r["cutoff"] > hi + 1e-9]
    assert len({r["recommended"] for r in above}) == 1, "cut-offs above the stable range no longer pick one design"
    below, seen = [], []
    for r in sorted((r for r in rows if r["cutoff"] < lo - 1e-9), key=lambda r: -r["cutoff"]):
        if r["recommended"] not in seen:
            seen.append(r["recommended"])
            below.append(r)
    return lo, hi, above, below


def pv_sel(c, topology=None):
    s = c["pv"]["summary"]
    return s["recommended"] if topology is None else s["per_topology"][topology]["recommended"]


def pv_params(c, topology=None):
    return c["pv_by_id"][pv_sel(c, topology)["design_id"]]["params_mm"]


def pv_rp(c, design_id):
    r = c["pv_by_id"][design_id]["result"]
    return [r["r_th_k_per_w"] * 1e3, r["dp_pa"] / 1e3]


def pv_pick_rp(c, topology=None):
    return pv_rp(c, pv_sel(c, topology)["design_id"])


def vs_straight(c, *topologies):
    rows = [v for v in c["pv"]["summary"]["vs_straight_front"] if v["topology"] in topologies]
    diffs = [abs(v["difference_pct"]) for v in rows]
    return [min(diffs), max(diffs), min(v["dp_pa"] for v in rows) / 1e3, max(v["dp_pa"] for v in rows) / 1e3]


def best_jet(c):
    b = min(c["jets"], key=lambda x: x["result"]["r_th_k_per_w"])["result"]
    return [b["r_th_k_per_w"] * 1e3, b["dp_pa"] / 1e3]


def duct_max_pct(c):
    finest = {}
    for r in c["duct"]:
        if r["alpha"] not in finest or r["h"] < finest[r["alpha"]]["h"]:
            finest[r["alpha"]] = r
    return max(max(abs(r["f_re_err_pct"]), abs(r["nu_err_pct"])) for r in finest.values())


def ns_energy_max_pct(c):
    return 100 * max(abs(c["pv"]["summary"]["max_energy_balance_error"]), abs(c["l2_entrance"]["energy_balance_error"]))


def bias(c):
    """Plan view below the Level 1 CFD on the recommended straight design, in % (R_th, ΔP)."""
    return [-100 * c["g3"]["g3_r_th_k_per_w"]["difference"], -100 * c["g3"]["g3_dp_pa"]["difference"]]


def rise_k(c):
    return c["req"].r_th_threshold_k_per_w * c["req"].heat_load_w


def up(v):
    """An error bound: compared after rounding up to the printed precision."""
    return ("ceil", v)


def _by_geometry(c, g):
    """Level 1 CFD design whose geometry equals the captured (w, t, h) strings at their printed precision."""
    hits = [d for d in c["by_id"].values() if all(agrees(s, v) for s, v in zip(g, geom(d)))]
    return hits[0] if len(hits) == 1 else None


def _flow_row(c, g):
    f = flow(c, float(g[0]))
    return [f["flow_lpm"], f["dp_pa"] / 1e3, f["r_th_k_per_w"] * 1e3, f["base_avg_temp_c"], f["base_max_temp_c"]]


def _walk_row(c, g):
    d = _by_geometry(c, g[:3])
    return None if d is None else geom(d) + rp(d) + [walk_return(c, d["cfd"]["design_id"])]


def _reference_row(key):
    return lambda c, g: geom(c["by_id"][c["l1"]["references"][key]])


def _cutoff_readme(c, g):
    lo, hi, above, _ = cutoff_groups(c)
    d = c["by_id"][above[0]["recommended"]]
    return [lo, hi, above[0]["cutoff"], geom(d)[0], rp(d)[1]]


def _cutoff_index(c, g):
    lo, hi, above, below = cutoff_groups(c)
    d = c["by_id"][above[0]["recommended"]]
    return ([rp(rec(c))[1], lo, hi, rp(d)[1], geom(d)[0], rp(d)[0], above[0]["cutoff"], above[-1]["cutoff"], lo]
            + [r["dp_pa"] / 1e3 for r in below])


def _walk_start(c, g):
    s = c["l1"]["selection_steps"]
    first, second = c["by_id"][s[0]["from"]], c["by_id"][s[1]["from"]]
    return [geom(first)[0]] + rp(first) + [geom(second)[0], step_return(s, s[0]["from"], s[0]["to"]),
                                           step_return(s, s[1]["from"], s[1]["to"])]


def _offset_first(c):
    return pv_rp(c, pv_sel(c, "offset_strip")["walk_stop"])


def _entrance_values(c):
    e = entrance_rec(c)
    return [e["r_th_fd"] * 1e3, e["r_th_combined"] * 1e3, e["dp_fd"] / 1e3, e["dp_combined_core"] / 1e3]


def _resolved_rec(c, *keys):
    return [abs(100 * entrance_resolved(c, "rec")["comparison"][k]) for k in keys]


# ---------------------------------------------------------------- claims
R, I, A, BOTH = ("README.md",), ("docs/index.html",), ("ASSUMPTIONS.md",), PUBLISHED

# (files, label, pattern with {N} for each number, builder(c, captured) -> expected values, minimum occurrences per file)
CLAIMS = [
    # Level 1 recommendation
    (R, "recommended design, summary table", r"Straight copper fins: {N} mm channels, {N} mm fins, {N} mm tall\. {N} mK/W at {N} kPa",
     lambda c, g: geom(rec(c)) + rp(rec(c)), 1),
    (R, "recommended design, results table",
     r"\| Recommended \(marginal return, meets the {N} mK/W objective\) \| {N} / {N} / {N} mm \| {N} mK/W \| {N} kPa \|",
     lambda c, g: [objective_mk(c)] + geom(rec(c)) + rp(rec(c)), 1),
    (I, "recommended design", r"Straight copper fins with {N} mm channels, {N} mm fins, {N} mm tall\. The CFD gives {N} mK/W at {N} kPa, "
     r"which puts the plate base at {N} °C on average and {N} °C at the outlet end",
     lambda c, g: geom(rec(c)) + rp(rec(c)) + [rec(c)["cfd"]["base_avg_temp_c"], rec(c)["cfd"]["base_max_temp_c"]], 1),
    (I, "recommended design, figure alt text", r"circled at {N} mK/W and {N} kPa", lambda c, g: rp(rec(c)), 1),
    (I, "objective margin", r"with a {N}% margin", lambda c, g: [margin(c, rec(c)["cfd"]["r_th_k_per_w"])], 2),
    (I, "flow at which the objective is missed", r"misses it at {N} LPM",
     lambda c, g: [next(f["flow_lpm"] for f in c["l1"]["flow_sensitivity"]
                        if f["r_th_k_per_w"] > c["l1"]["r_th_objective_k_per_w"])], 1),
    (I, "1D split of the recommended design", r"The 1D model gives {N} mK/W for this design and splits it into about {N} for fins and "
     r"convection, {N} for the coolant warming up as it crosses the plate, and {N} for conduction through the base",
     lambda c, g: [1e3 * rec(c)["one_d"][k] for k in ("r_th_k_per_w", "r_conv", "r_caloric", "r_base")], 1),
    (R, "capacity of the recommended design", r"carry at a {N} K rise \({N} kW for the recommended one\)",
     lambda c, g: [rise_k(c), rec(c)["ratios"]["capacity_kw"]], 1),
    (I, "allowable rise for the capacity ratio", r"Heat-removal capacity at the {N} K allowable base-to-inlet rise",
     lambda c, g: [rise_k(c)], 1),
    # flow sensitivity
    (R, "flow sensitivity", r"averages {N} °C at the chip interface; at {N} LPM it rises to {N} mK/W \({N} °C\), at {N} LPM it falls to "
     r"{N} mK/W \({N} °C\)",
     lambda c, g: [flow(c, 1.5)["base_avg_temp_c"], 1.0, flow(c, 1.0)["r_th_k_per_w"] * 1e3, flow(c, 1.0)["base_avg_temp_c"],
                   2.0, flow(c, 2.0)["r_th_k_per_w"] * 1e3, flow(c, 2.0)["base_avg_temp_c"]], 1),
    (I, "flow sensitivity table row", r'<td class="n">{N} LPM</td><td class="n">{N} kPa</td><td class="n">{N} mK/W</td>'
     r'<td class="n">{N} °C</td><td class="n">{N} °C</td>', _flow_row, 3),
    # counts
    (R, "designs screened and solved", r"screens {N} fin geometries, a conjugate CFD model checks all {N} on the Pareto front",
     lambda c, g: [c["n_sweep"], n_cfd(c)], 1),
    (R, "designs screened, flowchart", r"1D screen<br/>{N} designs", lambda c, g: [c["n_sweep"]], 1),
    (R, "designs solved, method", r"within the budget \({N} here\)", lambda c, g: [n_cfd(c)], 1),
    (I, "designs screened and feasible", r"It covers {N} designs in under 20 milliseconds, and {N} of them",
     lambda c, g: [c["n_sweep"], c["n_feasible"]], 1),
    (I, "designs solved", r"All {N} designs on the 1D front", lambda c, g: [n_cfd(c)], 1),
    (I, "designs solved, parity caption", r"1D against CFD for the {N} designs", lambda c, g: [n_cfd(c)], 1),
    (I, "designs solved, entrance", r"Across the {N} designs the order", lambda c, g: [n_cfd(c)], 1),
    (I, "designs solved, wetted area", r"For the {N} straight-fin designs checked with CFD", lambda c, g: [n_cfd(c)], 1),
    # parity, every occurrence in both files
    (BOTH, "1D vs CFD difference", r"{N}% on average,? (?:and )?{N}% at worst",
     lambda c, g: [c["l1"]["parity"]["mean_abs_pct"], c["l1"]["parity"]["max_abs_pct"]], 1),
    (R, "1D vs CFD, design count", r"{N}% at worst \({N} designs\)", lambda c, g: [c["l1"]["parity"]["max_abs_pct"], n_cfd(c)], 1),
    # reference designs (their values are checked by the geometry-row check; here each label is tied to its design)
    (R, "reference: lowest ΔP meeting the objective", r"\| Lowest ΔP meeting the 15 mK/W objective \| {N} / {N} / {N} mm \|",
     _reference_row("lowest_dp_meeting_objective"), 1),
    (R, "reference: lowest ΔP meeting the threshold", r"\| Lowest ΔP meeting the 20 mK/W threshold \| {N} / {N} / {N} mm \|",
     _reference_row("lowest_dp_meeting_threshold"), 1),
    (R, "reference: lowest resistance within budget", r"\| Lowest resistance within 5 kPa \| {N} / {N} / {N} mm \|",
     _reference_row("lowest_r_th_within_budget"), 1),
    # selection walk and cut-off sensitivity
    (I, "walk table row", r'<tr><td>{N} / {N} / {N} mm</td><td class="n">{N} mK/W</td><td class="n">{N} kPa</td><td class="n">{N}</td>',
     _walk_row, 3),
    (I, "start of the walk", r"the walk starts at {N} mm channels \({N} mK/W at {N} kPa\), the lowest pressure drop that meets the 20 mK/W "
     r"limit, and moves through {N} mm channels with returns of {N} and {N}", _walk_start, 1),
    (I, "best return at the stop", r"In the last row, {N} is the best return", lambda c, g: [walk_return(c, c["l1"]["recommended"])], 1),
    (I, "cut-off in use", r"each 1% of extra pressure drop buys at least {N}% lower resistance\. The {N} cut-off",
     lambda c, g: [c["l1"]["marginal_return_cutoff"]] * 2, 1),
    (R, "cut-off sensitivity", r"the pick holds from {N} to {N}; from {N} up it is the {N} mm design at {N} kPa", _cutoff_readme, 1),
    (I, "cut-off sensitivity", r"the pick stays at {N} kPa for cut-offs of {N} and {N} and at {N} kPa \({N} mm channels, {N} mK/W\) from "
     r"{N} to {N}\. Below {N} it moves to {N} kPa and then {N} kPa", _cutoff_index, 1),
    # entrance effects
    (R, "entrance-corrected design, summary table",
     r"entrance effects are quantified separately \({N} mK/W at {N} kPa, objective margin {N}%",
     lambda c, g: [_entrance_values(c)[1], _entrance_values(c)[3], margin(c, entrance_rec(c)["r_th_combined"])], 1),
    (R, "entrance-corrected design", r"lowers R_th from {N} to {N} mK/W and raises the channel-core ΔP from {N} to {N} kPa, so the "
     r"objective margin grows from {N}% to {N}%",
     lambda c, g: _entrance_values(c) + [margin(c, entrance_rec(c)["r_th_fd"]), margin(c, entrance_rec(c)["r_th_combined"])], 1),
    (I, "entrance-corrected design", r"the correction lowers the resistance from {N} to {N} mK/W and raises the pressure drop from {N} to "
     r"{N} kPa, so the margin on the objective grows to {N}%",
     lambda c, g: _entrance_values(c) + [margin(c, entrance_rec(c)["r_th_combined"])], 1),
    (I, "resolved entrance increments", r"raises the length-mean heat transfer coefficient by {N}% and the channel pressure drop by {N}%",
     lambda c, g: _resolved_rec(c, "h_increase_resolved", "dp_increase_resolved"), 1),
    (I, "entrance agreement, recommended", r"matches the resulting changes in resistance and pressure drop to within {N}% and {N}%",
     lambda c, g: _resolved_rec(c, "r_th_contribution_disagreement", "dp_contribution_disagreement"), 1),
    (I, "entrance agreement, widest channel", r"the widest channel misses by {N}% on resistance",
     lambda c, g: [abs(100 * entrance_resolved(c, 0)["comparison"]["r_th_contribution_disagreement"])], 1),
    (I, "entrance, earlier formulation", r"An earlier version of the correction missed them by {N} to {N}%",
     lambda c, g: sorted(_resolved_rec(c, "earlier_r_th_contribution_disagreement", "earlier_dp_contribution_disagreement")), 1),
    (I, "entrance, last-step return", r"the return on the last step of the selection falls from {N} to {N}",
     lambda c, g: last_step_returns(c), 1),
    # straight fins against a flat base; heated width
    (I, "wetted area ratio", r"heated wetted area \(channel floor plus both fin sides\) is {N} to {N} times the heated footprint, and {N} "
     r"for the recommended design",
     lambda c, g: [min(wetted_ratios(c).values()), max(wetted_ratios(c).values()), wetted_ratios(c)[c["l1"]["recommended"]]], 1),
    (A, "heated width", r"slightly under the {N} mm footprint: {N} mm \({N} × {N} mm\)",
     lambda c, g: [c["req"].footprint_width_m * 1e3, heated_width_mm(c), rec(c)["one_d"]["n_channels"],
                   heated_width_mm(c) / rec(c)["one_d"]["n_channels"]], 1),
    # slot jets
    (R, "best slot-jet array", r"the best reaches {N} mK/W at {N} kPa", lambda c, g: best_jet(c), 1),
    (I, "best slot-jet array", r"the best of {N} arrays reached {N} mK/W", lambda c, g: [len(c["jets"]), best_jet(c)[0]], 1),
    # plan view
    (R, "plan-view designs solved", r"\({N} plan-view designs solved", lambda c, g: [len(c["pv_by_id"])], 1),
    (R, "plan-view designs excluded", r"The {N} designs whose fins are taller", lambda c, g: [len(c["pv"]["summary"]["excluded"])], 1),
    (R, "cross-topology pick", r"pins of {N} mm diameter at {N} × {N} mm pitch \({N} mK/W at {N} kPa\) across all topologies",
     lambda c, g: [pv_params(c)[k] for k in ("d", "s_t", "s_l")] + pv_pick_rp(c), 1),
    (R, "straight-fin plan-view pick", r"and {N} mm straight fins and channels \({N} mK/W at {N} kPa\) among straight fins",
     lambda c, g: [pv_params(c, "straight")["w"]] + pv_pick_rp(c, "straight"), 1),
    (R, "offset-strip walk stop", r"walk stops at its first point \({N} mK/W\)", lambda c, g: _offset_first(c)[:1], 1),
    (R, "offset-strip pick", r"meets {N} mK/W \({N} mK/W at {N} kPa\)",
     lambda c, g: [objective_mk(c)] + pv_pick_rp(c, "offset_strip"), 1),
    (I, "offset-strip walk stop", r"walk stops at its first design \({N} mK/W at {N} kPa\), which misses the {N} mK/W objective",
     lambda c, g: _offset_first(c) + [objective_mk(c)], 1),
    (I, "plan-view table, circular pins", r'<tr><td>Circular pins</td><td>{N} mm pins, {N} × {N} mm pitch</td><td class="n">{N} mK/W</td>'
     r'<td class="n">{N} kPa</td></tr>',
     lambda c, g: [pv_params(c, "pin_circle")[k] for k in ("d", "s_t", "s_l")] + pv_pick_rp(c, "pin_circle"), 1),
    (I, "plan-view table, straight",
     r'<tr><td>Straight</td><td>{N} mm channels and fins</td><td class="n">{N} mK/W</td><td class="n">{N} kPa</td></tr>',
     lambda c, g: [pv_params(c, "straight")["w"]] + pv_pick_rp(c, "straight"), 1),
    (I, "plan-view table, zigzag",
     r'<tr><td>Zigzag</td><td>{N} mm gap, {N} mm amplitude</td><td class="n">{N} mK/W</td><td class="n">{N} kPa</td></tr>',
     lambda c, g: [pv_params(c, "zigzag")[k] for k in ("w", "amp")] + pv_pick_rp(c, "zigzag"), 1),
    (I, "plan-view table, sinusoidal",
     r'<tr><td>Sinusoidal</td><td>{N} mm gap, {N} mm amplitude</td><td class="n">{N} mK/W</td><td class="n">{N} kPa</td></tr>',
     lambda c, g: [pv_params(c, "sinusoidal")[k] for k in ("w", "amp")] + pv_pick_rp(c, "sinusoidal"), 1),
    (I, "plan-view table, diamond pins",
     r'<tr><td>Diamond pins</td><td>{N} mm side</td><td class="n">{N} mK/W</td><td class="n">{N} kPa</td></tr>',
     lambda c, g: [pv_params(c, "pin_diamond")["s"]] + pv_pick_rp(c, "pin_diamond"), 1),
    (I, "plan-view table, offset strip", r'<tr><td>Offset strip</td><td>{N} mm strips, {N} mm thick, {N} mm tall</td>'
     r'<td class="n">{N} mK/W</td><td class="n">{N} kPa</td></tr>',
     lambda c, g: [pv_params(c, "offset_strip")[k] for k in ("l", "t", "h_fin")] + pv_pick_rp(c, "offset_strip"), 1),
    (I, "against straight fins, circular pins", r"Circular pins: {N} to {N}% lower resistance, between {N} and {N} kPa",
     lambda c, g: vs_straight(c, "pin_circle"), 1),
    (I, "against straight fins, wavy", r"Zigzag and sinusoidal: {N} to {N}% higher", lambda c, g: vs_straight(c, "zigzag", "sinusoidal")[:2], 1),
    (I, "against straight fins, diamond pins", r"Diamond pins: {N} to {N}% higher", lambda c, g: vs_straight(c, "pin_diamond")[:2], 1),
    (I, "against straight fins, offset strips", r"Offset strips: {N} to {N}% higher", lambda c, g: vs_straight(c, "offset_strip")[:2], 1),
    (I, "pin lead", r"A {N} to {N}% lead for pins", lambda c, g: vs_straight(c, "pin_circle")[:2], 1),
    (R, "pin lead", r"only by {N}–{N}% at equal pressure drop", lambda c, g: vs_straight(c, "pin_circle")[:2], 1),
    (R, "plan-view tall-fin bias", r"reads ΔP {N}% and R_th {N}% below", lambda c, g: [bias(c)[1], bias(c)[0]], 1),
    (R, "plan-view tall-fin bias, headline", r"inside the {N}% this model reads low", lambda c, g: [bias(c)[0]], 1),
    (I, "plan-view tall-fin bias", r"gives resistance {N}% lower than the cross-section CFD", lambda c, g: [bias(c)[0]], 1),
    # validation; "within" figures are error bounds, rounded up
    (R, "duct toolchain error bound", r"f·Re and Nu to within {N}%", lambda c, g: [up(duct_max_pct(c))], 1),
    (I, "checks table: duct", r'four aspect ratios[^<]*</td><td>[^<]*</td><td>[^<]*</td><td class="n">within {N}%</td>',
     lambda c, g: [up(duct_max_pct(c))], 1),
    (I, "checks table: cross-section mesh",
     r'Cross-section mesh at twice the resolution</td>(?:<td>[^<]*</td>){{2}}<td class="n">{N}% change</td>',
     lambda c, g: [up(c["l1"]["confidence"]["mesh_study"]["r_th_change_pct"])], 1),
    (I, "checks table: developing flow", r'developing flow between plates</td>(?:<td>[^<]*</td>){{2}}<td class="n">within {N}%</td>',
     lambda c, g: [up(100 * abs(c["l2_entrance"]["nu_local_worst_error"]))], 1),
    (I, "checks table: plan-view solver", r'Plan-view solver, straight channel</td>(?:<td>[^<]*</td>){{2}}<td class="n">within {N}%</td>',
     lambda c, g: [up(100 * max(abs(c["g3"]["g1_dpdx_error"]), abs(c["g3"]["g1_nu_error"])))], 1),
    (I, "checks table: plan view against cross-section",
     r'same straight design</td>(?:<td>[^<]*</td>){{2}}<td class="n">{N}% and {N}% low</td>', lambda c, g: bias(c), 1),
    (I, "checks table: energy balance",
     r'Energy balance, every Navier-Stokes run</td>(?:<td>[^<]*</td>){{2}}<td class="n">within {N}%</td>',
     lambda c, g: [up(ns_energy_max_pct(c))], 1),
]


# ---------------------------------------------------------------- comparison and checks
def agrees(s: str, v) -> bool:
    """True when the printed token s equals v at the precision s is printed with (v = ("ceil", x) rounds x up)."""
    mode, v = v if isinstance(v, tuple) else ("round", v)
    if "," in s:
        return mode == "round" and f"{round(v):,}" == s
    d = len(s.split(".")[1]) if "." in s else 0
    if mode == "ceil":
        v = math.ceil(v * 10**d - 1e-9) / 10**d
    return f"{v:.{d}f}" == s


def _claim_failures(flat: dict, c) -> list[str]:
    out = []
    for files, label, pattern, build, minimum in CLAIMS:
        rx = re.compile(pattern.format(N=N))
        for name in files:
            matches = list(rx.finditer(flat[name]))
            if len(matches) < minimum:
                out.append(f"{name}: {label}: expected at least {minimum} occurrence(s), found {len(matches)}")
            for m in matches:
                expected = build(c, m.groups())
                if expected is None:
                    out.append(f"{name}: {label}: '{m.group(0)}' names no design in the results")
                    continue
                for k, (s, v) in enumerate(zip(m.groups(), expected), 1):
                    if not agrees(s, v):
                        how, shown = ("rounded up", v[1]) if isinstance(v, tuple) else ("rounded", v)
                        out.append(f"{name}: {label}: number {k} reads {s} in '{m.group(0)}'; expected {shown:.6g} from the "
                                   f"results, {how} to the printed precision")
    return out


def _known_pairs(c) -> list[tuple[float, float]]:
    """(R_th in mK/W, ΔP in kPa) of every computed design or operating point the text may quote."""
    pairs = [tuple(rp(d)) for d in c["by_id"].values()]
    pairs += [(f["r_th_k_per_w"] * 1e3, f["dp_pa"] / 1e3) for f in c["l1"]["flow_sensitivity"]]
    pairs += [(r["cfd_anchored"]["r_th_combined"] * 1e3, r["cfd_anchored"]["dp_combined_core"] / 1e3) for r in c["en"]["one_d"]]
    pairs += [tuple(pv_rp(c, i)) for i in c["pv_by_id"]]
    pairs += [(j["result"]["r_th_k_per_w"] * 1e3, j["result"]["dp_pa"] / 1e3) for j in c["jets"]]
    return pairs


PAIR_PATTERNS = (r"{N} mK/W (?:at|and) {N} kPa", r'{N} mK/W</td><td class="n">{N} kPa', r"{N} mK/W \| {N} kPa")


def _pair_failures(flat: dict, c) -> list[str]:
    """Every published 'R mK/W at P kPa' pair must be one computed design or operating point, at the printed precision."""
    out, known = [], _known_pairs(c)
    for name in PUBLISHED:
        for p in PAIR_PATTERNS:
            for m in re.finditer(p.format(N=N), flat[name]):
                r, dp = m.groups()
                if not any(agrees(r, kr) and agrees(dp, kp) for kr, kp in known):
                    out.append(f"{name}: '{m.group(0)}' matches no computed design or operating point in the results")
    return out


def _geometry_row_failures(texts: dict, flat: dict, c) -> list[str]:
    """Every table row naming a straight-fin geometry 'w / t / h mm' quotes that design's CFD R_th and ΔP."""
    rows = [("README.md", r) for r in re.findall(r"^\|.*\|$", texts["README.md"], flags=re.M)]
    rows += [("docs/index.html", r) for r in re.findall(r"<tr>.*?</tr>", flat["docs/index.html"])]
    out = []
    for name, row in rows:
        g = re.search(r"{N} / {N} / {N} mm".format(N=N), row)
        if not g:
            continue
        d = _by_geometry(c, g.groups())
        if d is None:
            out.append(f"{name}: row '{row}' names a geometry that is not a CFD-verified design")
            continue
        for unit, v in (("mK/W", rp(d)[0]), ("kPa", rp(d)[1])):
            m = re.search(f"{N} {re.escape(unit)}", row[g.end():])
            if m is None or not agrees(m.group(1), v):
                out.append(f"{name}: row '{row}': {unit} should be {v:.6g} for {d['cfd']['design_id']}")
    return out


def _fronts_capacity_falls(c) -> list[str]:
    """docs/index.html says capacity per kPa falls as ΔP rises on every front in this study; check each front."""
    fronts = {"Level 1 CFD": [(d["cfd"]["dp_pa"], d["ratios"]["capacity_kw_per_kpa"]) for d in c["by_id"].values()]}
    s = c["pv"]["summary"]
    for name, sel in [("cross-topology", s["recommended"])] + list((t, p["recommended"]) for t, p in s["per_topology"].items()):
        if sel:
            fronts[f"plan-view {name}"] = [(c["pv_by_id"][i]["result"]["dp_pa"], s["ratios"][i]["capacity_kw_per_kpa"])
                                           for i in sel["front"]]
    out = []
    for name, pts in fronts.items():
        v = [x for _, x in sorted(pts)]
        if any(b >= a for a, b in zip(v, v[1:])):
            out.append(f"docs/index.html: capacity per kPa does not fall along the {name} front, as the text says it does")
    return out


def failures(texts: dict) -> list[str]:
    """Every disagreement between the published text (file name -> raw text) and the committed results."""
    c = ctx()
    flat = {name: " ".join(t.split()) for name, t in texts.items()}
    out = _claim_failures(flat, c) + _pair_failures(flat, c) + _geometry_row_failures(texts, flat, c)
    if "On every front in this study it falls as pressure drop rises" in flat["docs/index.html"]:
        out += _fronts_capacity_falls(c)
    for name in PUBLISHED:
        for value in SUPERSEDED:
            if re.search(rf"(?<![\d.]){re.escape(value)}(?![\d])", flat[name]):
                out.append(f"{name}: superseded value '{value}' from an earlier study still appears")
    bad = [i + 1 for i, line in enumerate(texts["docs/index.html"].splitlines()) if "—" in line or "–" in line]
    if bad:
        out.append(f"docs/index.html: em or en dash on lines {bad}; expected none (write ranges as 'a to b')")
    return out


def load_texts() -> dict:
    return {name: path.read_text() for name, path in TEXT.items()}


def test_published_text_matches_results():
    found = failures(load_texts())
    assert not found, "\n".join(found)


def test_number_tokens_are_complete():
    assert re.findall(N, "114.0 and 14.0, 88.9% and 8.9%, 4,485") == ["114.0", "14.0", "88.9", "8.9", "4,485"]
    assert agrees("14.0", 14.0487) and not agrees("114.0", 14.0487) and agrees("14.05", 14.0487) and not agrees("14.04", 14.0487)
    assert agrees("0.008", ("ceil", 0.00735)) and not agrees("0.007", ("ceil", 0.00735))


# Corruptions of the published text that the checks must catch: (id, file, exact old text, new text).
MUTATIONS = [
    ("README hero R_th gains a leading digit", "README.md", "3.75 mm tall. 14.0 mK/W at 1.60 kPa", "3.75 mm tall. 114.0 mK/W at 1.60 kPa"),
    ("README entrance margin 88.9%", "README.md", "objective margin 8.9%", "objective margin 88.9%"),
    ("README second parity claim disagrees", "README.md", "to 1.0% on average and 2.9% at worst", "to 1.0% on average and 3.9% at worst"),
    ("README entrance dP 1.627 kPa", "README.md", "from 1.600 to 1.617 kPa", "from 1.600 to 1.627 kPa"),
    ("README entrance R_th 13.56 at two decimals", "README.md", "from 14.05 to 13.66 mK/W", "from 14.05 to 13.56 mK/W"),
    ("README pin pitch 0.8 x 1.0 mm", "README.md", "at 0.8 × 0.9 mm pitch", "at 0.8 × 1.0 mm pitch"),
    ("README screened count", "README.md", "screens 4,485 fin geometries", "screens 4,486 fin geometries"),
    ("README CFD count in hero table", "README.md", "2.9% at worst (18 designs)", "2.9% at worst (19 designs)"),
    ("README duct bound below the saved error", "README.md", "Nu to within 0.008%", "Nu to within 0.007%"),
    ("README reference row R_th", "README.md", "| 0.42 / 0.39 / 3.75 mm | 14.9 mK/W | 1.29 kPa |",
     "| 0.42 / 0.39 / 3.75 mm | 13.9 mK/W | 1.29 kPa |"),
    ("README 2.0 LPM flow sensitivity", "README.md", "falls to 12.8 mK/W (52.8 °C)", "falls to 12.3 mK/W (52.8 °C)"),
    ("README plan-view straight pick", "README.md", "(13.4 mK/W at 1.31 kPa)", "(13.4 mK/W at 1.21 kPa)"),
    ("index entrance dP 1.72 kPa", "docs/index.html", "from 1.60 to 1.62 kPa", "from 1.60 to 1.72 kPa"),
    ("index duct bound below the saved error", "docs/index.html", "within 0.008%", "within 0.007%"),
    ("index energy-balance bound below the saved error", "docs/index.html", "within 0.76%", "within 0.75%"),
    ("index mesh change", "docs/index.html", "0.002% change", "0.02% change"),
    ("index wetted-area range", "docs/index.html", "6.2 to 10.3 times", "6.2 to 8.3 times"),
    ("index widest-channel miss", "docs/index.html", "misses by 10.8%", "misses by 1.8%"),
    ("index walk table row", "docs/index.html", '0.46 / 0.44 / 4.0 mm</td><td class="n">15.9 mK/W',
     '0.46 / 0.44 / 4.0 mm</td><td class="n">16.9 mK/W'),
    ("index flow table 1.0 LPM resistance", "docs/index.html", '1.07 kPa</td><td class="n">16.5 mK/W',
     '1.07 kPa</td><td class="n">15.5 mK/W'),
    ("index pin pitch", "docs/index.html", "0.4 mm pins, 0.8 × 0.9 mm pitch", "0.4 mm pins, 0.9 × 0.9 mm pitch"),
    ("index plan-view bias in prose", "docs/index.html", "resistance 5.6% lower than", "resistance 4.6% lower than"),
    ("index 1D split part", "docs/index.html", "about 7.5 for fins", "about 7.8 for fins"),
    ("index pin advantage range", "docs/index.html", "2.6 to 5.8% lower resistance", "2.6 to 6.8% lower resistance"),
    ("index margin in the introduction", "docs/index.html", "with a 6.3% margin, but", "with a 6.8% margin, but"),
    ("index allowable rise for capacity", "docs/index.html", "at the 20 K allowable", "at the 25 K allowable"),
    ("index entrance gate agreement", "docs/index.html", "to within 3.1% and 2.9%", "to within 3.1% and 1.9%"),
    ("ASSUMPTIONS heated width", "ASSUMPTIONS.md", "59.64 mm (77 × 0.775 mm)", "59.84 mm (77 × 0.775 mm)"),
    # second review round (Codex, 2026-11-18): its ten probes, then more of the same classes
    ("probe: second parity maximum 9.9%", "README.md", "to 1.0% on average and 2.9% at worst.", "to 1.0% on average and 9.9% at worst."),
    ("probe: second parity unit K", "README.md", "to 1.0% on average and 2.9% at worst.", "to 1.0% on average and 2.9 K at worst."),
    ("probe: delete second parity sentence", "README.md",
     "The 1D screen and the CFD agree on thermal resistance to 1.0% on average and 2.9% at worst.\n", ""),
    ("probe: duct bound unit K", "README.md", "Nu to within 0.008%", "Nu to within 0.008 K"),
    ("probe: swap straight and zigzag pressure drops", "docs/index.html",
     '<td class="n">13.4 mK/W</td><td class="n">1.31 kPa</td></tr>\n      <tr><td>Zigzag</td><td>0.4 mm gap, 0.3 mm amplitude</td><td class="n">13.4 mK/W</td><td class="n">1.46 kPa</td>',
     '<td class="n">13.4 mK/W</td><td class="n">1.46 kPa</td></tr>\n      <tr><td>Zigzag</td><td>0.4 mm gap, 0.3 mm amplitude</td><td class="n">13.4 mK/W</td><td class="n">1.31 kPa</td>'),
    ("probe: walk return 0.33 to 0.93", "docs/index.html", '1.05 kPa</td><td class="n">0.33</td>', '1.05 kPa</td><td class="n">0.93</td>'),
    ("probe: pin caption pitch", "README.md", "0.4 mm diameter, 0.7 × 0.6 mm pitch", "0.4 mm diameter, 0.8 × 0.6 mm pitch"),
    ("probe: solver tolerance in README", "README.md", "5% on a limiting Nusselt number", "2% on a limiting Nusselt number"),
    ("probe: reverse capacity trend", "docs/index.html", "On every front in this study it falls as pressure drop rises",
     "On every front in this study it rises as pressure drop rises"),
    ("probe: parity average −1.0%", "README.md", "to 1.0% on average and 2.9% at worst.", "to −1.0% on average and 2.9% at worst."),
    ("parity average -1.0% with a hyphen", "docs/index.html", "to 1.0% on average and 2.9% at worst", "to -1.0% on average and 2.9% at worst"),
    ("solver tolerance in the case study", "docs/index.html", "(5% on a limiting Nusselt number)", "(1% on a limiting Nusselt number)"),
    ("f·Re tolerance in README", "README.md", "(0.2% on f·Re,", "(0.1% on f·Re,"),
    ("index pin caption diameter", "docs/index.html", "Circular pins, 0.4 mm across.", "Circular pins, 0.5 mm across."),
    ("index objective margin unit changed", "docs/index.html", "with a 6.3% margin at rated flow", "with a 6.3 K margin at rated flow"),
    ("index walk table row deleted", "docs/index.html",
     '      <tr><td>0.46 / 0.44 / 4.0 mm</td><td class="n">15.9 mK/W</td><td class="n">1.05 kPa</td><td class="n">0.33</td><td>Keep going</td></tr>\n', ""),
    ("index walk table action swapped", "docs/index.html", '<td class="n">0.19</td><td>Stop here</td>', '<td class="n">0.19</td><td>Keep going</td>'),
    ("README reference table row deleted", "README.md", "| Lowest resistance within 5 kPa | 0.27 / 0.26 / 2.5 mm | 11.9 mK/W | 4.66 kPa |\n", ""),
    ("README new unchecked number", "README.md", "Straight fins remain a sound choice;", "Straight fins remain a sound choice at 17.3 mK/W;"),
    ("index base thickness in caption", "docs/index.html", "The bottom 2 mm is the copper base", "The bottom 3 mm is the copper base"),
    ("index flow table row deleted", "docs/index.html",
     '      <tr><td class="n">2.0 LPM</td><td class="n">2.13 kPa</td><td class="n">12.8 mK/W</td><td class="n">52.8 °C</td><td class="n">56.5 °C</td></tr>\n', ""),
]


@pytest.mark.parametrize("name, old, new", [m[1:] for m in MUTATIONS], ids=[m[0] for m in MUTATIONS])
def test_checks_catch_mutation(name, old, new):
    texts = load_texts()
    assert texts[name].count(old) == 1, f"mutation target not found exactly once in {name}: {old!r}"
    texts[name] = texts[name].replace(old, new)
    assert failures(texts), f"{name}: corrupting {old!r} to {new!r} was not caught"


PRESENTATION_ONLY = [
    ("extra CSS class on walk-table cells", "docs/index.html", '<td class="n">15.9 mK/W</td><td class="n">1.05 kPa</td><td class="n">0.33</td>',
     '<td class="n num">15.9 mK/W</td><td class="n num" style="color: red">1.05 kPa</td><td class="n num">0.33</td>'),
    ("extra CSS class on a plan-view row", "docs/index.html", "<tr><td>Straight</td><td>0.4 mm channels and fins</td>",
     '<tr class="hl"><td class="name">Straight</td><td>0.4 mm channels and fins</td>'),
    ("README table cells padded", "README.md", "| 0.42 / 0.39 / 3.75 mm | 14.9 mK/W | 1.29 kPa |", "|  0.42 / 0.39 / 3.75 mm  |  14.9 mK/W |   1.29 kPa  |"),
]


@pytest.mark.parametrize("name, old, new", [m[1:] for m in PRESENTATION_ONLY], ids=[m[0] for m in PRESENTATION_ONLY])
def test_presentation_changes_do_not_fail(name, old, new):
    texts = load_texts()
    assert texts[name].count(old) == 1, f"target not found exactly once in {name}: {old!r}"
    texts[name] = texts[name].replace(old, new)
    found = failures(texts)
    assert not found, "\n".join(found)


def test_index_images_exist():
    index = TEXT["docs/index.html"]
    refs = re.findall(r'src="(img/[^"]+)"', index.read_text())
    assert refs, "docs/index.html: no img/ references found"
    missing = [r for r in refs if not (index.parent / r).exists()]
    assert not missing, f"docs/index.html: referenced images missing from docs/: {missing}"


def test_readme_images_exist():
    refs = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", TEXT["README.md"].read_text())
    missing = [r for r in refs if not (ROOT / r).exists()]
    assert not missing, f"README.md: referenced images missing: {missing}"
