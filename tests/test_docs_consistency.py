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
from html.parser import HTMLParser
from functools import lru_cache
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "examples" / "output"
TEXT = {"README.md": ROOT / "README.md", "docs/index.html": ROOT / "docs" / "index.html",
        "ASSUMPTIONS.md": ROOT / "ASSUMPTIONS.md"}
PUBLISHED = ("README.md", "docs/index.html")
SUPERSEDED = ("14.93", "672", "258")  # values of earlier studies that must not reappear in the published text

# A complete numeric token (4,485 or 0.39 or 18 or -1.0 / \u22121.0 with its sign), never part of a longer number or a
# name such as 1D or PG25.
N = r"(?<![\w.,\-\u2212])([-\u2212]?(?:\d{1,3}(?:,\d{3})+|\d+(?:\.\d+)?))(?![\d])(?!\.\d)(?![A-Za-z])"


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
            "fields": _json(OUT / "level2_planview" / "fields.json"),
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


def _walk_row(c, g):
    d = _by_geometry(c, g[:3])
    if d is None:
        return None
    stop = d["cfd"]["design_id"] == c["l1"]["walk_stop"]
    return geom(d) + rp(d) + [walk_return(c, d["cfd"]["design_id"]), "Stop here" if stop else "Keep going"]


def _reference(key):
    limit = {"lowest_dp_meeting_objective": objective_mk,
             "lowest_dp_meeting_threshold": lambda c: c["l1"]["r_th_threshold_k_per_w"] * 1e3,
             "lowest_r_th_within_budget": lambda c: c["l1"]["dp_budget_pa"] / 1e3}[key]

    def build(c, g):
        d = c["by_id"][c["l1"]["references"][key]]
        return [limit(c)] + geom(d) + rp(d)
    return build


def _walk_summary(c, g):
    steps = [st for st in c["l1"]["selection_steps"] if st["to"] is not None]
    last = c["by_id"][steps[-1]["from"]]
    r0, r1 = rp(last), rp(rec(c))
    return ([geom(c["by_id"][steps[0]["from"]])[0]] + [walk_return(c, st["from"]) for st in steps]
            + [geom(last)[0], geom(rec(c))[0], 100 * (1 - r1[0] / r0[0]), r1[1] / r0[1], walk_return(c, c["l1"]["recommended"])])


def _pin_field(c):
    return c["fields"]["pin_circle"]["params_mm"]


def production_tolerances() -> dict:
    """rel= tolerances of the production-solver tests in tests/test_run_cfd.py, keyed by the asserted attribute."""
    import ast
    tree = ast.parse((ROOT / "tests" / "test_run_cfd.py").read_text())
    out = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Compare) and isinstance(node.left, ast.Attribute) and len(node.comparators) == 1:
            call = node.comparators[0]
            if isinstance(call, ast.Call) and getattr(call.func, "attr", None) == "approx":
                for kw in call.keywords:
                    if kw.arg == "rel":
                        out[node.left.attr] = ast.literal_eval(kw.value)
    return out


def capacity_direction(c) -> str:
    """'falls' if capacity per kPa falls strictly with ΔP along every front in the study, 'rises' if it rises, else 'mixed'."""
    fronts = [[(d["cfd"]["dp_pa"], d["ratios"]["capacity_kw_per_kpa"]) for d in c["by_id"].values()]]
    s = c["pv"]["summary"]
    for sel in [s["recommended"]] + [p["recommended"] for p in s["per_topology"].values()]:
        if sel:
            fronts.append([(c["pv_by_id"][i]["result"]["dp_pa"], s["ratios"][i]["capacity_kw_per_kpa"]) for i in sel["front"]])
    series = [[x for _, x in sorted(f)] for f in fronts]
    if all(all(b < a for a, b in zip(v, v[1:])) for v in series):
        return "falls"
    if all(all(b > a for a, b in zip(v, v[1:])) for v in series):
        return "rises"
    return "mixed"


# ---------------------------------------------------------------- claims
R, I, A = "README.md", "docs/index.html", "ASSUMPTIONS.md"
TXT = r"([A-Za-z][A-Za-z ]*[A-Za-z])"  # a word or phrase captured for comparison, e.g. a walk-table action

# (label, {file: exact number of occurrences}, pattern, builder(c, captured) -> expected values).
# Patterns run on the document model of each file (prose on one line, then one line per table row, cells joined by " | "),
# with re.M, so ^ and $ anchor a table row. {N} is a complete, possibly signed, numeric token; units are part of the pattern.
CLAIMS = [
    # Level 1 recommendation
    ("recommended design, summary table", {R: 1},
     r"^\*\*Recommended design\*\* \| Straight copper fins: {N} mm channels, {N} mm fins, {N} mm tall\. {N} mK/W at {N} kPa$",
     lambda c, g: geom(rec(c)) + rp(rec(c))),
    ("recommended design, results table", {R: 1},
     r"^Recommended \(marginal return, meets the {N} mK/W objective\) \| {N} / {N} / {N} mm \| {N} mK/W \| {N} kPa$",
     lambda c, g: [objective_mk(c)] + geom(rec(c)) + rp(rec(c))),
    ("recommended design", {I: 1}, r"Straight copper fins with {N} mm channels, {N} mm fins, {N} mm tall\. The CFD gives {N} mK/W at "
     r"{N} kPa, which puts the plate base at {N} °C on average and {N} °C at the outlet end",
     lambda c, g: geom(rec(c)) + rp(rec(c)) + [rec(c)["cfd"]["base_avg_temp_c"], rec(c)["cfd"]["base_max_temp_c"]]),
    ("recommended design, figure alt text", {I: 1}, r"circled at {N} mK/W and {N} kPa", lambda c, g: rp(rec(c))),
    ("objective margin, introduction", {I: 1}, r"It meets the preferred thermal target with a {N}% margin, but",
     lambda c, g: [margin(c, rec(c)["cfd"]["r_th_k_per_w"])]),
    ("objective margin, recommended design", {I: 1}, r"This design meets the objective with a {N}% margin at rated flow, and misses it "
     r"at {N} LPM", lambda c, g: [margin(c, rec(c)["cfd"]["r_th_k_per_w"]),
                                 next(f["flow_lpm"] for f in c["l1"]["flow_sensitivity"]
                                      if f["r_th_k_per_w"] > c["l1"]["r_th_objective_k_per_w"])]),
    ("1D split of the recommended design", {I: 1}, r"The 1D model gives {N} mK/W for this design and splits it into about {N} for fins "
     r"and convection, {N} for the coolant warming up as it crosses the plate, and {N} for conduction through the base",
     lambda c, g: [1e3 * rec(c)["one_d"][k] for k in ("r_th_k_per_w", "r_conv", "r_caloric", "r_base")]),
    ("capacity of the recommended design", {R: 1}, r"carry at a {N} K rise \({N} kW for the recommended one\)",
     lambda c, g: [rise_k(c), rec(c)["ratios"]["capacity_kw"]]),
    ("capacity ratio: allowable rise", {I: 1}, r"Heat-removal capacity at the {N} K allowable base-to-inlet rise, per unit of pressure drop",
     lambda c, g: [rise_k(c)]),
    ("capacity ratio: trend along the fronts", {I: 1}, r"On every front in this study it (falls|rises) as pressure drop rises",
     lambda c, g: [capacity_direction(c)]),
    # flow sensitivity
    ("flow sensitivity", {R: 1}, r"averages {N} °C at the chip interface; at {N} LPM it rises to {N} mK/W \({N} °C\), at {N} LPM it falls "
     r"to {N} mK/W \({N} °C\)",
     lambda c, g: [flow(c, 1.5)["base_avg_temp_c"], 1.0, flow(c, 1.0)["r_th_k_per_w"] * 1e3, flow(c, 1.0)["base_avg_temp_c"],
                   2.0, flow(c, 2.0)["r_th_k_per_w"] * 1e3, flow(c, 2.0)["base_avg_temp_c"]]),
    ("flow sensitivity table row", {I: 3}, r"^{N} LPM \| {N} kPa \| {N} mK/W \| {N} °C \| {N} °C$", _flow_row),
    # counts
    ("designs screened and solved", {R: 1}, r"screens {N} fin geometries, a conjugate CFD model checks all {N} on the Pareto front",
     lambda c, g: [c["n_sweep"], n_cfd(c)]),
    ("designs screened, flowchart", {R: 1}, r"1D screen<br/>{N} designs", lambda c, g: [c["n_sweep"]]),
    ("designs solved, method", {R: 1}, r"within the budget \({N} here\)", lambda c, g: [n_cfd(c)]),
    ("designs screened and feasible", {I: 1}, r"It covers {N} designs in under 20 milliseconds, and {N} of them",
     lambda c, g: [c["n_sweep"], c["n_feasible"]]),
    ("designs solved", {I: 1}, r"All {N} designs on the 1D front", lambda c, g: [n_cfd(c)]),
    ("designs solved, parity caption", {I: 1}, r"1D against CFD for the {N} designs", lambda c, g: [n_cfd(c)]),
    ("designs solved, entrance", {I: 1}, r"Across the {N} designs the order", lambda c, g: [n_cfd(c)]),
    ("designs solved, wetted area", {I: 1}, r"For the {N} straight-fin designs checked with CFD", lambda c, g: [n_cfd(c)]),
    # parity: each occurrence anchored separately
    ("1D vs CFD, summary table", {R: 1},
     r"^\*\*1D screen vs CFD\*\* \| Thermal resistance agrees within {N}% on average, {N}% at worst \({N} designs\)$",
     lambda c, g: [c["l1"]["parity"]["mean_abs_pct"], c["l1"]["parity"]["max_abs_pct"], n_cfd(c)]),
    ("1D vs CFD, results section", {R: 1}, r"The 1D screen and the CFD agree on thermal resistance to {N}% on average and {N}% at worst\.",
     lambda c, g: [c["l1"]["parity"]["mean_abs_pct"], c["l1"]["parity"]["max_abs_pct"]]),
    ("1D vs CFD", {I: 1}, r"For these straight-fin designs the two models agree on thermal resistance to {N}% on average and {N}% at worst\.",
     lambda c, g: [c["l1"]["parity"]["mean_abs_pct"], c["l1"]["parity"]["max_abs_pct"]]),
    # reference designs, each row tied to its label
    ("reference: lowest ΔP meeting the objective", {R: 1},
     r"^Lowest ΔP meeting the {N} mK/W objective \| {N} / {N} / {N} mm \| {N} mK/W \| {N} kPa$", _reference("lowest_dp_meeting_objective")),
    ("reference: lowest ΔP meeting the threshold", {R: 1},
     r"^Lowest ΔP meeting the {N} mK/W threshold \| {N} / {N} / {N} mm \| {N} mK/W \| {N} kPa$", _reference("lowest_dp_meeting_threshold")),
    ("reference: lowest resistance within budget", {R: 1},
     r"^Lowest resistance within {N} kPa \| {N} / {N} / {N} mm \| {N} mK/W \| {N} kPa$", _reference("lowest_r_th_within_budget")),
    # selection walk and cut-off sensitivity
    ("walk table row", {I: 3}, r"^{N} / {N} / {N} mm \| {N} mK/W \| {N} kPa \| {N} \| " + TXT + "$", _walk_row),
    ("start of the walk", {I: 1}, r"the walk starts at {N} mm channels \({N} mK/W at {N} kPa\), the lowest pressure drop that meets the "
     r"20 mK/W limit, and moves through {N} mm channels with returns of {N} and {N}", _walk_start),
    ("best return at the stop", {I: 1}, r"In the last row, {N} is the best return", lambda c, g: [walk_return(c, c["l1"]["recommended"])]),
    ("cut-off in use", {I: 1}, r"each 1% of extra pressure drop buys at least {N}% lower resistance\. The {N} cut-off",
     lambda c, g: [c["l1"]["marginal_return_cutoff"]] * 2),
    ("cut-off sensitivity", {R: 1}, r"the pick holds from {N} to {N}; from {N} up it is the {N} mm design at {N} kPa, below {N} it moves "
     r"to higher pressure drop", lambda c, g: _cutoff_readme(c, g) + [c["l1"]["cutoff_stable_range"][0]]),
    ("walk summary", {R: 1}, r"The rule walks the front from {N} mm channels in four steps, with returns of {N}, {N}, {N} and finally {N} "
     r"from {N} mm to {N} mm channels \({N}% lower resistance for {N}× the pressure drop\); every further design returns {N} or less",
     _walk_summary),
    ("cut-off sensitivity", {I: 1}, r"the pick stays at {N} kPa for cut-offs of {N} and {N} and at {N} kPa \({N} mm channels, {N} mK/W\) "
     r"from {N} to {N}\. Below {N} it moves to {N} kPa and then {N} kPa", _cutoff_index),
    # entrance effects
    ("entrance-corrected design, summary table", {R: 1},
     r"entrance effects are quantified separately \({N} mK/W at {N} kPa, objective margin {N}%",
     lambda c, g: [_entrance_values(c)[1], _entrance_values(c)[3], margin(c, entrance_rec(c)["r_th_combined"])]),
    ("entrance-corrected design", {R: 1}, r"lowers R_th from {N} to {N} mK/W and raises the channel-core ΔP from {N} to {N} kPa, so the "
     r"objective margin grows from {N}% to {N}%",
     lambda c, g: _entrance_values(c) + [margin(c, entrance_rec(c)["r_th_fd"]), margin(c, entrance_rec(c)["r_th_combined"])]),
    ("entrance-corrected design", {I: 1}, r"the correction lowers the resistance from {N} to {N} mK/W and raises the pressure drop from "
     r"{N} to {N} kPa, so the margin on the objective grows to {N}%",
     lambda c, g: _entrance_values(c) + [margin(c, entrance_rec(c)["r_th_combined"])]),
    ("entrance check length", {I: 1}, r"also solved over its full {N} mm from a uniform inlet",
     lambda c, g: [c["en"]["settings"]["length_mm"]]),
    ("resolved entrance increments", {I: 1}, r"raises the length-mean heat transfer coefficient by {N}% and the channel pressure drop by {N}%",
     lambda c, g: _resolved_rec(c, "h_increase_resolved", "dp_increase_resolved")),
    ("entrance agreement, recommended", {I: 1}, r"matches the resulting changes in resistance and pressure drop to within {N}% and {N}%",
     lambda c, g: _resolved_rec(c, "r_th_contribution_disagreement", "dp_contribution_disagreement")),
    ("entrance agreement, widest channel", {I: 1}, r"the widest channel misses by {N}% on resistance",
     lambda c, g: [abs(100 * entrance_resolved(c, 0)["comparison"]["r_th_contribution_disagreement"])]),
    ("entrance-corrected selection", {I: 1}, r"would stop one step earlier, at {N} mm channels",
     lambda c, g: [geom(c["by_id"][c["en"]["selection_core"]["recommended"]])[0]]),
    ("entrance, earlier formulation", {I: 1}, r"An earlier version of the correction missed them by {N} to {N}%",
     lambda c, g: sorted(_resolved_rec(c, "earlier_r_th_contribution_disagreement", "earlier_dp_contribution_disagreement"))),
    ("entrance, last-step return", {I: 1}, r"the return on the last step of the selection falls from {N} to {N}",
     lambda c, g: last_step_returns(c)),
    # straight fins against a flat base; heated width; base
    ("wetted area ratio", {I: 1}, r"heated wetted area \(channel floor plus both fin sides\) is {N} to {N} times the heated footprint, and "
     r"{N} for the recommended design",
     lambda c, g: [min(wetted_ratios(c).values()), max(wetted_ratios(c).values()), wetted_ratios(c)[c["l1"]["recommended"]]]),
    ("heated width", {A: 1}, r"slightly under the {N} mm footprint: {N} mm \({N} × {N} mm\)",
     lambda c, g: [c["req"].footprint_width_m * 1e3, heated_width_mm(c), rec(c)["one_d"]["n_channels"],
                   heated_width_mm(c) / rec(c)["one_d"]["n_channels"]]),
    ("base thickness, figure caption", {I: 1}, r"The bottom {N} mm is the copper base", lambda c, g: [c["req"].base_thickness_m * 1e3]),
    ("base thickness, figure alt text", {I: 1}, r"between copper fins above a {N} mm base", lambda c, g: [c["req"].base_thickness_m * 1e3]),
    # slot jets
    ("best slot-jet array", {R: 1}, r"the best reaches {N} mK/W at {N} kPa", lambda c, g: best_jet(c)),
    ("slot-jet arrays", {R: 1}, r"\[comparison\]\([^)]*\): {N} slot-jet arrays", lambda c, g: [len(c["jets"])]),
    ("best slot-jet array", {I: 1}, r"the best of {N} arrays reached {N} mK/W", lambda c, g: [len(c["jets"]), best_jet(c)[0]]),
    # plan view
    ("plan-view designs solved", {R: 1}, r"\({N} plan-view designs solved", lambda c, g: [len(c["pv_by_id"])]),
    ("plan-view designs excluded", {R: 1}, r"The {N} designs whose fins are taller", lambda c, g: [len(c["pv"]["summary"]["excluded"])]),
    ("cross-topology pick", {R: 1}, r"pins of {N} mm diameter at {N} × {N} mm pitch \({N} mK/W at {N} kPa\) across all topologies",
     lambda c, g: [pv_params(c)[k] for k in ("d", "s_t", "s_l")] + pv_pick_rp(c)),
    ("straight-fin plan-view pick", {R: 1}, r"and {N} mm straight fins and channels \({N} mK/W at {N} kPa\) among straight fins",
     lambda c, g: [pv_params(c, "straight")["w"]] + pv_pick_rp(c, "straight")),
    ("offset-strip walk stop", {R: 1}, r"walk stops at its first point \({N} mK/W\)", lambda c, g: _offset_first(c)[:1]),
    ("offset-strip pick", {R: 1}, r"meets {N} mK/W \({N} mK/W at {N} kPa\)", lambda c, g: [objective_mk(c)] + pv_pick_rp(c, "offset_strip")),
    ("offset-strip walk stop", {I: 1}, r"walk stops at its first design \({N} mK/W at {N} kPa\), which misses the {N} mK/W objective",
     lambda c, g: _offset_first(c) + [objective_mk(c)]),
    ("pin field figure caption", {R: 1}, r"Staggered circular pins in plan view \({N} mm diameter, {N} × {N} mm pitch, {N} LPM\)",
     lambda c, g: [_pin_field(c)[k] for k in ("d", "s_t", "s_l")] + [c["pv"]["settings"]["flow_lpm"]]),
    ("pin field figure caption", {I: 1}, r"Circular pins, {N} mm across\.", lambda c, g: [_pin_field(c)["d"]]),
    ("plan-view table, circular pins", {I: 1}, r"^Circular pins \| {N} mm pins, {N} × {N} mm pitch \| {N} mK/W \| {N} kPa$",
     lambda c, g: [pv_params(c, "pin_circle")[k] for k in ("d", "s_t", "s_l")] + pv_pick_rp(c, "pin_circle")),
    ("plan-view table, straight", {I: 1}, r"^Straight \| {N} mm channels and fins \| {N} mK/W \| {N} kPa$",
     lambda c, g: [pv_params(c, "straight")["w"]] + pv_pick_rp(c, "straight")),
    ("plan-view table, zigzag", {I: 1}, r"^Zigzag \| {N} mm gap, {N} mm amplitude \| {N} mK/W \| {N} kPa$",
     lambda c, g: [pv_params(c, "zigzag")[k] for k in ("w", "amp")] + pv_pick_rp(c, "zigzag")),
    ("plan-view table, sinusoidal", {I: 1}, r"^Sinusoidal \| {N} mm gap, {N} mm amplitude \| {N} mK/W \| {N} kPa$",
     lambda c, g: [pv_params(c, "sinusoidal")[k] for k in ("w", "amp")] + pv_pick_rp(c, "sinusoidal")),
    ("plan-view table, diamond pins", {I: 1}, r"^Diamond pins \| {N} mm side \| {N} mK/W \| {N} kPa$",
     lambda c, g: [pv_params(c, "pin_diamond")["s"]] + pv_pick_rp(c, "pin_diamond")),
    ("plan-view table, offset strip", {I: 1}, r"^Offset strip \| {N} mm strips, {N} mm thick, {N} mm tall \| {N} mK/W \| {N} kPa$",
     lambda c, g: [pv_params(c, "offset_strip")[k] for k in ("l", "t", "h_fin")] + pv_pick_rp(c, "offset_strip")),
    ("against straight fins, circular pins", {I: 1}, r"Circular pins: {N} to {N}% lower resistance, between {N} and {N} kPa",
     lambda c, g: vs_straight(c, "pin_circle")),
    ("against straight fins, wavy", {I: 1}, r"Zigzag and sinusoidal: {N} to {N}% higher", lambda c, g: vs_straight(c, "zigzag", "sinusoidal")[:2]),
    ("against straight fins, diamond pins", {I: 1}, r"Diamond pins: {N} to {N}% higher", lambda c, g: vs_straight(c, "pin_diamond")[:2]),
    ("against straight fins, offset strips", {I: 1}, r"Offset strips: {N} to {N}% higher", lambda c, g: vs_straight(c, "offset_strip")[:2]),
    ("pin lead", {I: 1}, r"A {N} to {N}% lead for pins", lambda c, g: vs_straight(c, "pin_circle")[:2]),
    ("pin lead", {R: 1}, r"only by {N}–{N}% at equal pressure drop", lambda c, g: vs_straight(c, "pin_circle")[:2]),
    ("wavy-channel Reynolds numbers", {I: 1}, r"at Reynolds numbers of about {N} to {N}\.",
     lambda c, g: [("tens", f(x["result"]["reynolds_dh"] for x in c["pv_by_id"].values() if x["topology"] in ("sinusoidal", "zigzag")))
                   for f in (min, max)]),
    ("plan-view tall-fin bias", {R: 1}, r"reads ΔP {N}% and R_th {N}% below", lambda c, g: [bias(c)[1], bias(c)[0]]),
    ("plan-view tall-fin bias, headline", {R: 1}, r"inside the {N}% this model reads low", lambda c, g: [bias(c)[0]]),
    ("plan-view tall-fin bias", {I: 1}, r"gives resistance {N}% lower than the cross-section CFD", lambda c, g: [bias(c)[0]]),
    # validation; "within" figures are error bounds, rounded up
    ("duct toolchain error bound", {R: 1}, r"f·Re and Nu to within {N}%\.", lambda c, g: [up(duct_max_pct(c))]),
    ("production-solver test tolerances", {R: 1}, r"own tests are looser \({N}% on f·Re, {N}% on a limiting Nusselt number\)",
     lambda c, g: [100 * production_tolerances()["f_re_fanning"], 100 * production_tolerances()["nu_eff"]]),
    ("production-solver test tolerance", {I: 1}, r"whose own tests are looser \({N}% on a limiting Nusselt number\)",
     lambda c, g: [100 * production_tolerances()["nu_eff"]]),
    ("checks table: duct", {I: 1}, r"^Toolchain, fluid-only rectangular ducts at four aspect ratios [^|\n]* \| Friction and Nusselt number \| "
     r"Shah and London reference values \[4\] \| within {N}%$", lambda c, g: [up(duct_max_pct(c))]),
    ("checks table: cross-section mesh", {I: 1}, r"^Cross-section mesh at twice the resolution \| [^|\n]* \| Itself \| {N}% change$",
     lambda c, g: [up(c["l1"]["confidence"]["mesh_study"]["r_th_change_pct"])]),
    ("checks table: developing flow", {I: 1}, r"^Navier-Stokes solver, developing flow between plates \| [^|\n]* \| [^|\n]* \| within {N}%$",
     lambda c, g: [up(100 * abs(c["l2_entrance"]["nu_local_worst_error"]))]),
    ("checks table: plan-view solver", {I: 1}, r"^Plan-view solver, straight channel \| [^|\n]* \| [^|\n]* \| within {N}%$",
     lambda c, g: [up(100 * max(abs(c["g3"]["g1_dpdx_error"]), abs(c["g3"]["g1_nu_error"])))]),
    ("checks table: plan view against cross-section", {I: 1},
     r"^Plan view against cross-section, same straight design \| [^|\n]* \| [^|\n]* \| {N}% and {N}% low$", lambda c, g: bias(c)),
    ("checks table: energy balance", {I: 1}, r"^Energy balance, every Navier-Stokes run \| [^|\n]* \| [^|\n]* \| within {N}%$",
     lambda c, g: [up(ns_energy_max_pct(c))]),
]


# ---------------------------------------------------------------- document model
class _HtmlDoc(HTMLParser):
    """Visible prose (with image alt text) and table rows as lists of cell texts; tags and attributes are dropped."""
    BLOCK = {"p", "li", "h1", "h2", "h3", "figcaption", "figure", "div", "ol", "ul", "br", "hr", "article", "title"}

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.prose, self.rows, self.row, self.cell, self.skip, self.table = [], [], None, None, 0, 0

    def handle_starttag(self, tag, attrs):
        if tag in ("style", "script", "head"):
            self.skip += 1
        elif tag == "table":
            self.table += 1
        elif tag == "tr":
            self.row = []
        elif tag in ("td", "th"):
            self.cell = []
        elif tag == "img":
            self.prose.append(" " + (dict(attrs).get("alt") or "") + " ")
        if tag in self.BLOCK:
            (self.cell if self.cell is not None else self.prose).append(" ")

    def handle_endtag(self, tag):
        if tag in ("style", "script", "head"):
            self.skip -= 1
        elif tag == "table":
            self.table -= 1
        elif tag in ("td", "th") and self.cell is not None:
            self.row.append(" ".join("".join(self.cell).split()))
            self.cell = None
        elif tag == "tr" and self.row is not None:
            self.rows.append(self.row)
            self.row = None
        if tag in self.BLOCK and self.cell is None:
            self.prose.append(" ")

    def handle_data(self, data):
        if self.skip:
            return
        if self.cell is not None:
            self.cell.append(data)
        elif not self.table:
            self.prose.append(data)


def _markdown_doc(text: str) -> tuple[str, list[list[str]]]:
    prose, rows = [], []
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("|") and s.endswith("|") and len(s) > 1:
            cells = [x.strip() for x in s[1:-1].split("|")]
            if not all(re.fullmatch(r":?-{3,}:?", x) for x in cells):
                rows.append(cells)
        else:
            prose.append(line)
    return " ".join(" ".join(prose).split()), rows


def document(name: str, text: str) -> str:
    """Prose on the first line, then one line per table row with its cells joined by ' | '."""
    if name.endswith(".html"):
        p = _HtmlDoc()
        p.feed(text)
        prose, rows = " ".join("".join(p.prose).split()), p.rows
    else:
        prose, rows = _markdown_doc(text)
    return "\n".join([prose] + [" | ".join(r) for r in rows])


# ---------------------------------------------------------------- comparison and checks
def agrees(s: str, v) -> bool:
    """True when the printed token s equals v at the precision s is printed with (v = ("ceil", x) rounds x up; a str is
    compared as text). A signed token only agrees with a negative value."""
    if isinstance(v, str):
        return s == v
    mode, v = v if isinstance(v, tuple) else ("round", v)
    if mode == "tens":  # "about 90": the value rounded to the nearest ten
        return s == f"{round(v / 10) * 10:.0f}"
    s = s.replace("−", "-")
    if s.startswith("-") and (mode == "ceil" or v >= 0):
        return False
    if "," in s:
        return mode == "round" and f"{round(v):,}" == s
    d = len(s.split(".")[1]) if "." in s else 0
    if mode == "ceil":
        v = math.ceil(v * 10**d - 1e-9) / 10**d
    out = f"{v:.{d}f}"
    return out == s or (out == f"-{s}" and float(out) == 0)


def _claim_failures(docs: dict, c, covered: dict) -> list[str]:
    out = []
    for label, counts, pattern, build in CLAIMS:
        rx = re.compile(pattern.format(N=N), flags=re.M)
        for name, count in counts.items():
            matches = list(rx.finditer(docs[name]))
            covered[name] += [m.span(k) for m in matches for k in range(1, len(m.groups()) + 1)]
            if len(matches) != count:
                out.append(f"{name}: {label}: expected exactly {count} occurrence(s) of /{pattern}/, found {len(matches)} "
                           "(a number, its unit or the sentence around it changed, or the claim was deleted or repeated)")
            for m in matches:
                expected = build(c, m.groups())
                if expected is None:
                    out.append(f"{name}: {label}: '{m.group(0)}' names no design in the results")
                    continue
                if len(expected) != len(m.groups()):
                    out.append(f"{name}: {label}: checker built {len(expected)} values for {len(m.groups())} numbers")
                for k, (s, v) in enumerate(zip(m.groups(), expected), 1):
                    if not agrees(s, v):
                        modes = {"ceil": "rounded up to the printed precision", "tens": "rounded to the nearest ten"}
                        how, shown = ((modes[v[0]], v[1]) if isinstance(v, tuple) else ("as text", v) if isinstance(v, str)
                                      else ("rounded to the printed precision", v))
                        shown = shown if isinstance(shown, str) else f"{shown:.6g}"
                        out.append(f"{name}: {label}: item {k} reads {s} in '{m.group(0)}'; expected {shown} from the "
                                   f"results ({how})")
    return out


def _known_pairs(c) -> list[tuple[float, float]]:
    """(R_th in mK/W, ΔP in kPa) of every computed design or operating point the text may quote."""
    pairs = [tuple(rp(d)) for d in c["by_id"].values()]
    pairs += [(f["r_th_k_per_w"] * 1e3, f["dp_pa"] / 1e3) for f in c["l1"]["flow_sensitivity"]]
    pairs += [(r["cfd_anchored"]["r_th_combined"] * 1e3, r["cfd_anchored"]["dp_combined_core"] / 1e3) for r in c["en"]["one_d"]]
    pairs += [tuple(pv_rp(c, i)) for i in c["pv_by_id"]]
    pairs += [(j["result"]["r_th_k_per_w"] * 1e3, j["result"]["dp_pa"] / 1e3) for j in c["jets"]]
    return pairs


def _pair_failures(docs: dict, c) -> list[str]:
    """Every 'R mK/W at P kPa' pair, in prose or adjacent table cells, must be one computed design or operating point."""
    out, known = [], _known_pairs(c)
    for name in PUBLISHED:
        for m in re.finditer(r"{N} mK/W(?: at | and | \| ){N} kPa".format(N=N), docs[name]):
            r, dp = m.groups()
            if not any(agrees(r, kr) and agrees(dp, kp) for kr, kp in known):
                out.append(f"{name}: '{m.group(0)}' matches no computed design or operating point in the results")
    return out


def _geometry_row_failures(docs: dict, c) -> list[str]:
    """Every table row naming a straight-fin geometry 'w / t / h mm' carries that design's CFD R_th and ΔP."""
    out = []
    for name in PUBLISHED:
        for row in docs[name].split("\n")[1:]:
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


# Numbers in the published text that are not study results: requirement inputs are matched against the requirements
# file (value and unit); the rest are listed here with the reason they are not checked. Each entry is a regex for the
# context; every numeric token inside a match is accepted.
NOT_RESULTS = [
    (r"1D|2D|3D", "model names"),
    (r"\[\d+(?:, \d+)*\]", "citation numbers"),
    (r"issues? #?\d+|#\d+", "issue numbers"),
    (r"\b(?:19|20)\d\d\b", "years"),
    (r"November \d{4}", "date"),
    (r"PG\d+|\d+% propylene glycol", "coolant name"),
    (r"15 ms|under 20 milliseconds|about five seconds", "timings, machine dependent"),
    (r"a third of a 15 kPa whole-plate budget|one third of a 15 kPa", "budget allocation (ASSUMPTIONS R6)"),
    (r"1\.5 LPM per kW for a 10 K coolant rise", "OCP guidance [1]"),
    (r"under 20 mK/W case to inlet, interface material included, for a 1 kW module", "OCP figure [1]"),
    (r"reached 12 mK/W at 2\.6 LPM \[2\]\. 15 leaves about 5", "published plate [2] and the interface allowance"),
    (r"base of 2 to 3 mm assumed, the 1D model reads 6 to 16% low on thermal resistance from 1 LPM up", "issue 11 comparison, indicative"),
    (r"gaps of 0\.05 mm came close", "earlier slot-jet sweep (ASSUMPTIONS L8)"),
    (r"slots and gaps of 0\.1 mm or more|of at least 0\.1 mm", "slot-jet feature limit (ASSUMPTIONS L8)"),
    (r"a doubling of pressure drop must buy at least about 16% lower resistance", "meaning of the cut-off, arithmetic"),
    (r"less than 10% further in pressure drop", "selection guard (requirements file min_dp_ratio)"),
    (r"each 1% of extra pressure drop", "definition of the return"),
    (r"for cut-offs from 0\.05 to 0\.60|re-runs the rule from 0\.05 to 0\.60|cut-off from 0\.05 to 0\.60|from 0\.05 to 0\.60",
     "cut-off sweep range"),
    (r"a cut-off of 0\.20 picks a denser pin array and 0\.30 a sparser one|at 0\.20 and at 0\.30", "cut-off sensitivity, plan view"),
    (r"\b\d\. \*\*[A-Z]", "list numbering"),
    (r"Level \d", "fidelity level names"),
    (r"Scientific Reports \d+:\d+|arXiv:\d+\.\d+", "bibliographic numbers"),
    (r"about 0\.7 at 0\.2 kPa, 0\.28–0\.41 from 0\.7 to 1\.6 kPa, mostly below 0\.2 beyond", "local-return figure description"),
    (r"stopped at 1\.29 kPa because it had no design between 1\.29 and 2\.28 kPa", "earlier coarser grid, historical"),
    (r"falls from about 0\.7 near 0\.2 kPa to between 0\.28 and 0\.41 from 0\.7 to 1\.6 kPa, and is mostly below the 0\.25 cut-off",
     "local-return figure alt text"),
    (r"within a 5% band", "parity figure band"),
    (r"\b5% band|±5%", "parity figure band"),
    (r"six fin shapes|Six fin shapes|across all six|five designs|Five small panels", "counts in words"),
    (r"0\.15 mm pin or wavy features|the 0\.15 mm gap between strip rows", "manufacturing limit"),
    (r"about ten times taller", "aspect ratio in words"),
    (r"returns at least 0\.25 in log-log terms", "cut-off in use"),
]


def _requirement_tokens(c) -> set:
    """(value, unit) pairs the requirements file sets; a token with that unit and value at its precision is accepted."""
    q = c["req"]
    vals = {("kW", q.heat_load_w / 1e3), ("W", q.heat_load_w), ("mm", q.footprint_width_m * 1e3), ("mm", q.flow_length_m * 1e3),
            ("×", q.footprint_width_m * 1e3), ("°C", q.inlet_temp_c), ("LPM", flow(c, 1.5)["flow_lpm"]), ("kPa", q.dp_budget_pa / 1e3),
            ("mK/W", q.r_th_threshold_k_per_w * 1e3), ("mK/W", q.r_th_objective_k_per_w * 1e3), ("mm", q.base_thickness_m * 1e3),
            ("", q.marginal_return_cutoff), ("cut-off", q.marginal_return_cutoff)}
    vals |= {("LPM", f["flow_lpm"]) for f in c["l1"]["flow_sensitivity"]}
    import yaml
    cfg = yaml.safe_load((ROOT / "examples" / "gpu_coldplate.yaml").read_text())["constraints"]
    vals |= {("mm", cfg["min_feature_mm"]), ("×", cfg["max_aspect_ratio"]), ("times", cfg["max_aspect_ratio"])}
    return vals


UNIT = r"\s?(mK/W|kPa|kW|W|LPM|°C|mm|K|%|×|times|cut-off)?"


def _coverage_failures(docs: dict, c, covered: dict) -> list[str]:
    """Every number in the published prose and tables is checked by a claim, is a requirement input, or is listed in
    NOT_RESULTS with a reason. A new or reworded number therefore fails until a check covers it."""
    out, reqs = [], _requirement_tokens(c)
    for name in PUBLISHED:
        doc = docs[name]
        masked = re.sub(r"```.*?```|`[^`]*`|\]\([^)]*\)|https?://\S+", lambda m: " " * len(m.group(0)), doc, flags=re.S)
        allowed = [m.span() for p, _ in NOT_RESULTS for m in re.finditer(p, masked)]
        spans = covered[name] + allowed
        for m in re.finditer(N + UNIT, masked):
            a, b = m.span(1)
            if any(x <= a and b <= y for x, y in spans):
                continue
            unit = m.group(2) or ""
            if any(u == unit and agrees(m.group(1), v) for u, v in reqs):
                continue
            ctx_text = doc[max(0, a - 60):b + 40].replace("\n", " / ")
            out.append(f"{name}: number {m.group(1)}{(' ' + unit) if unit else ''} is not checked against the results: '…{ctx_text}…'")
    return out


def failures(texts: dict) -> list[str]:
    """Every disagreement between the published text (file name -> raw text) and the committed results."""
    c = ctx()
    docs = {name: document(name, t) for name, t in texts.items()}
    covered = {name: [] for name in docs}
    out = _claim_failures(docs, c, covered) + _pair_failures(docs, c) + _geometry_row_failures(docs, c)
    out += _coverage_failures(docs, c, covered)
    for name in PUBLISHED:
        for value in SUPERSEDED:
            if re.search(rf"(?<![\d.]){re.escape(value)}(?![\d])", docs[name]):
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
    assert re.findall(N, "114.0 and 14.0, 88.9% and 8.9%, 4,485, -1.0 and −1.0, PG25, 1D") == [
        "114.0", "14.0", "88.9", "8.9", "4,485", "-1.0", "−1.0"]
    assert agrees("14.0", 14.0487) and not agrees("114.0", 14.0487) and agrees("14.05", 14.0487) and not agrees("14.04", 14.0487)
    assert agrees("0.008", ("ceil", 0.00735)) and not agrees("0.007", ("ceil", 0.00735))
    assert not agrees("-1.0", 1.03) and not agrees("−1.0", 1.03) and agrees("-1.0", -1.03)


def test_production_tolerances_are_parsed():
    tol = production_tolerances()
    assert 0 < tol["nu_eff"] < 1 and 0 < tol["f_re_fanning"] < 1


def test_document_model_reads_cells_not_markup():
    d = document("x.html", '<table><tr><td class="n a" style="x">1.0 <b>kPa</b></td><td>b</td></tr></table><p>p <img alt="q"></p>')
    assert d.split("\n") == ["p q", "1.0 kPa | b"]
    assert document("x.md", "a\n\n|  x | y |\n|---|---|\n| 1 |  2 |\n").split("\n") == ["a", "x | y", "1 | 2"]


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
