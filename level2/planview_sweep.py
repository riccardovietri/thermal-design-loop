"""Level 2 plan view: checks on one design per topology, sweeps at rated flow, and the cross-topology comparison.

    uv run python -m level2.planview_sweep              # checks + sweeps + plot (about 10-20 minutes)
    uv run python -m level2.planview_sweep --checks     # Stage 2 only: mesh, periodicity, developing strip, field figures
    uv run python -m level2.planview_sweep --sweep      # Stage 3 only (keeps saved checks)
    uv run python -m level2.planview_sweep --extend     # solve only designs not yet in results.json, then summarise
    uv run python -m level2.planview_sweep --plot-only  # redraw from results.json
    uv run python -m level2.planview_sweep --fields     # re-solve the six field-figure designs and redraw their figures

All designs at 1.5 LPM PG25, 1 kW over 60 x 60 mm; compared at equal total flow and footprint, never at equal Re.
Results: examples/output/level2_planview/results.json (written after every design); per-design records in
runs/level2_planview/<design_id>/results.json.
"""
import argparse
import json
import sys
import time
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
from matplotlib.lines import Line2D

from tools import plotstyle as ps
from tools.fluids import fluid_props
from tools.pareto import capacity_ratios, pareto_indices, recommend_marginal_return
from tools.report import draw_requirements
from tools.types import Requirements

from .planar import NewtonDivergence
from .planforms import GeometryRejected, build
from .planview import describe, design_id, reynolds_b, run_design, run_developing

REQ = Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, 40.0, 5e3, "PG25", r_th_threshold_k_per_w=20e-3,
                   r_th_objective_k_per_w=15e-3)
FLOW_LPM = 1.5
CELLS, CELLS_FINE, CELLS_DEVELOPING = 10.0, 20.0, 8.0
OUT = Path("examples/output/level2_planview")
RUNS = Path("runs/level2_planview")
TOPOLOGIES = ("straight", "sinusoidal", "zigzag", "offset_strip", "pin_circle", "pin_diamond")
LABELS = {"straight": "Straight", "sinusoidal": "Sinusoidal", "zigzag": "Zigzag", "offset_strip": "Offset strip",
          "pin_circle": "Pins, circular", "pin_diamond": "Pins, diamond"}
COLOURS = {"straight": "#374151", "sinusoidal": "#0d9488", "zigzag": "#65a30d", "offset_strip": "#2563eb",
           "pin_circle": "#7c3aed", "pin_diamond": "#db2777"}
MARKERS = {"straight": "o", "sinusoidal": "s", "zigzag": "^", "offset_strip": "D", "pin_circle": "P", "pin_diamond": "X"}

# One representative design per topology for Stage 2, and the number of developing-strip periods to solve.
REPRESENTATIVE = {
    "straight": ({"w": 0.4, "t": 0.4, "h_fin": 4.0}, 12),
    "sinusoidal": ({"w": 0.4, "t": 0.4, "amp": 0.3, "lam": 2.0, "h_fin": 4.0}, 14),
    "zigzag": ({"w": 0.4, "t": 0.4, "amp": 0.3, "lam": 2.0, "h_fin": 4.0}, 14),
    "offset_strip": ({"l": 1.0, "t": 0.4, "p_t": 1.2, "g": 0.0, "h_fin": 4.0}, 10),
    "pin_circle": ({"d": 0.4, "s_t": 0.7, "s_l": 0.6, "h_fin": 4.0}, 10),
    "pin_diamond": ({"s": 0.4, "s_t": 0.9, "s_l": 0.7, "h_fin": 4.0}, 10),
}


def sweep_designs() -> dict[str, list[dict]]:
    """About ten designs per topology inside (or, where noted, at) the manufacturing limits; H 3 and 4 mm."""
    d = {t: [] for t in TOPOLOGIES}
    for w in (0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7):
        d["straight"].append({"w": w, "t": 0.4, "h_fin": 4.0})
    for w in (0.15, 0.2, 0.3, 0.4, 0.5):
        d["straight"].append({"w": w, "t": 0.3, "h_fin": 3.0})
    waves = [(0.4, 0.15, 2.0), (0.4, 0.3, 2.0), (0.4, 0.3, 4.0), (0.4, 0.6, 4.0), (0.3, 0.15, 2.0), (0.3, 0.3, 2.0),
             (0.3, 0.6, 4.0), (0.2, 0.3, 2.0), (0.5, 0.3, 2.0)]
    for topo in ("sinusoidal", "zigzag"):
        for w, amp, lam in waves:
            d[topo].append({"w": w, "t": 0.4, "amp": amp, "lam": lam, "h_fin": 4.0})
        d[topo].append({"w": 0.3, "t": 0.3, "amp": 0.3, "lam": 2.0, "h_fin": 3.0})
    for s in (0.15, 0.25, 0.4):
        for l in (0.5, 1.0, 2.0):
            d["offset_strip"].append({"l": l, "t": 0.4, "p_t": round(2 * (0.4 + s), 4), "g": 0.0, "h_fin": 4.0})
    d["offset_strip"].append({"l": 1.0, "t": 0.3, "p_t": 0.9, "g": 0.0, "h_fin": 3.0})
    for s_t in (0.6, 0.8, 1.0):
        for s_l in (0.5, 0.7, 0.9):
            d["pin_circle"].append({"d": 0.4, "s_t": s_t, "s_l": s_l, "h_fin": 4.0})
    d["pin_circle"] += [{"d": 0.6, "s_t": 0.9, "s_l": 0.7, "h_fin": 4.0}, {"d": 0.3, "s_t": 0.5, "s_l": 0.45, "h_fin": 3.0}]
    for s_t in (0.75, 0.9, 1.1):
        for s_l in (0.5, 0.7, 0.9):
            d["pin_diamond"].append({"s": 0.4, "s_t": s_t, "s_l": s_l, "h_fin": 4.0})
    d["pin_diamond"].append({"s": 0.3, "s_t": 0.6, "s_l": 0.45, "h_fin": 3.0})
    # Added 2026-11-11 (SPEC Amendment 3, item 6): H / min gap <= 10 left these topologies with fewer than 5 eligible designs.
    d["offset_strip"] += [{"l": l, "t": 0.3, "p_t": 1.2, "g": 0.0, "h_fin": 3.0} for l in (0.5, 1.0, 2.0)]
    d["offset_strip"] += [{"l": l, "t": 0.2, "p_t": 0.8, "g": 0.0, "h_fin": 2.0} for l in (0.5, 1.0)]
    d["pin_circle"] += [{"d": 0.4, "s_t": 0.9, "s_l": s_l, "h_fin": 4.0} for s_l in (0.7, 0.9)]
    d["pin_circle"] += [{"d": 0.3, "s_t": 0.6, "s_l": s_l, "h_fin": 3.0} for s_l in (0.55, 0.7)]
    d["pin_diamond"] += [{"s": 0.4, "s_t": 1.0, "s_l": s_l, "h_fin": 4.0} for s_l in (0.7, 0.9)]
    d["pin_diamond"] += [{"s": 0.4, "s_t": 1.2, "s_l": 0.6, "h_fin": 4.0}]
    d["pin_diamond"] += [{"s": 0.3, "s_t": 0.75, "s_l": s_l, "h_fin": 3.0} for s_l in (0.5, 0.7)]
    d["pin_diamond"] += [{"s": 0.3, "s_t": 0.85, "s_l": 0.5, "h_fin": 3.0}]
    return d


MAX_H_OVER_GAP = 10.0  # Level 1 skiving limit H / gap <= 10, applied to plan-view results (SPEC Amendment 3, item 6)


def aspect_exclusion(row: dict) -> str | None:
    """Reason a converged design is excluded by H / minimum fluid gap > 10, else None."""
    h, gap = row["params_mm"]["h_fin"], row["result"]["geometry"]["min_gap_mm"]
    return f"H / min gap = {h:g} / {gap:.3f} = {h / gap:.1f} > {MAX_H_OVER_GAP:g}" if h / gap > MAX_H_OVER_GAP + 1e-9 else None


def _load() -> dict:
    path = OUT / "results.json"
    return json.loads(path.read_text()) if path.exists() else {}


def _save(data: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "results.json").write_text(json.dumps(data, indent=2, default=float))


def _field_record(prod) -> dict:
    """Numbers quoted in the captions of the field figures (summary.md)."""
    return {"design_id": prod.design_id, "params_mm": prod.params_mm, "reynolds_dh": prod.reynolds_dh, "dh_mm": prod.dh_mm,
            "dpdx_pa_per_m": prod.dpdx_pa_per_m, "h_w_m2k": prod.h_w_m2k, "dp_pa": prod.dp_pa, "r_th_k_per_w": prod.r_th_k_per_w,
            "wall_flux_w_m2": REQ.heat_load_w / (prod.a_floor_m2 + prod.fin_efficiency * prod.a_fin_m2), "n_unknowns": prod.n_unknowns}


def _field_design(fluid, topo: str):
    """Solve the representative design of one topology at production settings and draw its field figure."""
    prod = run_design(build(topo, REPRESENTATIVE[topo][0]), REQ, fluid, FLOW_LPM, cells=CELLS, out_dir=OUT / "fields" / topo, plot=True)
    (OUT / "fields" / topo / "fields.png").rename(OUT / f"fields_{topo}.png")
    return prod


def _clean_fields_dir() -> None:
    for p in (OUT / "fields").glob("*/results.json"):
        p.unlink()
    for p in (OUT / "fields").glob("*"):
        p.rmdir()
    (OUT / "fields").rmdir()


def run_fields(fluid) -> dict:
    """Re-solve only the six field-figure designs (fields are not saved) and write fields.json for the captions."""
    rec = {}
    for topo in REPRESENTATIVE:
        prod = _field_design(fluid, topo)
        rec[topo] = _field_record(prod)
        print(f"[fields] {topo:13s} dP {prod.dp_pa / 1e3:.2f} kPa R {prod.r_th_k_per_w * 1e3:.2f} mK/W -> fields_{topo}.png", flush=True)
    _clean_fields_dir()
    (OUT / "fields.json").write_text(json.dumps(rec, indent=2))
    return rec


def run_checks(fluid, data: dict) -> dict:
    checks, fields = {}, {}
    for topo, (params, n_dev) in REPRESENTATIVE.items():
        pf = build(topo, params)
        t0 = time.perf_counter()
        prod = _field_design(fluid, topo)
        fields[topo] = _field_record(prod)
        fine = run_design(pf, REQ, fluid, FLOW_LPM, cells=CELLS_FINE)
        two = run_design(pf, REQ, fluid, FLOW_LPM, cells=CELLS, n_periods=2)
        pf_dev = replace(pf, period_mm=5.0) if topo == "straight" else pf  # straight: 5 mm reporting windows
        dev = run_developing(pf_dev, REQ, fluid, n_dev, FLOW_LPM, cells=CELLS_DEVELOPING)
        within = [p["period"] for p in dev["periods"] if abs(p["h_w_m2k"] / prod.h_w_m2k - 1) < 0.02
                  and abs(p["dpdx_pa_per_m"] / prod.dpdx_pa_per_m - 1) < 0.02]
        w1, w2 = two.windows
        c = {
            "design_id": prod.design_id, "params_mm": params, "reynolds_dh": prod.reynolds_dh, "dh_mm": prod.dh_mm,
            "reynolds_pin_max": prod.reynolds_pin_max,
            "mesh": {"cells": [CELLS, CELLS_FINE], "n_unknowns": [prod.n_unknowns, fine.n_unknowns],
                     "dpdx_change_pct": 100 * (fine.dpdx_pa_per_m / prod.dpdx_pa_per_m - 1),
                     "h_change_pct": 100 * (fine.h_w_m2k / prod.h_w_m2k - 1),
                     "r_th_change_pct": 100 * (fine.r_th_k_per_w / prod.r_th_k_per_w - 1),
                     "energy_balance_error": [prod.energy_balance_error, fine.energy_balance_error],
                     "undershoot_rel": [prod.undershoot_rel, fine.undershoot_rel]},
            "periodic_two_period_cell": {"dpdx_change_pct": 100 * (w2["dpdx_pa_per_m"] / w1["dpdx_pa_per_m"] - 1),
                                         "h_change_pct": 100 * (w2["h_w_m2k"] / w1["h_w_m2k"] - 1),
                                         "vs_one_period_h_pct": 100 * (two.h_w_m2k / prod.h_w_m2k - 1)},
            "developing": {**{k: v for k, v in dev.items() if k != "periods"},
                           "h_over_periodic": [p["h_w_m2k"] / prod.h_w_m2k for p in dev["periods"]],
                           "dpdx_over_periodic": [p["dpdx_pa_per_m"] / prod.dpdx_pa_per_m for p in dev["periods"]],
                           "first_period_within_2pct": min(within) if within else None,
                           "period_mm": pf_dev.period_mm},
            "seconds": time.perf_counter() - t0,
        }
        c["pass"] = (abs(c["mesh"]["dpdx_change_pct"]) < 2 and abs(c["mesh"]["h_change_pct"]) < 2
                     and abs(dev["last_two_dpdx_change_pct"]) < 2 and abs(dev["last_two_h_change_pct"]) < 2)
        checks[topo] = c
        print(f"[check] {topo:13s} Re_Dh {prod.reynolds_dh:6.1f}  mesh dp/dx {c['mesh']['dpdx_change_pct']:+.3f}% "
              f"h {c['mesh']['h_change_pct']:+.3f}%  2-period cell h {c['periodic_two_period_cell']['h_change_pct']:+.3f}%  "
              f"developing last two: dp/dx {dev['last_two_dpdx_change_pct']:+.2f}% h {dev['last_two_h_change_pct']:+.2f}%, "
              f"within 2% of periodic from period {c['developing']['first_period_within_2pct']} "
              f"({pf_dev.period_mm:g} mm each); undershoot {prod.undershoot_rel:+.1e}  {c['seconds']:.0f}s", flush=True)
        data["checks"] = checks
        _save(data)
    _clean_fields_dir()
    (OUT / "fields.json").write_text(json.dumps(fields, indent=2))
    return checks


def run_sweep(fluid, data: dict, extend: bool = False) -> list[dict]:
    """Solve every sweep design; with extend, keep converged entries already in results.json and solve only new ones."""
    done = {r["design_id"]: r for r in data.get("designs", []) if r["status"] == "converged"} if extend else {}
    rows = []
    for topo, designs in sweep_designs().items():
        for params in designs:
            entry = {"topology": topo, "params_mm": params, "flow_lpm": FLOW_LPM}
            prev = next((r for r in done.values() if r["topology"] == topo and r["params_mm"] == params), None)
            if prev is not None:
                rows.append(prev)
                continue
            try:
                pf = build(topo, params)
            except GeometryRejected as err:
                entry.update(status="rejected", reasons=err.reasons, design_id=design_id(build(topo, params, strict=False)))
                rows.append(entry)
                print(f"[sweep] {entry['design_id']}: REJECTED {err}", flush=True)
                data["designs"] = rows
                _save(data)
                continue
            did = design_id(pf)
            try:
                r = run_design(pf, REQ, fluid, FLOW_LPM, cells=CELLS, out_dir=RUNS / did)
                entry.update(status="converged", design_id=did, result=asdict(r))
                print(f"[sweep] {did:55s} Re_Dh {r.reynolds_dh:6.1f} dP {r.dp_pa / 1e3:6.2f} kPa R {r.r_th_k_per_w * 1e3:6.2f} mK/W "
                      f"h {r.h_w_m2k / 1e3:5.1f} eta {r.fin_efficiency:.2f} energy {100 * r.energy_balance_error:+.3f}% "
                      f"{r.seconds:.0f}s", flush=True)
            except NewtonDivergence as err:
                entry.update(status="newton_failed", design_id=did, reynolds_b=reynolds_b(REQ, fluid, FLOW_LPM, params["h_fin"]),
                             failed_at_reynolds=err.reynolds, last_step=err.last_step)
                print(f"[sweep] {did}: NEWTON FAILED {err}", flush=True)
            rows.append(entry)
            data["designs"] = rows
            _save(data)
    return rows


def _converged(rows):
    return [r for r in rows if r["status"] == "converged"]


def _brief(r: dict) -> dict:
    return {"design_id": r["design_id"], "topology": r["topology"], "r_th_k_per_w": r["result"]["r_th_k_per_w"],
            "dp_pa": r["result"]["dp_pa"]}


def _select(rows: list[dict]) -> dict:
    """Marginal-return rule (tools.pareto.recommend_marginal_return) over eligible converged designs, with design ids."""
    sel = recommend_marginal_return([(r["result"]["dp_pa"], r["result"]["r_th_k_per_w"]) for r in rows], REQ.dp_budget_pa,
                                    REQ.r_th_threshold_k_per_w, REQ.r_th_objective_k_per_w, REQ.marginal_return_cutoff)
    rid = lambda i: rows[i]["design_id"]
    return {**_brief(rows[sel["index"]]), "rule": sel["rule"], "meets_threshold": sel["meets_threshold"],
            "meets_objective": sel["meets_objective"], "front": [rid(i) for i in sel["front"]],
            "steps": [{"from": rid(st["from"]), "to": rid(st["to"]) if st["to"] is not None else None,
                       "returns": [{"to": rid(e["to"]), "return": e["return"]} for e in st["returns"]]} for st in sel["steps"]]}


def summarise(data: dict) -> dict:
    rows = data["designs"]
    conv_all = _converged(rows)
    excluded = [{"design_id": r["design_id"], "topology": r["topology"], "reason": aspect_exclusion(r)}
                for r in conv_all if aspect_exclusion(r)]
    ok = [r for r in conv_all if not aspect_exclusion(r)]  # eligible: converged and H / gap <= 10
    budget, threshold, objective = REQ.dp_budget_pa, REQ.r_th_threshold_k_per_w, REQ.r_th_objective_k_per_w
    flow_m3_s = FLOW_LPM / 6e4
    ratios = {r["design_id"]: capacity_ratios(r["result"]["dp_pa"], r["result"]["r_th_k_per_w"], REQ.heat_load_w, flow_m3_s, threshold)
              for r in ok}
    per_topo = {}
    for topo in TOPOLOGIES:
        sel = [r for r in rows if r["topology"] == topo]
        conv = [r for r in ok if r["topology"] == topo]
        inb = [r for r in conv if r["result"]["dp_pa"] <= budget]
        best = min(inb, key=lambda r: r["result"]["r_th_k_per_w"]) if inb else None
        low = lambda lim: min((r for r in inb if r["result"]["r_th_k_per_w"] <= lim), key=lambda r: r["result"]["dp_pa"], default=None)
        per_topo[topo] = {
            "n_designs": len(sel), "n_converged": len(_converged(sel)), "n_rejected": sum(r["status"] == "rejected" for r in sel),
            "n_newton_failed": sum(r["status"] == "newton_failed" for r in sel),
            "n_excluded_h_over_gap": sum(e["topology"] == topo for e in excluded), "n_eligible": len(conv),
            "reynolds_dh_range": [min(r["result"]["reynolds_dh"] for r in conv), max(r["result"]["reynolds_dh"] for r in conv)] if conv else None,
            "recommended": _select(conv) if inb else None,
            "best_in_budget": _brief(best) if best else None,
            "lowest_dp_meeting_threshold": _brief(low(threshold)) if low(threshold) else None,
            "lowest_dp_meeting_objective": _brief(low(objective)) if low(objective) else None,
        }
    equal_dp = []
    for level in (1e3, 2e3, 3e3, 4e3, 5e3):
        def best_at(cond):
            c = [r for r in ok if cond(r) and r["result"]["dp_pa"] <= level]
            return min(c, key=lambda r: r["result"]["r_th_k_per_w"]) if c else None
        s, n = best_at(lambda r: r["topology"] == "straight"), best_at(lambda r: r["topology"] != "straight")
        equal_dp.append({"dp_max_pa": level, "best_straight": _brief(s) if s else None, "best_other": _brief(n) if n else None})
    # every eligible non-straight design on its topology's front within budget vs the eligible straight front (log-log interpolation)
    st = [r for r in ok if r["topology"] == "straight"]
    sf = [st[i] for i in pareto_indices([r["result"]["dp_pa"] for r in st], [r["result"]["r_th_k_per_w"] for r in st])]
    lx = np.log([r["result"]["dp_pa"] for r in sf])
    ly = np.log([r["result"]["r_th_k_per_w"] for r in sf])
    vs_straight = []
    for topo in TOPOLOGIES[1:]:
        sel = [r for r in ok if r["topology"] == topo]
        front = [sel[i] for i in pareto_indices([r["result"]["dp_pa"] for r in sel], [r["result"]["r_th_k_per_w"] for r in sel])]
        for r in front:
            dp = r["result"]["dp_pa"]
            if dp <= budget and lx[0] <= np.log(dp) <= lx[-1]:
                r_st = float(np.exp(np.interp(np.log(dp), lx, ly)))
                vs_straight.append({"design_id": r["design_id"], "topology": topo, "dp_pa": dp,
                                    "r_th_k_per_w": r["result"]["r_th_k_per_w"], "straight_r_th_at_same_dp": r_st,
                                    "difference_pct": 100 * (r["result"]["r_th_k_per_w"] / r_st - 1)})
    return {"requirements": {"dp_budget_pa": budget, "r_th_threshold_k_per_w": threshold, "r_th_objective_k_per_w": objective,
                             "marginal_return_cutoff": REQ.marginal_return_cutoff, "max_h_over_gap": MAX_H_OVER_GAP},
            "recommended": _select(ok), "per_topology": per_topo, "excluded": excluded, "ratios": ratios,
            "vs_straight_front": vs_straight, "equal_dp": equal_dp,
            "max_energy_balance_error": max(abs(r["result"]["energy_balance_error"]) for r in conv_all),
            "worst_undershoot_rel": min(r["result"]["undershoot_rel"] for r in conv_all),
            "total_seconds": sum(r["result"]["seconds"] for r in conv_all)}


X_LIM, Y_LIM = (0.3, 1.5 * REQ.dp_budget_pa / 1e3), (8.0, 1.5 * REQ.r_th_threshold_k_per_w * 1e3)
VIEW_X, VIEW_Y = (0.3, 7.5), (9.5, 25.0)  # drawn window; every front point lies inside it
# Overview labels placed in clear spots next to their lines (x kPa, y mK/W, ha, va); the rest go in the right margin.
INLINE_LABELS = {"straight": (0.37, 21.3, "left", "bottom"), "offset_strip": (2.5, 14.0, "left", "bottom"),
                 "pin_circle": (2.3, 10.6, "left", "top")}


def _xy(rows):
    return [r["result"]["dp_pa"] / 1e3 for r in rows], [r["result"]["r_th_k_per_w"] * 1e3 for r in rows]


def _front(rows):
    x, y = _xy(rows)
    idx = sorted(pareto_indices(x, y), key=lambda i: x[i])
    return [x[i] for i in idx], [y[i] for i in idx]


def _log_axes(ax, x_ticks=None, y_ticks=None):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*VIEW_X)
    ax.set_ylim(*VIEW_Y)
    ps.plain_log_ticks(ax, *VIEW_X, "x", ticks=x_ticks)
    ps.plain_log_ticks(ax, *VIEW_Y, "y", ticks=y_ticks)


def _spread(values: list[float], min_gap: float) -> list[float]:
    """Move sorted label positions apart to at least min_gap, keeping their mean."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    pos = [values[i] for i in order]
    for _ in range(200):
        moved = False
        for k in range(1, len(pos)):
            if pos[k] - pos[k - 1] < min_gap:
                d = (min_gap - (pos[k] - pos[k - 1])) / 2
                pos[k - 1] -= d
                pos[k] += d
                moved = True
        if not moved:
            break
    out = [0.0] * len(values)
    for k, i in enumerate(order):
        out[i] = pos[k]
    return out


def _eligible_and_excluded(data):
    conv = _converged(data["designs"])
    return [r for r in conv if not aspect_exclusion(r)], [r for r in conv if aspect_exclusion(r)]


def _outside(rows) -> list[str]:
    outside = []
    for topo in TOPOLOGIES:
        x, y = _xy([r for r in rows if r["topology"] == topo])
        n_out = sum(not (X_LIM[0] <= a <= X_LIM[1] and Y_LIM[0] <= b <= Y_LIM[1]) for a, b in zip(x, y))
        if n_out:
            outside.append(f"{LABELS[topo]}: {n_out} of {len(x)} eligible outside the window")
    return outside


def plot(data: dict, out: Path) -> tuple[Path, list[str]]:
    """Overview: each topology's Pareto front over its eligible designs, labelled directly, and the cross-topology pick."""
    rows, _ = _eligible_and_excluded(data)
    budget = REQ.dp_budget_pa / 1e3
    rec = data["summary"]["recommended"]
    fig, ax = ps.figure(5.8)
    draw_requirements(ax, VIEW_X[0], VIEW_Y[0], budget, REQ.r_th_threshold_k_per_w * 1e3, REQ.r_th_objective_k_per_w * 1e3,
                      VIEW_Y[1], labels_right=True)
    ends = {}
    for topo in TOPOLOGIES:
        fx, fy = _front([r for r in rows if r["topology"] == topo])
        if not fx:
            continue
        ax.plot(fx, fy, color=COLOURS[topo], marker=MARKERS[topo], ms=6, lw=2.0 if topo == "straight" else 1.8,
                mec="white", mew=0.6, zorder=3)
        ends[topo] = (fx[-1], fy[-1])
    # direct labels: in clear spots inside the axes where a line has one, otherwise in the right margin at the height of the
    # line's last point (pushed apart where they would overlap)
    for t, (x, y, ha, va) in INLINE_LABELS.items():
        if t in ends:
            ax.text(x, y, LABELS[t], color=COLOURS[t], fontsize=ps.SMALL + 0.5, ha=ha, va=va, zorder=6)
    names = [t for t in ends if t not in INLINE_LABELS]
    ly = _spread([np.log10(ends[t][1]) for t in names], 0.032)
    for t, y in zip(names, ly):
        ax.annotate(LABELS[t], xy=ends[t], xycoords="data", xytext=(1.02, 10 ** y), textcoords=("axes fraction", "data"),
                    color=COLOURS[t], fontsize=ps.SMALL + 0.5, va="center", ha="left", annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color=COLOURS[t], lw=0.8, alpha=0.6, shrinkA=0, shrinkB=4))
    r = next(r for r in rows if r["design_id"] == rec["design_id"])
    px, py = r["result"]["dp_pa"] / 1e3, r["result"]["r_th_k_per_w"] * 1e3
    ax.scatter([px], [py], s=260, facecolors="none", edgecolors=ps.INK, linewidths=2.0, zorder=5)
    pick = describe(r["topology"], r["params_mm"]).split(", ", 1)
    ax.annotate(f"Pick: {pick[0].lower()}\n{pick[1]}\n{py:.2f} mK/W at {px:.2f} kPa", xy=(px, py),
                xytext=(0.62, 10.9), textcoords="data", fontsize=ps.SMALL + 0.5, color=ps.INK, ha="center", va="center",
                arrowprops=dict(arrowstyle="-", color=ps.INK, lw=1.0, shrinkA=2, shrinkB=9))
    _log_axes(ax, x_ticks=(0.3, 0.5, 1, 2, 3, 5, 7))
    ax.set_xlabel("Pressure drop ΔP (kPa)")
    ax.set_ylabel("Thermal resistance R_th (mK/W)")
    ax.set_title("Pareto front of each fin planform at equal flow")
    ps.save(fig, out)
    return out, _outside(rows)


def plot_panels(data: dict, out: Path) -> Path:
    """One panel per non-straight topology: its eligible designs and front in colour against the straight-fin front in grey."""
    rows, excl = _eligible_and_excluded(data)
    summ = data["summary"]
    budget = REQ.dp_budget_pa / 1e3
    sx, sy = _front([r for r in rows if r["topology"] == "straight"])
    fig, axes = ps.figure(6.4, 2, 3, sharex=True, sharey=True)
    for ax, topo in zip(axes.flat, TOPOLOGIES[1:]):
        ax.axvline(budget, color=ps.ACCENT, ls="--", lw=1.0, zorder=1)
        ax.axhline(REQ.r_th_threshold_k_per_w * 1e3, color=ps.THRESHOLD, ls="--", lw=1.0, zorder=1)
        ax.axhline(REQ.r_th_objective_k_per_w * 1e3, color=ps.OBJECTIVE, ls=":", lw=1.2, zorder=1)
        ax.plot(sx, sy, color=ps.FIN_GREY, marker="o", ms=4, lw=2.2, zorder=2)
        ex, ey = _xy([r for r in excl if r["topology"] == topo])
        ax.scatter(ex, ey, s=16, facecolors="none", edgecolors=COLOURS[topo], linewidths=0.7, alpha=0.35, zorder=2)
        sel = [r for r in rows if r["topology"] == topo]
        x, y = _xy(sel)
        ax.scatter(x, y, s=30, color=COLOURS[topo], marker=MARKERS[topo], edgecolors="white", linewidths=0.5, zorder=3)
        fx, fy = _front(sel)
        ax.plot(fx, fy, color=COLOURS[topo], lw=1.8, zorder=3)
        pk = summ["per_topology"][topo]["recommended"]
        if pk:
            ax.scatter([pk["dp_pa"] / 1e3], [pk["r_th_k_per_w"] * 1e3], s=150, facecolors="none", edgecolors=ps.INK,
                       linewidths=1.6, zorder=4)
        ax.set_title(f"{LABELS[topo]} ({len(sel)} designs)", fontsize=11.5, loc="left", color=COLOURS[topo])
        _log_axes(ax, x_ticks=(0.5, 1, 2, 5), y_ticks=(10, 15, 20))
    axes[0, 2].tick_params(labelbottom=True)
    for ax in (*axes[1, :2], axes[0, 2]):
        ax.set_xlabel("ΔP (kPa)")
    for ax in axes[:, 0]:
        ax.set_ylabel("R_th (mK/W)")
    key = axes.flat[-1]
    key.axis("off")
    handles = [
        Line2D([], [], color="#6b7280", marker="P", ms=7, lw=1.8, label="Planform front"),
        Line2D([], [], color=ps.FIN_GREY, marker="o", ms=4, lw=2.2, label="Straight-fin front"),
        Line2D([], [], ls="none", marker="o", ms=11, mfc="none", mec=ps.INK, mew=1.6, label="Planform's pick"),
        Line2D([], [], ls="none", marker="o", ms=5, mfc="none", mec="#6b7280", alpha=0.5, label="Too tall for its gap"),
        Line2D([], [], color=ps.THRESHOLD, ls="--", lw=1.0, label=f"Threshold {REQ.r_th_threshold_k_per_w * 1e3:g} mK/W"),
        Line2D([], [], color=ps.OBJECTIVE, ls=":", lw=1.2, label=f"Objective {REQ.r_th_objective_k_per_w * 1e3:g} mK/W"),
        Line2D([], [], color=ps.ACCENT, ls="--", lw=1.0, label=f"ΔP budget {budget:g} kPa"),
    ]
    key.legend(handles=handles, loc="center left", fontsize=ps.SMALL + 0.5, borderaxespad=0, handlelength=2.0)
    fig.suptitle("Each fin planform against straight fins")
    ps.save(fig, out)
    return out


def write_summary(data: dict, out: Path, outside: list[str] | None = None) -> Path:
    """Short markdown summary of the plan-view selection: per-topology and cross-topology picks with the rule's steps."""
    summ, rows = data["summary"], {r["design_id"]: r for r in _converged(data["designs"])}
    outside = outside or []
    req = summ["requirements"]
    cutoff = req["marginal_return_cutoff"]
    pt = lambda did: f"{did} ({rows[did]['result']['r_th_k_per_w']*1e3:.2f} mK/W at {rows[did]['result']['dp_pa']/1e3:.2f} kPa)"

    def steps(sel):
        out = []
        for k, st in enumerate(sel["steps"], 1):
            rets = "; ".join((f"**{e['return']:.3f}**" if e["return"] >= cutoff else f"{e['return']:.3f}") + f" to `{e['to']}`"
                             for e in st["returns"]) or "none"
            out.append(f"{k}. From {pt(st['from'])}: returns {rets} → " + (f"move to `{st['to']}`." if st["to"] else "stop."))
        return out

    lines = ["# Plan-view topology comparison: selection", "",
             f"Requirements (SPEC Amendment 3): ΔP ≤ {req['dp_budget_pa']/1e3:g} kPa, threshold R_th ≤ {req['r_th_threshold_k_per_w']*1e3:g} mK/W "
             f"(must meet), objective ≤ {req['r_th_objective_k_per_w']*1e3:g} mK/W (want), marginal-return cut-off {cutoff:g}, "
             f"H / min gap ≤ {req['max_h_over_gap']:g}. All numbers are plan-view (tall-fin limit), which reads R_th about 6% low on "
             "straight fins.", "",
             "Rule: from the lowest-ΔP point of the eligible front (within budget and threshold), move to the nearest further point "
             "whose return e = −ln(R_j/R_i) / ln(ΔP_j/ΔP_i) is at least the cut-off; stop when none is. "
             "Capacity = threshold × heat load / R_th; kW per pumping W = heat load / (ΔP × flow). Both ratios always favour the lowest "
             "ΔP inside the window, which is why they are reported and not used to choose.", "",
             "![Pareto front of each fin planform](pareto_topologies.png)", "",
             f"Each line is one planform's Pareto front over its eligible designs, at {FLOW_LPM:g} LPM PG25 at 40 °C, 1 kW over "
             "60 × 60 mm, copper, 2 mm base, features and gaps ≥ 0.15 mm, H / t ≤ 10 and H / gap ≤ 10. The ring marks the design "
             "the rule picks across all planforms. Not drawn: designs off the fronts, and the "
             f"{len(summ['excluded'])} designs excluded because their fins are taller than 10 × the narrowest gap. Shading: within "
             "the budget and the threshold (light), within the budget and the objective (darker). The model is 2D Navier-Stokes in "
             "the plane of the fin layer (tall-fin limit, developed periodic flow) plus a resistance network. It has no floor or "
             "cover: on straight fins it reads ΔP 7% and R_th 6% below the Level 1 cross-section CFD, so a gap between planforms "
             "smaller than that is not resolved." + (" " + "; ".join(outside) + "." if outside else ""), "",
             "![Each fin planform against straight fins](pareto_topologies_panels.png)", "",
             "One panel per planform: its eligible designs (markers) and Pareto front in colour, the straight-fin front in grey, "
             "and the planform's own pick (ring). Faint hollow markers are designs of that planform excluded by H / gap > 10. "
             "Same axes, budget, threshold and objective in every panel.", "",
             "## Per topology", "",
             "| Topology | Eligible / converged | Excluded (H/gap > 10) | Pick | R_th (mK/W) | ΔP (kPa) | Objective met | Capacity (kW) | kW per kPa | kW per pumping W |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for t in TOPOLOGIES:
        p = summ["per_topology"][t]
        sel = p["recommended"]
        if sel is None:
            lines.append(f"| {LABELS[t]} | {p['n_eligible']} / {p['n_converged']} | {p['n_excluded_h_over_gap']} | none in budget | | | | | | |")
            continue
        rt = summ["ratios"][sel["design_id"]]
        lines.append(f"| {LABELS[t]} | {p['n_eligible']} / {p['n_converged']} | {p['n_excluded_h_over_gap']} | `{sel['design_id']}` "
                     f"(`{sel['rule']}`) | {sel['r_th_k_per_w']*1e3:.2f} | {sel['dp_pa']/1e3:.2f} | {'yes' if sel['meets_objective'] else 'no'} | "
                     f"{rt['capacity_kw']:.2f} | {rt['capacity_kw_per_kpa']:.2f} | {rt['kw_per_pumping_w']:.1f} |")
    for t in TOPOLOGIES:
        sel = summ["per_topology"][t]["recommended"]
        if sel:
            lines += ["", f"### {LABELS[t]}", ""] + steps(sel)
    rec = summ["recommended"]
    lines += ["", "## Across all topologies", "",
              f"Pick: `{rec['design_id']}` ({LABELS[rec['topology']]}), {rec['r_th_k_per_w']*1e3:.2f} mK/W at {rec['dp_pa']/1e3:.2f} kPa, "
              f"rule `{rec['rule']}`; meets the objective: {'yes' if rec['meets_objective'] else 'no'}.", ""] + steps(rec)
    lines += ["", "## Against the straight-fin front at equal ΔP (eligible designs)", "",
              "| Design | ΔP (kPa) | R_th (mK/W) | Straight at same ΔP | Difference |", "|---|---|---|---|---|"]
    lines += [f"| `{v['design_id']}` | {v['dp_pa']/1e3:.2f} | {v['r_th_k_per_w']*1e3:.2f} | {v['straight_r_th_at_same_dp']*1e3:.2f} | "
              f"{v['difference_pct']:+.1f}% |" for v in summ["vs_straight_front"]]
    lines += ["", "## Excluded by H / min gap > 10 (solved, kept in results.json, not in fronts or selection)", "",
              "| Design | Reason |", "|---|---|"] + [f"| `{e['design_id']}` | {e['reason']} |" for e in summ["excluded"]]
    fields_path = OUT / "fields.json"
    if fields_path.exists():
        fields = json.loads(fields_path.read_text())
        lines += ["", "## Field plots", "",
                  "One design per planform, solved as a developed periodic cell and tiled over a few periods; fins are grey. Top: "
                  "speed with streamlines. Bottom: temperature above the local mean bulk temperature (velocity-weighted) under a "
                  f"uniform wall heat flux, the periodic part only. All at {FLOW_LPM:g} LPM PG25. Pin pitch is transverse × streamwise."]
        for topo, f in fields.items():
            lines += ["", f"![{describe(topo, f['params_mm'])}](fields_{topo}.png)", "",
                      f"{describe(topo, f['params_mm'])}, {f['params_mm']['h_fin']:g} mm tall: Re_Dh = {f['reynolds_dh']:.0f} "
                      f"(Dh = {f['dh_mm']:.3f} mm), dp/dx = {f['dpdx_pa_per_m'] / 1e3:.1f} kPa/m, h = {f['h_w_m2k'] / 1e3:.1f} kW/m²K, "
                      f"wall flux {f['wall_flux_w_m2'] / 1e4:.1f} W/cm² → ΔP {f['dp_pa'] / 1e3:.2f} kPa, "
                      f"R_th {f['r_th_k_per_w'] * 1e3:.1f} mK/W."]
    out.write_text("\n".join(lines) + "\n")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Level 2 plan-view checks, sweeps and topology comparison.")
    ap.add_argument("--checks", action="store_true")
    ap.add_argument("--sweep", action="store_true")
    ap.add_argument("--extend", action="store_true", help="solve only sweep designs not already converged in results.json")
    ap.add_argument("--plot-only", action="store_true")
    ap.add_argument("--fields", action="store_true", help="re-solve the six field-figure designs only, then redraw everything")
    a = ap.parse_args()
    fluid = fluid_props(REQ.coolant, REQ.inlet_temp_c)
    data = _load()
    t0 = time.perf_counter()
    if a.extend:
        run_sweep(fluid, data, extend=True)
    elif a.fields:
        run_fields(fluid)
    elif not a.plot_only:
        data.setdefault("settings", {})
        data["settings"] = {"flow_lpm": FLOW_LPM, "cells": CELLS, "cells_fine": CELLS_FINE, "cells_developing": CELLS_DEVELOPING,
                            "reference_length_mm": 1.0, "velocity_scale": "superficial velocity Q / (W H)",
                            "re_dh": "rho (U_s / porosity) Dh / mu, Dh = 4 A_fluid / P_wetted (plan)"}
        if a.checks or not a.sweep:
            run_checks(fluid, data)
        if a.sweep or not a.checks:
            run_sweep(fluid, data)
        data["wall_seconds_last_run"] = time.perf_counter() - t0
    if "designs" in data:
        data["summary"] = summarise(data)
        _save(data)
        out, outside = plot(data, OUT / "pareto_topologies.png")
        print(out, plot_panels(data, OUT / "pareto_topologies_panels.png"), write_summary(data, OUT / "summary.md", outside))
        rec = data["summary"]["recommended"]
        print(f"cross-topology pick {rec['design_id']} ({rec['rule']}): {rec['r_th_k_per_w']*1e3:.2f} mK/W at {rec['dp_pa']/1e3:.2f} kPa")
        for o in outside:
            print("outside window:", o)
    return 0


if __name__ == "__main__":
    sys.exit(main())
