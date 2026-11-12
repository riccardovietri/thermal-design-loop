"""Level 2 plan view: checks on one design per topology, sweeps at rated flow, and the cross-topology comparison.

    uv run python -m level2.planview_sweep              # checks + sweeps + plot (about 10-20 minutes)
    uv run python -m level2.planview_sweep --checks     # Stage 2 only: mesh, periodicity, developing strip, field figures
    uv run python -m level2.planview_sweep --sweep      # Stage 3 only (keeps saved checks)
    uv run python -m level2.planview_sweep --plot-only  # redraw from results.json

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

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tools.fluids import fluid_props
from tools.pareto import pareto_indices, recommend
from tools.report import plain_log_ticks
from tools.types import Requirements

from .planar import NewtonDivergence
from .planforms import GeometryRejected, build
from .planview import design_id, reynolds_b, run_design, run_developing

REQ = Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, 40.0, 6e3, "PG25", r_th_target_k_per_w=15e-3)
FLOW_LPM = 1.5
CELLS, CELLS_FINE, CELLS_DEVELOPING = 10.0, 20.0, 8.0
OUT = Path("examples/output/level2_planview")
RUNS = Path("runs/level2_planview")
TOPOLOGIES = ("straight", "sinusoidal", "zigzag", "offset_strip", "pin_circle", "pin_diamond")
LABELS = {"straight": "Straight", "sinusoidal": "Sinusoidal", "zigzag": "Zigzag", "offset_strip": "Offset strip",
          "pin_circle": "Pins, circular", "pin_diamond": "Pins, diamond"}
COLOURS = {"straight": "#1d4ed8", "sinusoidal": "#0d9488", "zigzag": "#65a30d", "offset_strip": "#c2410c",
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
    return d


def _load() -> dict:
    path = OUT / "results.json"
    return json.loads(path.read_text()) if path.exists() else {}


def _save(data: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "results.json").write_text(json.dumps(data, indent=2, default=float))


def run_checks(fluid, data: dict) -> dict:
    checks = {}
    for topo, (params, n_dev) in REPRESENTATIVE.items():
        pf = build(topo, params)
        t0 = time.perf_counter()
        prod = run_design(pf, REQ, fluid, FLOW_LPM, cells=CELLS, out_dir=OUT / "fields" / topo, plot=True)
        (OUT / "fields" / topo / "fields.png").rename(OUT / f"fields_{topo}.png")
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
    for p in (OUT / "fields").glob("*/results.json"):
        p.unlink()
    for p in (OUT / "fields").glob("*"):
        p.rmdir()
    (OUT / "fields").rmdir()
    return checks


def run_sweep(fluid, data: dict) -> list[dict]:
    rows = []
    for topo, designs in sweep_designs().items():
        for params in designs:
            entry = {"topology": topo, "params_mm": params, "flow_lpm": FLOW_LPM}
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


def summarise(data: dict) -> dict:
    rows = data["designs"]
    ok = _converged(rows)
    budget, target = REQ.dp_budget_pa, REQ.r_th_target_k_per_w
    pick = recommend([(r["result"]["dp_pa"], r["result"]["r_th_k_per_w"]) for r in ok], budget, target)
    rec = ok[pick["index"]]
    per_topo = {}
    for topo in TOPOLOGIES:
        sel = [r for r in rows if r["topology"] == topo]
        conv = _converged(sel)
        inb = [r for r in conv if r["result"]["dp_pa"] <= budget]
        best = min(inb, key=lambda r: r["result"]["r_th_k_per_w"]) if inb else None
        meets = [r for r in inb if r["result"]["r_th_k_per_w"] <= target]
        low = min(meets, key=lambda r: r["result"]["dp_pa"]) if meets else None
        per_topo[topo] = {
            "n_designs": len(sel), "n_converged": len(conv), "n_rejected": sum(r["status"] == "rejected" for r in sel),
            "n_newton_failed": sum(r["status"] == "newton_failed" for r in sel),
            "reynolds_dh_range": [min(r["result"]["reynolds_dh"] for r in conv), max(r["result"]["reynolds_dh"] for r in conv)] if conv else None,
            "best_in_budget": {"design_id": best["design_id"], "r_th_k_per_w": best["result"]["r_th_k_per_w"],
                               "dp_pa": best["result"]["dp_pa"]} if best else None,
            "meets_target": bool(meets),
            "lowest_dp_meeting_target": {"design_id": low["design_id"], "r_th_k_per_w": low["result"]["r_th_k_per_w"],
                                         "dp_pa": low["result"]["dp_pa"]} if low else None,
        }
    equal_dp = []
    for level in (1e3, 2e3, 4e3, 6e3):
        def best_at(cond):
            c = [r for r in ok if cond(r) and r["result"]["dp_pa"] <= level]
            return min(c, key=lambda r: r["result"]["r_th_k_per_w"]) if c else None
        s, n = best_at(lambda r: r["topology"] == "straight"), best_at(lambda r: r["topology"] != "straight")
        equal_dp.append({"dp_max_pa": level,
                         "best_straight": {"design_id": s["design_id"], "r_th_k_per_w": s["result"]["r_th_k_per_w"],
                                           "dp_pa": s["result"]["dp_pa"]} if s else None,
                         "best_other": {"design_id": n["design_id"], "r_th_k_per_w": n["result"]["r_th_k_per_w"],
                                        "dp_pa": n["result"]["dp_pa"]} if n else None})
    # every non-straight design on its topology's front within budget vs the straight front interpolated (log-log) at its dP
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
    return {"vs_straight_front": vs_straight, "recommended": {"design_id": rec["design_id"], "topology": rec["topology"], "rule": pick["rule"],
                            "dp_pa": rec["result"]["dp_pa"], "r_th_k_per_w": rec["result"]["r_th_k_per_w"]},
            "per_topology": per_topo, "equal_dp": equal_dp,
            "max_energy_balance_error": max(abs(r["result"]["energy_balance_error"]) for r in ok),
            "worst_undershoot_rel": min(r["result"]["undershoot_rel"] for r in ok),
            "total_seconds": sum(r["result"]["seconds"] for r in ok)}


X_LIM, Y_LIM = (0.3, 9.0), (8.0, 30.0)


def plot(data: dict, out: Path) -> tuple[Path, list[str]]:
    rows = _converged(data["designs"])
    level1 = json.loads(Path("examples/output/gpu_coldplate/results.json").read_text())["designs"]
    l1 = [(d["cfd"]["dp_pa"] / 1e3, d["cfd"]["r_th_k_per_w"] * 1e3) for d in level1]
    f1 = pareto_indices([p[0] for p in l1], [p[1] for p in l1])
    budget, target = REQ.dp_budget_pa / 1e3, REQ.r_th_target_k_per_w * 1e3
    rec = data["summary"]["recommended"]
    outside = []
    fig, ax = plt.subplots(figsize=(10.5, 6.0), constrained_layout=True)
    ax.fill_between([X_LIM[0], budget], Y_LIM[0], target, color="#16a34a", alpha=0.08, lw=0, zorder=0,
                    label="Meets budget and target")
    ax.plot([l1[i][0] for i in f1], [l1[i][1] for i in f1], color="#111827", lw=0.8, ls=":", zorder=1,
            label="Straight fins, Level 1 CFD front (floor + cover)")
    for topo in TOPOLOGIES:
        sel = [r for r in rows if r["topology"] == topo]
        if not sel:
            continue
        x = [r["result"]["dp_pa"] / 1e3 for r in sel]
        y = [r["result"]["r_th_k_per_w"] * 1e3 for r in sel]
        n_out = sum(not (X_LIM[0] <= a <= X_LIM[1] and Y_LIM[0] <= b <= Y_LIM[1]) for a, b in zip(x, y))
        if n_out:
            outside.append(f"{LABELS[topo]}: {n_out} of {len(sel)} outside the plotted window")
        ax.scatter(x, y, s=34, color=COLOURS[topo], marker=MARKERS[topo], edgecolors="white", linewidths=0.5, zorder=3, alpha=0.9)
        front = pareto_indices(x, y)
        ax.plot([x[i] for i in front], [y[i] for i in front], color=COLOURS[topo], lw=1.6, zorder=2,
                label=f"{LABELS[topo]} (plan view, {len(sel)} designs)")
    r = next(r for r in rows if r["design_id"] == rec["design_id"])
    ax.scatter([r["result"]["dp_pa"] / 1e3], [r["result"]["r_th_k_per_w"] * 1e3], s=170, facecolors="none",
               edgecolors="#111827", linewidths=1.4, zorder=4, label=f"Recommended ({rec['rule']})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*X_LIM)
    ax.set_ylim(*Y_LIM)
    plain_log_ticks(ax, *X_LIM, "x")
    plain_log_ticks(ax, *Y_LIM, "y")
    ax.axvline(budget, color="#c2410c", ls="--", lw=1)
    ax.text(budget * 0.97, Y_LIM[1] * 0.97, f"ΔP budget {budget:g} kPa", color="#c2410c", ha="right", va="top", fontsize=9)
    ax.axhline(target, color="#15803d", ls="--", lw=1)
    ax.text(X_LIM[0] * 1.04, target * 0.98, f"R_th target {target:g} mK/W", color="#15803d", ha="left", va="top", fontsize=9)
    ax.set_xlabel("Channel-core pressure drop ΔP (kPa)")
    ax.set_ylabel("Thermal resistance R_th (mK/W)")
    fig.suptitle("Fin planforms at equal flow: 2D plan-view Navier-Stokes (tall-fin limit, developed periodic) + resistance network\n"
                 f"{FLOW_LPM} LPM PG25 at 40 °C, 1 kW over 60 × 60 mm, copper, 2 mm base; features and gaps ≥ 0.15 mm, H/t ≤ 10",
                 fontsize=10.5, x=0.02, ha="left")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    note = "Lines: each topology's own Pareto front.\nPlan view ignores floor and cover:\nfor straight fins it reads ΔP 7% and\nR_th 6% below Level 1 CFD."
    if outside:
        note += "\n" + "\n".join(outside)
    ax.text(1.02, 0.02, note, transform=ax.transAxes, fontsize=8, va="bottom", ha="left", color="#374151")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out, outside


def main() -> int:
    ap = argparse.ArgumentParser(description="Level 2 plan-view checks, sweeps and topology comparison.")
    ap.add_argument("--checks", action="store_true")
    ap.add_argument("--sweep", action="store_true")
    ap.add_argument("--plot-only", action="store_true")
    a = ap.parse_args()
    fluid = fluid_props(REQ.coolant, REQ.inlet_temp_c)
    data = _load()
    t0 = time.perf_counter()
    if not a.plot_only:
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
        print(out)
        print(json.dumps(data["summary"]["recommended"], indent=2))
        for o in outside:
            print("outside window:", o)
    return 0


if __name__ == "__main__":
    sys.exit(main())
