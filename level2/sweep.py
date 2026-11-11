"""Level 2 sweep: slot-jet arrays at rated flow over manufacturable slot width, gap and pitch, compared with the
Level 1 straight-fin designs at the same flow; plus a mesh-independence study on one slot-jet case.

    uv run python -m level2.sweep                   # mesh study, 18-case sweep, comparison plot
    uv run python -m level2.sweep --mesh-study-only
"""
import argparse
import itertools
import json
import sys
import time
from dataclasses import asdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tools.fluids import fluid_props
from tools.pareto import pareto_indices
from tools.report import plain_log_ticks
from tools.types import Requirements

from .slot_jet import SlotJetDesign, run_slot_jet

REQ = Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, 40.0, 6e3, "PG25", r_th_target_k_per_w=15e-3)
# Minimum jet gap and slot width 0.1 mm (ASSUMPTIONS.md L8); rated flow 1.5 LPM (R4).
SLOTS_MM, GAPS_MM, PITCHES_MM, FLOW_LPM = (0.1, 0.2), (0.1, 0.15, 0.2), (1.0, 2.0, 4.0), 1.5
MESH_CASE_MM, MESH_LEVELS, PRODUCTION_CELLS = (0.2, 0.1, 2.0), (10, 20, 40), 20
OUT = Path("examples/output/level2_slot_jet")


def mesh_study(fluid) -> dict:
    slot, gap, pitch = MESH_CASE_MM
    d = SlotJetDesign(slot / 1e3, gap / 1e3, pitch / 1e3, FLOW_LPM)
    levels = []
    for cells in MESH_LEVELS:
        t0 = time.perf_counter()
        r = run_slot_jet(d, REQ, fluid, Path("runs/level2/_mesh_study") / f"cells_{cells:02d}", cells_across=cells)
        levels.append({"cells_across": cells, "n_unknowns": r.n_unknowns, "r_th_k_per_w": r.r_th_k_per_w, "dp_pa": r.dp_pa,
                       "energy_balance_error": r.energy_balance_error, "seconds": time.perf_counter() - t0})
        print(f"mesh {cells:2d} cells: {r.n_unknowns} unknowns, R={r.r_th_k_per_w*1e3:.4f} mK/W, dP={r.dp_pa:.3f} Pa", flush=True)
    prod = next(lv for lv in levels if lv["cells_across"] == PRODUCTION_CELLS)
    fine = levels[-1]
    study = {"design": asdict(d), "production_cells_across": PRODUCTION_CELLS, "levels": levels,
             "r_th_change_pct_20_to_40": 100 * (fine["r_th_k_per_w"] / prod["r_th_k_per_w"] - 1),
             "dp_change_pct_20_to_40": 100 * (fine["dp_pa"] / prod["dp_pa"] - 1)}
    (OUT / "mesh_study.json").write_text(json.dumps(study, indent=2))
    print(f"mesh study 20 -> 40 cells: R_th {study['r_th_change_pct_20_to_40']:+.3f}%, ΔP {study['dp_change_pct_20_to_40']:+.3f}%")
    return study


def sweep(fluid) -> list[dict]:
    rows = []
    for slot, gap, pitch in itertools.product(SLOTS_MM, GAPS_MM, PITCHES_MM):
        d = SlotJetDesign(slot / 1e3, gap / 1e3, pitch / 1e3, FLOW_LPM)
        t0 = time.perf_counter()
        try:
            r = run_slot_jet(d, REQ, fluid, Path("runs/level2") / d.design_id, cells_across=PRODUCTION_CELLS)
        except RuntimeError as err:
            print(f"{d.design_id}: FAILED {err}", flush=True)
            continue
        rows.append({"design": asdict(d), "result": asdict(r), "seconds": time.perf_counter() - t0})
        print(f"{d.design_id}: Re={r.reynolds_slot:5.1f} dP={r.dp_pa/1e3:8.3f} kPa R={r.r_th_k_per_w*1e3:7.2f} mK/W "
              f"energy err={100*r.energy_balance_error:+.3f}% {rows[-1]['seconds']:.0f}s", flush=True)
        (OUT / "results.json").write_text(json.dumps(rows, indent=2))
    return rows


def plot_compare(rows: list[dict], out: Path) -> Path:
    level1 = json.loads(Path("examples/output/gpu_coldplate/results.json").read_text())["designs"]
    fins = sorted(((d["cfd"]["dp_pa"] / 1e3, d["cfd"]["r_th_k_per_w"] * 1e3) for d in level1))
    jets = [(r["result"]["dp_pa"] / 1e3, r["result"]["r_th_k_per_w"] * 1e3) for r in rows]
    budget, target = REQ.dp_budget_pa / 1e3, REQ.r_th_target_k_per_w * 1e3
    xs, ys = [p[0] for p in fins + jets], [p[1] for p in fins + jets]
    x_lo, x_hi = 0.7 * min(xs), 1.5 * max(max(xs), budget)
    y_lo, y_hi = 0.7 * min(ys + [target]), 1.3 * max(ys)

    fig, ax = plt.subplots(figsize=(10.5, 5.6), constrained_layout=True)
    ax.fill_between([x_lo, budget], y_lo, target, color="#16a34a", alpha=0.08, lw=0, zorder=0, label="Meets budget and target")
    ax.plot([p[0] for p in fins], [p[1] for p in fins], "o-", color="#1d4ed8", label=f"Straight fins (Level 1 CFD, {FLOW_LPM} LPM)")
    for gap, colour in zip(GAPS_MM, ("#7c2d12", "#c2410c", "#f59e0b")):
        for slot, marker in zip(SLOTS_MM, ("D", "s")):
            sel = [r for r in rows if abs(r["design"]["gap_m"] * 1e3 - gap) < 1e-9 and abs(r["design"]["slot_width_m"] * 1e3 - slot) < 1e-9]
            ax.scatter([r["result"]["dp_pa"] / 1e3 for r in sel], [r["result"]["r_th_k_per_w"] * 1e3 for r in sel], s=46, color=colour,
                       marker=marker, edgecolors="#374151", linewidths=0.6, zorder=3,
                       label=f"Slot jets, {gap:.2f} mm gap, {slot:.1f} mm slot")
    front = pareto_indices([r["result"]["dp_pa"] for r in rows], [r["result"]["r_th_k_per_w"] for r in rows])
    ax.plot([rows[i]["result"]["dp_pa"] / 1e3 for i in front], [rows[i]["result"]["r_th_k_per_w"] * 1e3 for i in front],
            color="#c2410c", lw=1, ls="--", label="Slot-jet Pareto front")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    plain_log_ticks(ax, x_lo, x_hi, "x")
    plain_log_ticks(ax, y_lo, y_hi, "y")
    ax.axvline(budget, color="#c2410c", ls="--", lw=1)
    ax.text(budget * 0.95, y_hi * 0.97, f"ΔP budget {budget:g} kPa", color="#c2410c", ha="right", va="top", fontsize=9)
    ax.axhline(target, color="#15803d", ls="--", lw=1)
    ax.text(x_lo * 1.05, target * 0.98, f"R_th target {target:g} mK/W", color="#15803d", ha="left", va="top", fontsize=9)
    ax.set_xlabel("Pressure drop ΔP (kPa)")
    ax.set_ylabel("Thermal resistance R_th (mK/W)")
    ax.set_title(f"Slot-jet arrays on a flat base vs straight fins\n{FLOW_LPM} LPM PG25, 1 kW, 60 × 60 mm; jet gap and slot ≥ 0.1 mm, fins and channels ≥ 0.15 mm", fontsize=11)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Level 2 slot-jet mesh study, sweep and comparison plot.")
    ap.add_argument("--mesh-study-only", action="store_true")
    ap.add_argument("--plot-only", action="store_true", help="redraw the comparison plot from the saved results.json")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    fluid = fluid_props(REQ.coolant, REQ.inlet_temp_c)
    if a.plot_only:
        rows = json.loads((OUT / "results.json").read_text())
    else:
        mesh_study(fluid)
        if a.mesh_study_only:
            return 0
        rows = sweep(fluid)
    print(plot_compare(rows, OUT / "pareto_compare.png"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
