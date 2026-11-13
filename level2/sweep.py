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

from matplotlib.lines import Line2D

from tools import plotstyle as ps
from tools.fluids import fluid_props
from tools.pareto import pareto_indices
from tools.report import draw_requirements, plain_log_ticks
from tools.types import Requirements

from .slot_jet import SlotJetDesign, run_slot_jet

REQ = Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, 40.0, 5e3, "PG25", r_th_threshold_k_per_w=20e-3,
                   r_th_objective_k_per_w=15e-3)
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
    budget = REQ.dp_budget_pa / 1e3
    threshold, objective = REQ.r_th_threshold_k_per_w * 1e3, REQ.r_th_objective_k_per_w * 1e3
    xs, ys = [p[0] for p in fins + jets], [p[1] for p in fins + jets]
    x_lo, x_hi = 0.7 * min(xs), 1.5 * max(max(xs), budget)
    y_lo, y_hi = 0.7 * min(ys + [objective]), 1.3 * max(ys)

    fig, ax = ps.figure(5.6)
    draw_requirements(ax, x_lo, y_lo, budget, threshold, objective, y_hi)
    ax.plot([p[0] for p in fins], [p[1] for p in fins], "o-", color=ps.BLUE, ms=5.5, zorder=3)
    i_lbl = len(fins) // 2
    ax.annotate("Straight fins", (fins[i_lbl][0], fins[i_lbl][1]), xytext=(10, 4), textcoords="offset points", color=ps.BLUE,
                fontsize=ps.SMALL, ha="left", va="bottom")
    gap_colours = dict(zip(GAPS_MM, ("#7c2d12", "#ea580c", "#fbbf24")))
    slot_markers = dict(zip(SLOTS_MM, ("D", "s")))
    for gap, colour in gap_colours.items():
        for slot, marker in slot_markers.items():
            sel = [r for r in rows if abs(r["design"]["gap_m"] * 1e3 - gap) < 1e-9 and abs(r["design"]["slot_width_m"] * 1e3 - slot) < 1e-9]
            ax.scatter([r["result"]["dp_pa"] / 1e3 for r in sel], [r["result"]["r_th_k_per_w"] * 1e3 for r in sel], s=48, color=colour,
                       marker=marker, edgecolors=ps.INK, linewidths=0.6, zorder=4)
    front = pareto_indices([r["result"]["dp_pa"] for r in rows], [r["result"]["r_th_k_per_w"] for r in rows])
    ax.plot([rows[i]["result"]["dp_pa"] / 1e3 for i in front], [rows[i]["result"]["r_th_k_per_w"] * 1e3 for i in front],
            color=ps.ACCENT, lw=1.1, ls="--", zorder=2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    plain_log_ticks(ax, x_lo, x_hi, "x", ticks=(0.005, 0.01, 0.03, 0.1, 0.3, 1, 3, 5))
    plain_log_ticks(ax, y_lo, y_hi, "y")
    ax.set_xlabel("Pressure drop ΔP (kPa)")
    ax.set_ylabel("Thermal resistance R_th (mK/W)")
    ax.set_title("Slot-jet arrays vs straight fins")
    gaps = [Line2D([], [], ls="none", marker="o", ms=8, color=c, mec=ps.INK, mew=0.6, label=f"{g:.2f} mm")
            for g, c in gap_colours.items()]
    slots = [Line2D([], [], ls="none", marker=m, ms=7, color="white", mec=ps.INK, mew=1.0, label=f"{s:.1f} mm")
             for s, m in slot_markers.items()]
    front_h = [Line2D([], [], color=ps.ACCENT, lw=1.1, ls="--", label="Slot-jet front")]
    leg = ax.legend(handles=gaps, title="Jet gap (colour)", loc="lower left", bbox_to_anchor=(0.0, 0.0), alignment="left")
    ax.add_artist(leg)
    leg2 = ax.legend(handles=slots + front_h, title="Slot width (shape)", loc="lower left", bbox_to_anchor=(0.2, 0.0),
                     alignment="left")
    for lg in (leg, leg2):
        lg.get_title().set_fontsize(ps.SMALL + 0.5)
    ps.save(fig, out)
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
