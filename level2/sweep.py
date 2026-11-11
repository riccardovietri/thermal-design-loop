"""Level 2 sweep: slot-jet arrays over slot width, gap and pitch, compared with the Level 1 straight-fin front."""
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
from tools.types import Requirements

from .slot_jet import SlotJetDesign, run_slot_jet

REQ = Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, 40.0, 30e3, "PG25")
SLOTS_MM, GAPS_MM, PITCHES_MM, FLOW_LPM = (0.15, 0.3), (0.05, 0.1, 0.2, 0.4), (2.0, 4.0, 8.0), 1.5


def main() -> int:
    out = Path("examples/output/level2_slot_jet")
    out.mkdir(parents=True, exist_ok=True)
    fluid = fluid_props(REQ.coolant, REQ.inlet_temp_c)
    rows = []
    for slot, gap, pitch in itertools.product(SLOTS_MM, GAPS_MM, PITCHES_MM):
        d = SlotJetDesign(slot / 1e3, gap / 1e3, pitch / 1e3, FLOW_LPM)
        t0 = time.perf_counter()
        try:
            r = run_slot_jet(d, REQ, fluid, Path("runs/level2") / d.design_id)
        except RuntimeError as err:
            print(f"{d.design_id}: FAILED {err}", flush=True)
            continue
        rows.append({"design": asdict(d), "result": asdict(r), "seconds": time.perf_counter() - t0})
        print(f"{d.design_id}: Re={r.reynolds_slot:5.1f} dP={r.dp_pa/1e3:8.3f} kPa R={r.r_th_k_per_w*1e3:7.2f} mK/W "
              f"energy err={100*r.energy_balance_error:+.3f}% {rows[-1]['seconds']:.0f}s", flush=True)
        (out / "results.json").write_text(json.dumps(rows, indent=2))

    level1 = json.loads(Path("examples/output/gpu_coldplate/results.json").read_text())["designs"]
    fig, ax = plt.subplots(figsize=(7.8, 5.4), constrained_layout=True)
    ax.plot([d["cfd"]["dp_pa"] / 1e3 for d in level1], [d["cfd"]["r_th_k_per_w"] * 1e3 for d in level1], "o-", color="#1d4ed8",
            label="Straight fins (Level 1 CFD front, 1–2 LPM)")
    for gap, colour in zip(GAPS_MM, ("#7c2d12", "#c2410c", "#f59e0b", "#fcd34d")):
        sel = [r for r in rows if abs(r["design"]["gap_m"] * 1e3 - gap) < 1e-9]
        ax.scatter([r["result"]["dp_pa"] / 1e3 for r in sel], [r["result"]["r_th_k_per_w"] * 1e3 for r in sel], s=46, color=colour,
                   marker="s" if gap >= 0.15 else "D", edgecolors="#374151", linewidths=0.6, label=f"Slot jets, {gap:.2f} mm gap"
                   + ("" if gap >= 0.15 else " (below 0.15 mm feature limit)"))
    front = pareto_indices([r["result"]["dp_pa"] for r in rows], [r["result"]["r_th_k_per_w"] for r in rows])
    ax.plot([rows[i]["result"]["dp_pa"] / 1e3 for i in front], [rows[i]["result"]["r_th_k_per_w"] * 1e3 for i in front],
            color="#c2410c", lw=1, ls="--", label="Slot-jet Pareto front")
    ax.axvline(REQ.dp_budget_pa / 1e3, color="#9ca3af", ls=":", lw=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Pressure drop ΔP (kPa)")
    ax.set_ylabel("Thermal resistance R_th (mK/W)")
    ax.set_title(f"Slot-jet arrays on a flat base vs straight fins ({FLOW_LPM} LPM PG25, 1 kW, 60 × 60 mm)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(out / "pareto_compare.png", dpi=200)
    print(out / "pareto_compare.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
