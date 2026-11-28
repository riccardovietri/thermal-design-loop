"""Level 3 stage-two gates (issue #10): circular pins in 3D (OpenFOAM), checked on the committed results.

    uv run python validation/level3_pins.py      # seconds, no solver

Recomputes every gate from the run values saved in examples/output/level3_pins_vs_straight/results.json (not from the
recorded pass flags) and exits non-zero if one fails, or if the plan-view references there no longer match
examples/output/level2_planview/results.json. Gates of detached jobs that have not finished are reported as pending.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from level3.pins_study import GATES, OUT, PLANVIEW, run_gates  # noqa: E402


def check(results: Path = OUT / "results.json", planview: Path = PLANVIEW) -> int:
    data = json.loads(Path(results).read_text())
    pv = {r["design_id"]: r for r in json.loads(Path(planview).read_text())["designs"]}
    fails, lines = [], []

    def gate(name, value, limit, signed=True):
        ok = abs(value) < limit
        lines.append(f"{name}: {100 * value:+.3g}% (limit {100 * limit:g}%) {'pass' if ok else 'FAIL'}" if signed else
                     f"{name}: {value:.2g} (limit {limit:g}) {'pass' if ok else 'FAIL'}")
        if not ok:
            fails.append(name)

    # P0 against the plan view, and its mesh refinement
    p0 = data.get("p0")
    if not p0:
        fails.append("P0 missing")
    else:
        prod, fine = p0["runs"]["2"], p0["runs"]["4"]
        ref = pv[data["pick"]]["result"]
        for k in ("dpdx_pa_per_m", "h_w_m2k"):
            if abs(p0["planview"][k] / ref[k] - 1) > 1e-9:
                fails.append(f"STALE P0 plan-view {k}")
            gate(f"P0 {k} vs plan view", prod[k] / ref[k] - 1, GATES["p0"])
            gate(f"P0 {k} production vs 2x finer", prod[k] / fine[k] - 1, GATES["p0_mesh"])
    pc = data.get("parallel_check")
    if pc:
        for k in ("r_th_k_per_w", "dp_pa"):
            gate(f"parallel vs serial {k}", pc["parallel_np4"][k] / pc["serial"][k] - 1, GATES["parallel"])
    # per-run gates of the production pin cases
    for e in data["pins"]:
        if "level3" not in e:
            fails.append(f"{e['design_id']} not run")
            continue
        ref = pv[e["design_id"]]["result"]
        for k in ("r_th_k_per_w", "dp_pa"):
            if abs(e["planview"][k] / ref[k] - 1) > 1e-9:
                fails.append(f"STALE {e['design_id']} plan-view {k}")
        g = run_gates(e["level3"])
        gate(f"{e['design_id']} residual", g["residual"]["value"], GATES["residual"], signed=False)
        gate(f"{e['design_id']} R_th and ΔP flat over the last 10%", g["flat"]["value"], GATES["flat"])
        gate(f"{e['design_id']} energy balance", g["energy"]["value"], GATES["energy"])
    # long jobs
    lj = data.get("long_jobs", {})
    mg = lj.get("mesh_gate", {})
    pick = next(e for e in data["pins"] if e["is_pick"]).get("level3")
    if "result" in mg and pick:
        r = mg["result"]
        gate("mesh gate R_th production vs 2x finer", pick["r_th_k_per_w"] / r["r_th_k_per_w"] - 1, GATES["mesh"])
        gate("mesh gate ΔP production vs 2x finer", pick["dp_pa"] / r["dp_pa"] - 1, GATES["mesh"])
    else:
        lines.append(f"mesh gate: pending ({mg.get('state', 'not started')})")
    tr = lj.get("transient", {})
    lines.append("transient: " + ("finished, see summary" if "result" in tr else f"pending ({tr.get('state', 'not started')})"))
    print("\n".join(lines))
    if fails:
        print("FAILED: " + "; ".join(fails))
        return 1
    print("all recorded gates pass")
    return 0


if __name__ == "__main__":
    sys.exit(check())
