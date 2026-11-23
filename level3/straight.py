"""Level 3 stage one: 3D conjugate straight-fin unit cell (OpenFOAM, streamwise-periodic) against the Level 1 cross-section CFD.

    uv run python -m level3.straight                       # recommended + widest-channel design, three mesh levels each
    uv run python -m level3.straight --design w0387_t0387_h3750_q150 --levels 1 2
    uv run python -m level3.straight --summary-only        # rewrite the committed summary from runs/level3/

Inputs, logs, residual histories and results per run: runs/level3/<design_id>/level<k>/ (gitignored). Committed summary:
examples/output/level3_straight_validation/ (results.json, summary.md, mesh_convergence.png).
Needs Docker and the pinned OpenFOAM image (`level3.foam.IMAGE`). Gates G1-G4 are defined in docs/work-plan.md (issue #10).
"""
import argparse
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
from pathlib import Path

import numpy as np

from tools import plotstyle as ps
from tools.fluids import fluid_props
from tools.types import Design, Requirements

from .foam import IMAGE, parse_check_mesh, parse_dat, parse_solver_log, run_docker, write_case
from .geometry import MESH_LEVELS, UnitCell, duct_f_re_exact, grid, r_th_from_temperatures, unit_cell

LEVEL1 = Path("examples/output/gpu_coldplate/results.json")
RUNS = Path("runs/level3")
OUT = Path("examples/output/level3_straight_validation")
N_ITER = {0: 1000, 1: 1500, 2: 3000}  # SIMPLE iterations per mesh level (fixed count; see G4)
PRODUCTION, FINE = 1, 2
GATES = {"agreement": 0.02, "mesh": 0.02, "energy": 0.01, "residual": 1e-8, "flat": 0.001}


def requirements(level1: dict) -> Requirements:
    return Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, 40.0, level1["dp_budget_pa"], "PG25")


def reference(level1: dict, design_id: str) -> tuple[Design, dict]:
    d = next(x for x in level1["designs"] if x["cfd"]["design_id"] == design_id)
    return Design(**d["one_d"]["design"]), d["cfd"]


def widest(level1: dict) -> str:
    """The committed Level 1 design with the lowest fin-height to channel-width ratio (the widest channel)."""
    return min(level1["designs"], key=lambda d: d["one_d"]["design"]["h_fin_m"] / d["one_d"]["design"]["w_ch_m"])["cfd"]["design_id"]


def flatness(values: list[float], fraction: float = 0.1) -> float:
    """Largest relative departure from the final value over the last `fraction` of the iterations."""
    v = np.asarray(values, dtype=float)
    tail = v[-max(2, int(round(fraction * len(v)))):]
    return float(np.max(np.abs(tail / v[-1] - 1)))


def richardson(values: list[float]) -> dict:
    """Observed order and extrapolated value from three solutions on meshes refined by 2 in every direction."""
    f0, f1, f2 = values
    order = float(np.log((f0 - f1) / (f1 - f2)) / np.log(2))
    return {"order": order, "extrapolated": f2 + (f2 - f1) / (2**order - 1)}


def histories(case: Path, cell: UnitCell) -> dict:
    """Per-iteration R_th, ΔP and initial residuals from the run's logs and function-object files."""
    log = parse_solver_log((case / "log.solver").read_text())
    base = [r[1] for r in parse_dat(case / "postProcessing/solid/baseT/0/surfaceFieldValue.dat")]
    n = min(len(base), len(log["bulk_t"]))
    r_th = [r_th_from_temperatures(cell, base[i], log["bulk_t"][i]) for i in range(n)]
    grad = log["dpdx_pa_per_m"][:n]
    fluid = parse_dat(case / "postProcessing/fluid/residualsFluid/0/solverInfo.dat")
    solid = parse_dat(case / "postProcessing/solid/residualsSolid/0/solverInfo.dat")
    return {"r_th_k_per_w": r_th, "dp_pa": [cell.dp_core_pa(g) for g in grad], "flow_m3_s": log["flow_m3_s"][:n],
            "res_ux": [r[2] for r in fluid], "res_h_fluid": [r[13] for r in fluid], "res_p_rgh": [r[18] for r in fluid],
            "res_h_solid": [r[2] for r in solid]}


def _integral(path: Path, patch: str) -> float:
    """Last integral of wallHeatFlux on a patch (columns: time, patch, min, max, integral)."""
    rows = [ln.split() for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    return float([r for r in rows if r[1] == patch][-1][4])


def analyse(case: Path, cell: UnitCell, factor: int, seconds: float) -> dict:
    case = Path(case)
    h = histories(case, cell)
    g = grid(cell, factor)
    q_cell = cell.cell_heat_w(g["lx_m"])
    q_fluid = _integral(case / "postProcessing/fluid/heatFluid/0/wallHeatFlux.dat", "fluid_to_solid")
    q_bottom = _integral(case / "postProcessing/solid/heatSolid/0/wallHeatFlux.dat", "bottom")
    dpdx = h["dp_pa"][-1] / cell.flow_length_m
    f_re = dpdx * cell.hydraulic_diameter_m**2 / (2 * cell.mu * cell.u_bar)
    mesh = {r: parse_check_mesh((case / f"log.checkMesh.{r}").read_text()) for r in ("fluid", "solid")}
    peak = case / "memory.peak"
    return {
        "design_id": cell.design_id, "level": next(k for k, v in MESH_LEVELS.items() if v == factor), "grid": g,
        "r_th_k_per_w": h["r_th_k_per_w"][-1], "dp_pa": h["dp_pa"][-1], "f_re_fanning": f_re,
        "f_re_exact_duct": duct_f_re_exact(2 * cell.half_channel_m / cell.height_m),
        "reynolds": cell.reynolds, "u_bar_m_s": cell.u_bar,
        "flow_cell_m3_s": h["flow_m3_s"][-1], "flow_cell_target_m3_s": cell.cell_flow_m3_s(),
        "heat_applied_w": q_cell, "heat_into_fluid_w": q_fluid, "heat_recovered_at_base_w": q_bottom,
        "energy_balance_error": q_fluid / q_cell - 1, "base_flux_error": q_bottom / q_cell - 1,
        "iterations": len(h["r_th_k_per_w"]),
        "final_residuals": {"Ux": h["res_ux"][-1], "h_fluid": h["res_h_fluid"][-1], "h_solid": h["res_h_solid"][-1],
                            "p_rgh": h["res_p_rgh"][-1]},
        "flat_last_10pct": {"r_th": flatness(h["r_th_k_per_w"]), "dp": flatness(h["dp_pa"])},
        "mesh_quality": mesh, "seconds": seconds,
        "peak_memory_mb": int(peak.read_text().split()[0]) / 2**20 if peak.exists() else None,
    }


def run_case(cell: UnitCell, level: int, root: Path) -> dict:
    case = root / cell.design_id / f"level{level}"
    factor = MESH_LEVELS[level]
    write_case(case, cell, factor, N_ITER[level])
    t0 = time.perf_counter()
    run_docker(case)
    res = analyse(case, cell, factor, time.perf_counter() - t0)
    (case / "results.json").write_text(json.dumps(res, indent=2))
    (case / "history.json").write_text(json.dumps(histories(case, cell)))
    print(f"{cell.design_id} level {level}: R_th {res['r_th_k_per_w']*1e3:.4f} mK/W, ΔP {res['dp_pa']:.2f} Pa, "
          f"{res['grid']['n_cells_fluid'] + res['grid']['n_cells_solid']} cells, {res['seconds']:.0f} s", flush=True)
    return res


def gates(levels: dict, ref: dict) -> dict:
    prod, fine = levels[PRODUCTION], levels[FINE]
    d = lambda a, b: a / b - 1
    g1r, g1p = d(prod["r_th_k_per_w"], ref["r_th_k_per_w"]), d(prod["dp_pa"], ref["dp_pa"])
    g2r, g2p = d(prod["r_th_k_per_w"], fine["r_th_k_per_w"]), d(prod["dp_pa"], fine["dp_pa"])
    worst_energy = max(abs(x["energy_balance_error"]) for x in (prod, fine))  # base flux is imposed: not a check
    worst_res = max(v for x in (prod, fine) for k, v in x["final_residuals"].items() if k != "p_rgh")
    worst_flat = max(v for x in (prod, fine) for v in x["flat_last_10pct"].values())
    pct = lambda v: f"{100 * v:+.2f}%"
    return {
        "G1 R_th vs Level 1": {"value": pct(g1r), "limit": "±2%", "pass": abs(g1r) < GATES["agreement"]},
        "G1 ΔP vs Level 1": {"value": pct(g1p), "limit": "±2%", "pass": abs(g1p) < GATES["agreement"]},
        "G2 R_th production vs fine": {"value": pct(g2r), "limit": "±2%", "pass": abs(g2r) < GATES["mesh"]},
        "G2 ΔP production vs fine": {"value": pct(g2p), "limit": "±2%", "pass": abs(g2p) < GATES["mesh"]},
        "G3 energy balance: heat into the coolant across the interface vs applied": {"value": f"{100 * worst_energy:.1e}%", "limit": "1%",
                                                          "pass": worst_energy < GATES["energy"]},
        "G4 final initial residual (Ux, h fluid, h solid)": {"value": f"{worst_res:.1e}", "limit": "< 1e-8",
                                                            "pass": worst_res < GATES["residual"]},
        "G4 R_th and ΔP flat over the last 10% of iterations": {"value": f"{100 * worst_flat:.1e}%", "limit": "0.1%",
                                                               "pass": worst_flat < GATES["flat"]},
    }


def collect(level1: dict, ids: list[str], root: Path = RUNS) -> dict:
    out = {"image": IMAGE, "solver": "chtMultiRegionSimpleFoam (OpenFOAM v2406)", "level1": str(LEVEL1),
           "recommended": level1["recommended"], "n_iter": N_ITER, "mesh_levels": MESH_LEVELS, "gate_limits": GATES, "designs": {}}
    for did in ids:
        design, ref = reference(level1, did)
        levels = {}
        for lv in MESH_LEVELS:
            p = root / did / f"level{lv}" / "results.json"
            if p.exists():
                levels[lv] = json.loads(p.read_text())
        entry = {"design": asdict(design), "level1_cfd": {k: ref[k] for k in ("r_th_k_per_w", "dp_pa", "f_re_fanning",
                                                                              "reynolds", "n_channels")},
                 "levels": {str(k): v for k, v in levels.items()}}
        if PRODUCTION in levels and FINE in levels:
            entry["gates"] = gates(levels, ref)
        if len(levels) == 3:
            entry["richardson"] = {k: richardson([levels[lv][k] for lv in sorted(levels)]) for k in ("r_th_k_per_w", "dp_pa")}
        out["designs"][did] = entry
    return out


def plot(data: dict, out: Path) -> Path:
    fig, axes = ps.figure(3.6, 1, 2)
    colours = [ps.ACCENT, ps.BLUE]
    for (did, e), col in zip(data["designs"].items(), colours):
        lv = sorted(e["levels"].values(), key=lambda r: r["level"])
        cells = [r["grid"]["n_cells_fluid"] + r["grid"]["n_cells_solid"] for r in lv]
        g = e["design"]
        label = f"{g['w_ch_m']*1e3:.2f} / {g['t_fin_m']*1e3:.2f} / {g['h_fin_m']*1e3:.2f} mm"
        for ax, key in zip(axes, ("r_th_k_per_w", "dp_pa")):
            ax.plot(cells, [100 * (r[key] / e["level1_cfd"][key] - 1) for r in lv], "o-", color=col, label=label)
    for ax, name in zip(axes, ("R_th", "ΔP")):
        ax.axhspan(-2, 2, color=ps.FAINT, alpha=0.5, lw=0)
        ax.axhline(0, color=ps.MUTED, lw=0.8)
        ax.set_xscale("log")
        ax.set_xlabel("cells (fluid + solid)")
        ax.set_ylabel(f"{name}, 3D − Level 1 (%)")
        ax.set_ylim(-3, 3)
    axes[0].legend(loc="lower right")
    fig.suptitle("3D periodic cell against Level 1 cross-section CFD")
    ps.save(fig, out)
    return out


def write_summary(data: dict, out: Path) -> Path:
    rec = data["recommended"]
    a = []
    a.append("# Level 3 stage one: 3D straight-fin validation\n")
    a.append(f"OpenFOAM v2406 `chtMultiRegionSimpleFoam` in Docker (`{data['image']}`), steady laminar, constant properties. "
             "One half fin and one half channel with symmetry planes on both sides, 2 mm copper base under uniform flux, "
             "adiabatic cover, no tip clearance. Streamwise-periodic (fully developed): mean velocity held at the Level 1 value, "
             "temperature solved as a linear rise plus a periodic part. R_th = (base-underside area average − local bulk) / Q "
             "plus half the coolant rise, ΔP = dp/dx × 60 mm, as Level 1. Reference values are read from "
             f"`{data['level1']}`.\n")
    for did, e in data["designs"].items():
        ref, g = e["level1_cfd"], e["design"]
        tag = " (Level 1 recommended)" if did == rec else " (widest channel, lowest H / w)"
        a.append(f"## `{did}`{tag}\n")
        a.append(f"{g['w_ch_m']*1e3:.3f} mm channel, {g['t_fin_m']*1e3:.3f} mm fin, {g['h_fin_m']*1e3:.2f} mm tall, "
                 f"{ref['n_channels']} channels, Re {ref['reynolds']:.0f}. Level 1: {ref['r_th_k_per_w']*1e3:.3f} mK/W, "
                 f"{ref['dp_pa']:.1f} Pa.\n")
        a.append("| Level | Cells (fluid + solid) | R_th (mK/W) | vs Level 1 | ΔP (Pa) | vs Level 1 | f·Re (exact duct) | "
                 "Energy balance | Iterations | Final residual | Run time (s) | Peak memory (MB) |")
        a.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for k in sorted(e["levels"], key=int):
            r = e["levels"][k]
            cells = r["grid"]["n_cells_fluid"] + r["grid"]["n_cells_solid"]
            mark = " (production)" if int(k) == PRODUCTION else (" (fine)" if int(k) == FINE else "")
            res = max(v for kk, v in r["final_residuals"].items() if kk != "p_rgh")
            a.append(f"| {k}{mark} | {cells:,} | {r['r_th_k_per_w']*1e3:.4f} | {100*(r['r_th_k_per_w']/ref['r_th_k_per_w']-1):+.2f}% | "
                     f"{r['dp_pa']:.2f} | {100*(r['dp_pa']/ref['dp_pa']-1):+.2f}% | {r['f_re_fanning']:.3f} ({r['f_re_exact_duct']:.3f}) | "
                     f"{100*r['energy_balance_error']:+.1e}% | {r['iterations']} | {res:.1e} | {r['seconds']:.0f} | "
                     f"{r['peak_memory_mb']:.0f} |")
        if "richardson" in e:
            rr, rp = e["richardson"]["r_th_k_per_w"], e["richardson"]["dp_pa"]
            a.append(f"\nRichardson extrapolation over the three levels (observed order {rr['order']:.2f} on R_th, {rp['order']:.2f} "
                     f"on ΔP): R_th {rr['extrapolated']*1e3:.4f} mK/W ({100*(rr['extrapolated']/ref['r_th_k_per_w']-1):+.2f}% from "
                     f"Level 1), ΔP {rp['extrapolated']:.2f} Pa ({100*(rp['extrapolated']/ref['dp_pa']-1):+.2f}%). Not used by any gate.")
        if "gates" in e:
            a.append("\n| Gate | Value | Limit | Result |\n|---|---|---|---|")
            a += [f"| {k} | {v['value']} | {v['limit']} | {'pass' if v['pass'] else 'FAIL'} |" for k, v in e["gates"].items()]
        a.append("")
    a.append("![Mesh convergence](mesh_convergence.png)\n")
    a.append("R_th (left) and ΔP (right) of the 3D periodic cell relative to the Level 1 cross-section CFD at three mesh levels "
             "(each 2× finer in every direction); the band is ±2% (gate G1). Level 1 is the production mesh of the committed "
             "study; the fine level is the refinement check (gate G2).\n")
    a.append("Notes. The energy balance compares the heat that crosses the fluid-solid interface, computed from the fluid-side "
             "temperature gradient, with the heat applied at the base (the base flux itself is a boundary condition, so it is not a check); "
             "in the periodic formulation the coolant's enthalpy rise is the sink term, which removes exactly the applied heat, "
             "so G3 is a conservation check of the coupled solve rather than an independent caloric measurement. The f·Re "
             "column compares the 3D flow solve with the exact series solution for a rectangular duct (Shah & London 1978), "
             "independent of Level 1. ΔP is a pressure gradient held by `meanVelocityForce` times 60 mm; normalised residuals "
             "of the cross-stream velocity and pressure are round-off (|u_y|, |u_z| ≈ 1e-12 m/s) and are not used. "
             "Only straight fins are run here; pins are stage two.\n")
    out.write_text("\n".join(a))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--design", action="append", help="Level 1 design id (default: recommended and widest channel)")
    ap.add_argument("--levels", type=int, nargs="+", default=list(MESH_LEVELS))
    ap.add_argument("--level1", type=Path, default=LEVEL1)
    ap.add_argument("--summary-only", action="store_true", help="rewrite the committed summary from runs/level3/ (no solves)")
    ap.add_argument("--jobs", type=int, default=6, help="cases run at once, one CPU each")
    a = ap.parse_args(argv)
    level1 = json.loads(a.level1.read_text())
    ids = a.design or [level1["recommended"], widest(level1)]
    if not a.summary_only:
        req = requirements(level1)
        fluid = fluid_props(req.coolant, req.inlet_temp_c)
        cells = {did: unit_cell(reference(level1, did)[0], req, fluid) for did in ids}
        with ThreadPoolExecutor(a.jobs) as pool:
            list(pool.map(lambda job: run_case(cells[job[0]], job[1], RUNS), [(d, lv) for d in ids for lv in a.levels]))
    data = collect(level1, ids)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "results.json").write_text(json.dumps(data, indent=2))
    plot(data, OUT / "mesh_convergence.png")
    write_summary(data, OUT / "summary.md")
    failed = [f"{did}: {k}" for did, e in data["designs"].items() if did == level1["recommended"]
              for k, v in e.get("gates", {}).items() if not v["pass"]]
    for did, e in data["designs"].items():
        for k, v in e.get("gates", {}).items():
            print(f"{did} {k}: {v['value']} (limit {v['limit']}) {'pass' if v['pass'] else 'FAIL'}")
    if failed:
        print("GATE FAILED: " + "; ".join(failed), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
