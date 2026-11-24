"""Level 3 stage two: staggered circular pins against straight fins in 3D (OpenFOAM, issue #10).

    uv run python -m level3.pins_study --p0               # P0: pick as a slip slab, fluid only, against the plan view
    uv run python -m level3.pins_study --parallel-check   # coarse conjugate pick, serial against 4-way parallel
    uv run python -m level3.pins_study --steady           # three pin designs + plan-view straight pick (production)
    uv run python -m level3.pins_study --summary-only     # rewrite results.json, summary and figures from runs/level3/
    uv run python -m level3.pins_study --launch-long      # start the detached mesh-gate and transient jobs
    uv run python -m level3.status                        # state and latest monitored values of the detached jobs

Runs: runs/level3/pins/<design_id>/<kind>_f<factor>/ (gitignored); straight fins reuse level3.straight under
runs/level3/<design_id>/level<k>/. Committed output: examples/output/level3_pins_vs_straight/.
"""
import argparse
import json
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from tools import plotstyle as ps
from tools.fluids import fluid_props
from tools.pareto import pareto_indices
from tools.types import Design

from . import pins as lp
from .foam import IMAGE, parse_check_mesh, parse_dat, parse_solver_log
from .straight import RUNS as STRAIGHT_RUNS, flatness, requirements

PLANVIEW = Path("examples/output/level2_planview/results.json")
LEVEL1 = Path("examples/output/gpu_coldplate/results.json")
STRAIGHT3D = Path("examples/output/level3_straight_validation/results.json")
RUNS = Path("runs/level3/pins")
LONG = Path("runs/level3/long")
OUT = Path("examples/output/level3_pins_vs_straight")
PRODUCTION, FINE, COARSE = 2, 4, 1  # plan-mesh factors; z layers scale with them (lp.Z_CELLS)
P0_FACTORS = (PRODUCTION, FINE)
PROD_PROCS = 3  # processes per production pin case (three cases at once on the 10-CPU Docker VM)
N_ITER = {COARSE: 3000, PRODUCTION: 5000, FINE: 10000, "p0": 4000}
GATES = {"p0": 0.03, "p0_mesh": 0.02, "mesh": 0.02, "energy": 0.01, "residual": 1e-8, "flat": 0.001, "parallel": 1e-3}
STRAIGHT_PICK = Design(0.4e-3, 0.4e-3, 4.0e-3, 1.5)  # plan-view straight pick, straight_w0400_t0400_h_fin4000


def setup():
    level1 = json.loads(LEVEL1.read_text())
    req = requirements(level1)
    return level1, req, fluid_props(req.coolant, req.inlet_temp_c)


def pin_front(pv: dict) -> list[dict]:
    """Eligible circular-pin plan-view designs on the topology's Pareto front, by ΔP (as level2.planview_sweep)."""
    from level2.planview_sweep import aspect_exclusion
    rows = [r for r in pv["designs"] if r["topology"] == "pin_circle" and r["status"] == "converged" and not aspect_exclusion(r)]
    idx = pareto_indices([r["result"]["dp_pa"] for r in rows], [r["result"]["r_th_k_per_w"] for r in rows])
    return sorted((rows[i] for i in idx), key=lambda r: r["result"]["dp_pa"])


def pin_designs(pv: dict) -> list[dict]:
    """The cross-topology pick and its two neighbours on the eligible circular-pin front."""
    pick = pv["summary"]["recommended"]["design_id"]
    front = pin_front(pv)
    i = [r["design_id"] for r in front].index(pick)
    return front[max(0, i - 1): i + 2]


def straight_planview_pick(pv: dict) -> dict:
    return next(r for r in pv["designs"] if r["design_id"] == pv["summary"]["per_topology"]["straight"]["recommended"]["design_id"])


# --------------------------------------------------------------------------------------------------------------------
# case set-up and running


def case_dir(design_id: str, kind: str, factor: int, n_procs: int = 1, root: Path = RUNS) -> Path:
    return root / design_id / (f"{kind}_f{factor}" + (f"_np{n_procs}" if n_procs > 1 else ""))


def prepare(params: dict, kind: str, factor: int, n_iter: int, n_procs: int = 1, full_pitch: bool = False,
            case: Path | None = None) -> Path:
    _, req, fluid = setup()
    cell = lp.pin_cell(params, req, fluid, full_pitch=full_pitch)
    case = Path(case or case_dir(cell.design_id, kind, factor, n_procs))
    shutil.rmtree(case, ignore_errors=True)
    mode = "p0" if kind == "p0" else "conjugate"
    mesh = lp.write_mesh(case / "mesh.msh", params, factor, mode, full_pitch)
    info = lp.write_pin_case(case, cell, mesh, mode, n_iter, n_procs=n_procs)
    (case / "Allrun").write_text(lp.allrun(info["regions"], n_procs))
    (case / "setup.json").write_text(json.dumps({"params_mm": params, "kind": kind, "factor": factor, "n_iter": n_iter,
                                                 "n_procs": n_procs, "full_pitch": full_pitch, "mesh": mesh, "case": info},
                                                indent=2))
    return case


def finish(case: Path) -> dict:
    """Analyse a case whose container ran to the end without the Python driver (results.json from the logs; the run time
    is the solver's last ClockTime)."""
    import re
    case = Path(case)
    res = analyse(case)
    res["seconds"] = float(re.findall(r"ClockTime = (\S+) s", (case / "log.solver").read_text())[-1])
    (case / "results.json").write_text(json.dumps(res, indent=2))
    return res


def run(case: Path, cpus: int = 1) -> dict:
    t0 = time.perf_counter()
    lp.run_docker(case, cpus)
    s = json.loads((case / "setup.json").read_text())
    res = analyse(case)
    res["seconds"] = time.perf_counter() - t0
    (case / "results.json").write_text(json.dumps(res, indent=2))
    print(f"{case}: " + ", ".join(f"{k} {res[k]:.6g}" for k in ("r_th_k_per_w", "dp_pa", "h_w_m2k", "dpdx_pa_per_m")
                                    if res.get(k) is not None) + f", {res['seconds']:.0f} s", flush=True)
    return res


# --------------------------------------------------------------------------------------------------------------------
# analysis


def _cell(s: dict) -> lp.PinCell:
    _, req, fluid = setup()
    return lp.pin_cell(s["params_mm"], req, fluid, full_pitch=s["full_pitch"])


def residual_histories(case: Path) -> dict:
    out = {}
    for region in ("fluid", "solid"):
        f = case / f"postProcessing/{region}/residuals{region.capitalize()}/0/solverInfo.dat"
        if not f.exists():
            continue
        lines = f.read_text().splitlines()
        head = lines[1].lstrip("#").split()
        rows = [ln.split() for ln in lines[2:] if ln.strip()]
        for name in head:
            if name.endswith("_initial"):
                j = head.index(name)
                out[f"{name[:-8]}_{region}"] = [float(r[j]) for r in rows]
    return out


def _wall_integral(path: Path, patch: str) -> float:
    rows = [ln.split() for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    return float([r for r in rows if r[1] == patch][-1][4])


def histories(case: Path, cell: lp.PinCell, kind: str) -> dict:
    log = parse_solver_log((case / "log.solver").read_text())
    if kind == "p0":
        wall = [r[1] for r in parse_dat(case / "postProcessing/fluid/wallT/0/surfaceFieldValue.dat")]
        s = json.loads((case / "setup.json").read_text())
        q = s["case"]["q_wall_w_m2"]
        n = min(len(wall), len(log["bulk_t"]))
        h = [lp.planview_h(q, wall[i], log["bulk_t"][i]) for i in range(n)]
        return {"h_w_m2k": h, "dpdx_pa_per_m": log["dpdx_pa_per_m"][:n], "flow_m3_s": log["flow_m3_s"][:n]}
    base = [r[1] for r in parse_dat(case / "postProcessing/solid/baseT/0/surfaceFieldValue.dat")]
    n = min(len(base), len(log["bulk_t"]))
    return {"r_th_k_per_w": [lp.r_th_pins(cell, base[i], log["bulk_t"][i]) for i in range(n)],
            "dp_pa": [cell.dp_core_pa(g) for g in log["dpdx_pa_per_m"][:n]], "flow_m3_s": log["flow_m3_s"][:n]}


def analyse(case: Path) -> dict:
    case = Path(case)
    s = json.loads((case / "setup.json").read_text())
    kind, cell = s["kind"], _cell(s)
    h = histories(case, cell, kind)
    res_h = residual_histories(case)
    # cross-stream velocity in the P0 slab is identically zero: its normalised residual is round-off, not used
    used = [k for k in res_h if not (kind == "p0" and k.startswith("Uz"))]
    final_res = {k: res_h[k][-1] for k in res_h}
    regions = s["case"]["regions"]
    mesh_q = {r: parse_check_mesh((case / f"log.checkMesh.{r}").read_text()) for r in regions
              if (case / f"log.checkMesh.{r}").exists()}
    peak = case / "memory.peak"
    out = {"design_id": cell.design_id, "kind": kind, "factor": s["factor"], "n_procs": s["n_procs"],
           "params_mm": s["params_mm"], "n_cells": s["mesh"]["n_cells"], "plan": {k: s["mesh"][k] for k in
           ("fluid_area_mm2", "pin_area_mm2", "pin_perimeter_mm", "n_plan_cells", "wall_mm", "core_mm")},
           "ubar_m_s": s["case"]["ubar_m_s"], "u_superficial_m_s": cell.u_superficial, "porosity": cell.porosity,
           "re_d_max": cell.re_d_max, "iterations": len(next(iter(h.values()))),
           "flow_cell_m3_s": h["flow_m3_s"][-1], "flow_cell_target_m3_s": (
               cell.u_superficial * lp.P0_SLAB_MM * 1e-3 * cell.lateral_m if kind == "p0" else cell.cell_flow_m3_s),
           "final_residuals": final_res, "max_final_residual_used": max(final_res[k] for k in used),
           "residuals_used": used, "mesh_quality": mesh_q,
           "peak_memory_mb": int(peak.read_text().split()[0]) / 2**20 if peak.exists() and peak.read_text().strip() else None}
    q_cell = s["case"]["q_cell_w"]
    if kind == "p0":
        q_fluid = _wall_integral(case / "postProcessing/fluid/heatFluid/0/wallHeatFlux.dat", "pins")
        out.update({"h_w_m2k": h["h_w_m2k"][-1], "dpdx_pa_per_m": h["dpdx_pa_per_m"][-1],
                    "flat_last_10pct": {"h": flatness(h["h_w_m2k"]), "dpdx": flatness(h["dpdx_pa_per_m"])},
                    "heat_applied_w": q_cell, "heat_into_fluid_w": q_fluid,
                    # the wall flux is imposed; the check is that the bulk sink and the wall input agree
                    "energy_balance_error": q_fluid / q_cell - 1})
    else:
        q_fluid = _wall_integral(case / "postProcessing/fluid/heatFluid/0/wallHeatFlux.dat", "fluid_to_solid")
        q_bottom = _wall_integral(case / "postProcessing/solid/heatSolid/0/wallHeatFlux.dat", "bottom")
        out.update({"r_th_k_per_w": h["r_th_k_per_w"][-1], "dp_pa": h["dp_pa"][-1],
                    "dpdx_pa_per_m": h["dp_pa"][-1] / cell.flow_length_m,
                    "flat_last_10pct": {"r_th": flatness(h["r_th_k_per_w"]), "dp": flatness(h["dp_pa"])},
                    "heat_applied_w": q_cell, "heat_into_fluid_w": q_fluid, "heat_recovered_at_base_w": q_bottom,
                    # heat entering the coolant across the fluid-solid interface (fluid-side gradient) vs applied
                    "energy_balance_error": -q_fluid / q_cell - 1 if q_fluid < 0 else q_fluid / q_cell - 1})
    return out


def run_gates(r: dict) -> dict:
    flat = max(r["flat_last_10pct"].values())
    return {"residual": {"value": r["max_final_residual_used"], "limit": GATES["residual"],
                         "pass": r["max_final_residual_used"] < GATES["residual"]},
            "flat": {"value": flat, "limit": GATES["flat"], "pass": flat < GATES["flat"]},
            "energy": {"value": abs(r["energy_balance_error"]), "limit": GATES["energy"],
                       "pass": abs(r["energy_balance_error"]) < GATES["energy"]}}


# --------------------------------------------------------------------------------------------------------------------
# steps


def step_p0(pv: dict, jobs: int = 2) -> None:
    pick = next(r for r in pv["designs"] if r["design_id"] == pv["summary"]["recommended"]["design_id"])
    cases = [prepare(pick["params_mm"], "p0", f, N_ITER["p0"]) for f in P0_FACTORS]
    with ThreadPoolExecutor(jobs) as pool:
        list(pool.map(run, cases))


def step_parallel_check(pv: dict) -> None:
    pick = next(r for r in pv["designs"] if r["design_id"] == pv["summary"]["recommended"]["design_id"])
    cases = [(prepare(pick["params_mm"], "conj", COARSE, N_ITER[COARSE]), 1),
             (prepare(pick["params_mm"], "conj", COARSE, N_ITER[COARSE], n_procs=4), 4)]
    with ThreadPoolExecutor(2) as pool:
        list(pool.map(lambda c: run(*c), cases))


def step_steady(pv: dict, n_procs: int = None) -> None:
    n_procs = n_procs or PROD_PROCS
    from .straight import run_case
    from .geometry import unit_cell
    _, req, fluid = setup()
    cases = [prepare(r["params_mm"], "conj", PRODUCTION, N_ITER[PRODUCTION], n_procs=n_procs) for r in pin_designs(pv)]
    cell = unit_cell(STRAIGHT_PICK, req, fluid)
    with ThreadPoolExecutor(len(cases) + 2) as pool:
        futs = [pool.submit(run, c, n_procs) for c in cases]
        futs += [pool.submit(run_case, cell, lv, STRAIGHT_RUNS) for lv in (1, 2)
                 if not (STRAIGHT_RUNS / cell.design_id / f"level{lv}" / "results.json").exists()]
        for f in futs:
            f.result()


# --------------------------------------------------------------------------------------------------------------------
# long detached jobs


def long_jobs(pv: dict) -> dict:
    pick = next(r for r in pv["designs"] if r["design_id"] == pv["summary"]["recommended"]["design_id"])
    return {"mesh_gate": {"params_mm": pick["params_mm"], "factor": FINE, "n_procs": 5,
                          "source": str(case_dir(pick["design_id"], "conj", PRODUCTION, PROD_PROCS)),
                          "what": "pick at 2x the production mesh in every direction (steady, started from the production solution)"},
            "transient": {"params_mm": pick["params_mm"], "factor": PRODUCTION, "n_procs": 5,
                          "what": "pick on the full lateral pitch (periodic sides), steady start, perturbed, then transient"}}


# about 15 shedding periods if the wake sheds at St 0.2 on d and the gap velocity, and 10 passes of the flow through the period
TRANSIENT = {"steady_iter": 3000, "end_time": 0.15, "max_co": 1.5, "write_interval": 0.01, "dt0": 1e-6}
MESH_GATE_ITER = 5000  # started from the production solution mapped onto the finer mesh


def _transient_allrun(regions, n_procs, steady_iter) -> str:
    app_s, app_t = "chtMultiRegionSimpleFoam", "chtMultiRegionFoam"
    par = lambda app, log: (f"mpirun --oversubscribe -np {n_procs} {app} -parallel > {log} 2>&1 || true")
    lines = ["set -e", "cd /case", "gmshToFoam mesh.msh > log.gmshToFoam 2>&1",
             "transformPoints -scale '(0.001 0.001 0.001)' > log.transformPoints 2>&1",
             "splitMeshRegions -cellZonesOnly -overwrite > log.splitMeshRegions 2>&1", "rm -rf constant/polyMesh"]
    for r in regions:
        lines += [f"createPatch -region {r} -overwrite > log.createPatch.{r} 2>&1", f"checkMesh -region {r} > log.checkMesh.{r} 2>&1"]
    lines += ["rm -rf 0 && cp -r 0.orig 0", "rm -rf processor*", "decomposePar -allRegions > log.decomposePar.steady 2>&1",
              par(app_s, "log.solver.steady"), "reconstructPar -allRegions -latestTime > log.reconstructPar.steady 2>&1",
              f"rm -rf 0 processor* && mv {steady_iter} 0 && rm -rf 0/uniform 0/*/uniform && mv postProcessing postProcessing.steady",
              "mv system system.steady && cp -r system.transient system",
              "setFields -region fluid > log.setFields 2>&1",
              "decomposePar -allRegions > log.decomposePar 2>&1", "touch transient.started",
              par(app_t, "log.solver"), "reconstructPar -allRegions -latestTime > log.reconstructPar 2>&1 || true",
              "cat /sys/fs/cgroup/memory.peak > memory.peak 2>/dev/null || true"]
    return "\n".join(lines) + "\n"


def prepare_long(name: str, spec: dict, case: Path | None = None, transient: dict | None = None) -> Path:
    """Case directory of a detached job; `case` and `transient` override the defaults (used for a short smoke test)."""
    case = Path(case or LONG / name)
    transient = transient or TRANSIENT
    if name == "mesh_gate":
        prepare(spec["params_mm"], "conj", spec["factor"], MESH_GATE_ITER, n_procs=spec["n_procs"], case=case)
        # start from the converged production solution (mapped, same geometry); fall back to the uniform start
        mp = ("mapFields /src -consistent -sourceTime latestTime -sourceRegion {r} -targetRegion {r} > log.mapFields.{r} 2>&1")
        allrun = (case / "Allrun").read_text().replace(
            "rm -rf 0 && cp -r 0.orig 0\n",
            "rm -rf 0 && cp -r 0.orig 0\n(" + " && ".join(mp.format(r=r) for r in ("fluid", "solid"))
            + " && touch mapped) || (rm -rf 0 && cp -r 0.orig 0)\n")
        (case / "Allrun").write_text(allrun)
        return case
    p = spec["params_mm"]
    n_steady = transient["steady_iter"]
    prepare(p, "conj", spec["factor"], n_steady, n_procs=spec["n_procs"], full_pitch=True, case=case)
    tmp = case.parent / (case.name + "_tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    s = json.loads((case / "setup.json").read_text())
    _, req, fluid = setup()
    cell = lp.pin_cell(p, req, fluid, full_pitch=True)
    lp.write_pin_case(tmp, cell, s["mesh"], "conjugate", 0, transient=transient, n_procs=spec["n_procs"])
    shutil.rmtree(case / "system.transient", ignore_errors=True)
    shutil.copytree(tmp / "system", case / "system.transient")
    shutil.rmtree(tmp)
    # perturbation: a lateral velocity of 20% of the mean in a box in the wake of pin A, on one side only
    d, sl = p["d"] * 1e-3, p["s_l"] * 1e-3
    ub = s["case"]["ubar_m_s"]
    box = f"({d / 2:.6e} 0 0) ({d / 2 + 0.3 * sl:.6e} {d / 2:.6e} 1)"
    (case / "system.transient/fluid/setFieldsDict").write_text(lp._header("setFieldsDict") + f"""
defaultFieldValues ();
regions
(
    boxToCell {{ box {box}; fieldValues ( volVectorFieldValue U ({ub!r} {0.2 * ub!r} 0) ); }}
);
""")
    # wake probes on the pin-A centreline (y = 0) at mid-height: u_y there is zero in a steady symmetric wake
    zmid = cell.base_m + cell.height_m / 2
    pts = " ".join(f"({d / 2 + f * (sl - d):.6e} 0 {zmid:.6e})" for f in (0.25, 0.5, 0.75))
    cd = (case / "system.transient/controlDict").read_text()
    cd = cd.replace("functions\n{", "functions\n{\n    wake\n    {\n        type probes; libs (sampling); region fluid; "
                    f"fields (U T); probeLocations ({pts}); writeControl timeStep; writeInterval 1;\n    }}\n    baseTt\n    {{\n"
                    "        type surfaceFieldValue; libs (fieldFunctionObjects); region solid; regionType patch; name bottom;\n"
                    "        operation areaAverage; fields (T); writeFields false; writeControl timeStep; writeInterval 1; log false;\n    }", 1)
    (case / "system.transient/controlDict").write_text(cd)
    (case / "Allrun").write_text(_transient_allrun(s["case"]["regions"], spec["n_procs"], n_steady))
    s["transient"] = transient
    (case / "setup.json").write_text(json.dumps(s, indent=2))
    return case


def launch_long(pv: dict) -> None:
    LONG.mkdir(parents=True, exist_ok=True)
    jobs = {}
    for name, spec in long_jobs(pv).items():
        case = LONG / name
        if not ((case / "Allrun").exists() and (case / "mesh.msh").exists()) or (case / "log.solver").exists() \
                or (case / "log.solver.steady").exists():
            case = prepare_long(name, spec)  # fresh case unless one was prepared and never run
        cname = f"l3pins_{name}"
        subprocess.run(["docker", "rm", "-f", cname], capture_output=True)
        cid = subprocess.run(lp.docker_cmd(case, spec["n_procs"], detach=True, name=cname, source=spec.get("source")),
                             check=True, capture_output=True, text=True).stdout.strip()
        # keep the Mac awake while the container runs; survives this session
        with open(case / "caffeinate.log", "w") as fh:
            subprocess.Popen(["nohup", "caffeinate", "-i", "docker", "wait", cname], stdout=fh, stderr=subprocess.STDOUT,
                             start_new_session=True)
        jobs[name] = {**spec, "case": str(case), "container": cname, "container_id": cid,
                      "started": time.strftime("%Y-%m-%d %H:%M:%S")}
        print(f"started {name}: container {cname}, case {case}", flush=True)
    (LONG / "jobs.json").write_text(json.dumps(jobs, indent=2))


# --------------------------------------------------------------------------------------------------------------------
# collect, compare, write


def _load(case: Path) -> dict | None:
    p = Path(case) / "results.json"
    return json.loads(p.read_text()) if p.exists() else None


def straight_front(level1: dict) -> tuple[np.ndarray, np.ndarray]:
    """Level 1 CFD front (ΔP, R_th), ascending ΔP; the straight-fin reference between the 3D straight points."""
    pts = [(d["cfd"]["dp_pa"], d["cfd"]["r_th_k_per_w"]) for d in level1["designs"]]
    idx = pareto_indices([p[0] for p in pts], [p[1] for p in pts])
    f = sorted(pts[i] for i in idx)
    return np.array([p[0] for p in f]), np.array([p[1] for p in f])


def interp_loglog(x, xs, ys) -> float | None:
    lx, ly = np.log(xs), np.log(ys)
    if not (lx[0] <= np.log(x) <= lx[-1]):
        return None
    return float(np.exp(np.interp(np.log(x), lx, ly)))


def _foam_list(text: str, start: int) -> tuple[list[str], int]:
    """Entries of an ascii OpenFOAM list 'N ( ... )' beginning at or after `start`."""
    import re
    m = re.compile(r"(\d+)\s*\(").search(text, start)
    n, i = int(m.group(1)), m.end()
    depth, j, items, cur = 1, i, [], i
    while True:
        c = text[j]
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
            if depth == 0:
                break
        j += 1
    body = text[i:j]
    if "(" in body:
        items = re.findall(r"\(([^()]*)\)", body)
    else:
        items = body.split()
    assert len(items) == n, (len(items), n)
    return items, j


def heat_split(case: Path) -> dict:
    """Heat into the coolant through the floor (base top between the pins) and through the pin walls, from the
    wallHeatFlux field on fluid_to_solid (written by `chtMultiRegionSimpleFoam -postProcess`) and the face areas."""
    import re
    case = Path(case)
    s = json.loads((case / "setup.json").read_text())
    mesh = case / "constant/fluid/polyMesh"
    pts_t = (mesh / "points").read_text()
    pts = np.array([[float(v) for v in p.split()] for p in _foam_list(pts_t, pts_t.index("// *"))[0]])
    faces_t = (mesh / "faces").read_text()
    bnd = (mesh / "boundary").read_text()
    m = re.search(r"fluid_to_solid\s*\{[^}]*?nFaces\s+(\d+);\s*startFace\s+(\d+);", bnd, re.S)
    n, start = int(m.group(1)), int(m.group(2))
    faces = _foam_list(faces_t, faces_t.index("// *"))[0][start:start + n]
    times = sorted((p for p in case.iterdir() if p.name.replace(".", "").isdigit() and float(p.name) > 0), key=lambda p: float(p.name))
    wf = (times[-1] / "fluid/wallHeatFlux").read_text()
    k = wf.index("fluid_to_solid")
    q = np.array([float(v) for v in _foam_list(wf, wf.index("value", k))[0]])
    area, nz = np.zeros(n), np.zeros(n)
    for i, f in enumerate(faces):
        p = pts[[int(v) for v in f.split()]]
        a = 0.5 * np.cross(p[1:-1] - p[0], p[2:] - p[0]).sum(0)
        area[i] = np.linalg.norm(a)
        nz[i] = abs(a[2]) / area[i]
    floor = nz > 0.99
    q_floor, q_pin = float((q * area)[floor].sum()), float((q * area)[~floor].sum())
    tot = q_floor + q_pin
    return {"q_floor_w": q_floor, "q_pins_w": q_pin, "floor_fraction": q_floor / tot, "floor_area_m2": float(area[floor].sum()),
            "pin_wall_area_m2": float(area[~floor].sum()), "q_cell_w": s["case"]["q_cell_w"]}


def run_heat_split(case: Path) -> dict:
    case = Path(case).resolve()
    (case / "Allrun.heat").write_text("set -e\ncd /case\nchtMultiRegionSimpleFoam -postProcess -func wallHeatFlux -region fluid "
                                      "-latestTime > log.heatSplit 2>&1\n")
    cmd = lp.docker_cmd(case)
    cmd[-1] = "bash /case/Allrun.heat"
    subprocess.run(cmd, check=True)
    out = heat_split(case)
    (case / "heat_split.json").write_text(json.dumps(out, indent=2))
    return out


def collect() -> dict:
    pv = json.loads(PLANVIEW.read_text())
    level1, req, fluid = setup()
    st3d = json.loads(STRAIGHT3D.read_text())
    fx, fy = straight_front(level1)
    pv_st = [r for r in pv["designs"] if r["topology"] == "straight" and r["status"] == "converged"]
    from level2.planview_sweep import aspect_exclusion
    pv_st = [r for r in pv_st if not aspect_exclusion(r)]
    sidx = pareto_indices([r["result"]["dp_pa"] for r in pv_st], [r["result"]["r_th_k_per_w"] for r in pv_st])
    pvf = sorted((pv_st[i]["result"]["dp_pa"], pv_st[i]["result"]["r_th_k_per_w"]) for i in sidx)
    pvx, pvy = np.array([p[0] for p in pvf]), np.array([p[1] for p in pvf])
    out = {"image": IMAGE, "solver": "chtMultiRegionSimpleFoam (OpenFOAM v2406)", "planview": str(PLANVIEW),
           "level1": str(LEVEL1), "gate_limits": GATES, "factors": {"coarse": COARSE, "production": PRODUCTION, "fine": FINE},
           "n_iter": {str(k): v for k, v in N_ITER.items()}, "mesh_sizes": {"plan": lp.PLAN, "z": lp.Z_CELLS, "grading": lp.GRADING},
           "pick": pv["summary"]["recommended"]["design_id"]}
    # P0
    pick = next(r for r in pv["designs"] if r["design_id"] == out["pick"])
    p0 = {f: _load(case_dir(out["pick"], "p0", f)) for f in P0_FACTORS}
    if all(p0.values()):
        prod, fine = p0[PRODUCTION], p0[FINE]
        ref = pick["result"]
        dv = lambda a, b: a / b - 1
        g = {"dpdx vs plan view": dv(prod["dpdx_pa_per_m"], ref["dpdx_pa_per_m"]),
             "h vs plan view": dv(prod["h_w_m2k"], ref["h_w_m2k"]),
             "dpdx production vs fine": dv(prod["dpdx_pa_per_m"], fine["dpdx_pa_per_m"]),
             "h production vs fine": dv(prod["h_w_m2k"], fine["h_w_m2k"])}
        lim = {k: GATES["p0"] if "plan view" in k else GATES["p0_mesh"] for k in g}
        out["p0"] = {"planview": {k: ref[k] for k in ("dpdx_pa_per_m", "h_w_m2k", "nu_dh", "reynolds_b", "cells")},
                     "runs": {str(f): r for f, r in p0.items()},
                     "gates": {k: {"value": v, "limit": lim[k], "pass": abs(v) < lim[k]} for k, v in g.items()},
                     "run_gates": {str(f): run_gates(r) for f, r in p0.items()}}
        out["p0"]["pass"] = all(v["pass"] for v in out["p0"]["gates"].values())
    # parallel check
    ser, par = _load(case_dir(out["pick"], "conj", COARSE)), _load(case_dir(out["pick"], "conj", COARSE, 4))
    if ser and par:
        d = {k: par[k] / ser[k] - 1 for k in ("r_th_k_per_w", "dp_pa")}
        out["parallel_check"] = {"serial": ser, "parallel_np4": par, "difference": d,
                                 "pass": all(abs(v) < GATES["parallel"] for v in d.values())}
    # straight references in 3D
    st = {}
    rec = st3d["recommended"]
    st[rec] = {"source": str(STRAIGHT3D), "level3": st3d["designs"][rec]["levels"]["1"],
               "level1_cfd": st3d["designs"][rec]["level1_cfd"], "role": "Level 1 pick"}
    pvs = straight_planview_pick(pv)
    sid = STRAIGHT_PICK.design_id
    lv = {k: json.loads((STRAIGHT_RUNS / sid / f"level{k}" / "results.json").read_text())
          for k in (1, 2) if (STRAIGHT_RUNS / sid / f"level{k}" / "results.json").exists()}
    if 1 in lv:
        e = {"source": str(STRAIGHT_RUNS / sid), "level3": lv[1], "role": "plan-view straight pick",
             "planview": {k: pvs["result"][k] for k in ("r_th_k_per_w", "dp_pa")},
             "ratio_3d_to_planview": {"r_th": lv[1]["r_th_k_per_w"] / pvs["result"]["r_th_k_per_w"],
                                      "dp": lv[1]["dp_pa"] / pvs["result"]["dp_pa"]},
             "level1_front_at_same_dp": interp_loglog(lv[1]["dp_pa"], fx, fy)}
        if 2 in lv:
            e["fine"] = lv[2]
            e["mesh_change"] = {k: lv[1][k] / lv[2][k] - 1 for k in ("r_th_k_per_w", "dp_pa")}
        if e["level1_front_at_same_dp"]:
            e["vs_level1_front_pct"] = 100 * (lv[1]["r_th_k_per_w"] / e["level1_front_at_same_dp"] - 1)
        st[sid] = e
    out["straight"] = st
    # 3D straight-fin reference: Pareto front of the Level 1 CFD designs and the 3D straight points (Level 1 CFD and the
    # 3D model agree within 0.3% on straight fins, stage one)
    pts = [(d["cfd"]["dp_pa"], d["cfd"]["r_th_k_per_w"]) for d in level1["designs"]]
    pts += [(e["level3"]["dp_pa"], e["level3"]["r_th_k_per_w"]) for e in st.values()]
    ridx = pareto_indices([p[0] for p in pts], [p[1] for p in pts])
    ref = sorted(pts[i] for i in ridx)
    sx, sy = np.array([p[0] for p in ref]), np.array([p[1] for p in ref])
    out["straight_front_level1"] = {"dp_pa": fx.tolist(), "r_th_k_per_w": fy.tolist()}
    out["straight_front_3d_reference"] = {"dp_pa": sx.tolist(), "r_th_k_per_w": sy.tolist()}
    out["straight_front_planview"] = {"dp_pa": pvx.tolist(), "r_th_k_per_w": pvy.tolist()}
    # pins, production
    pins = []
    for r in pin_designs(pv):
        res = _load(case_dir(r["design_id"], "conj", PRODUCTION, PROD_PROCS))
        e = {"design_id": r["design_id"], "params_mm": r["params_mm"], "is_pick": r["design_id"] == out["pick"],
             "planview": {k: r["result"][k] for k in ("r_th_k_per_w", "dp_pa", "dpdx_pa_per_m", "h_w_m2k",
                                                        "reynolds_pin_max", "fin_efficiency", "a_floor_m2", "a_fin_m2",
                                                        "r_base", "r_conv", "r_caloric")}}
        pvr = e["planview"]
        e["planview"]["floor_fraction"] = pvr["a_floor_m2"] / (pvr["a_floor_m2"] + pvr["fin_efficiency"] * pvr["a_fin_m2"])
        if res:
            e["level3"] = res
            e["gates"] = run_gates(res)
            e["ratio_3d_to_planview"] = {"r_th": res["r_th_k_per_w"] / r["result"]["r_th_k_per_w"],
                                          "dp": res["dp_pa"] / r["result"]["dp_pa"]}
            s3 = interp_loglog(res["dp_pa"], sx, sy)
            e["straight_3d_front_at_same_dp"] = s3
            e["vs_straight_3d_pct"] = None if s3 is None else 100 * (res["r_th_k_per_w"] / s3 - 1)
            s1 = interp_loglog(res["dp_pa"], fx, fy)
            e["vs_level1_front_pct"] = None if s1 is None else 100 * (res["r_th_k_per_w"] / s1 - 1)
            if (Path(case_dir(r["design_id"], "conj", PRODUCTION, PROD_PROCS)) / "heat_split.json").exists():
                e["heat_split"] = json.loads((case_dir(r["design_id"], "conj", PRODUCTION, PROD_PROCS) / "heat_split.json").read_text())
        stp = interp_loglog(r["result"]["dp_pa"], pvx, pvy)
        e["planview_vs_straight_pct"] = None if stp is None else 100 * (r["result"]["r_th_k_per_w"] / stp - 1)
        pins.append(e)
    out["pins"] = pins
    # equal-ΔP table: pin 3D front (log-log through the three designs) against the 3D straight reference
    done = [(e["level3"]["dp_pa"], e["level3"]["r_th_k_per_w"]) for e in pins if "level3" in e]
    done = sorted(done[i] for i in pareto_indices([d[0] for d in done], [d[1] for d in done])) if done else []
    out["pins_3d_front_dp_pa"] = [d[0] for d in done]
    table = []
    if len(done) >= 2:
        px, py = np.array([d[0] for d in done]), np.array([d[1] for d in done])
        for lvl in np.linspace(px[0], px[-1], 3):
            rp, rs = interp_loglog(lvl, px, py), interp_loglog(lvl, sx, sy)
            table.append({"dp_pa": float(lvl), "pins_3d": rp, "straight_3d_reference": rs,
                          "difference_pct": 100 * (rp / rs - 1) if rp and rs else None})
    out["equal_dp_3d"] = table
    out["long_jobs"] = collect_long()
    return out


def collect_long() -> dict:
    from .status import job_state
    jf = LONG / "jobs.json"
    if not jf.exists():
        return {}
    jobs = json.loads(jf.read_text())
    out = {}
    for name, spec in jobs.items():
        st = job_state(name, spec)
        out[name] = {"what": spec["what"], "state": st["state"], "latest": st.get("latest")}
        case = Path(spec["case"])
        if st["state"] == "finished":
            try:
                if name == "mesh_gate":
                    r = analyse(case)
                    out[name]["result"] = r
                    out[name]["gates"] = run_gates(r)
                else:
                    out[name]["result"] = analyse_transient(case)
            except Exception as exc:  # report, do not hide
                out[name]["error"] = repr(exc)
    return out


def analyse_transient(case: Path) -> dict:
    """Time series of R_th, ΔP and wake u_y after the perturbation; means over the last half of the run."""
    case = Path(case)
    s = json.loads((case / "setup.json").read_text())
    cell = _cell(s)
    text = (case / "log.solver").read_text()
    log = parse_solver_log(text)
    times = [float(t) for t in __import__("re").findall(r"^Time = (\S+)", text, __import__("re").M)]
    base = parse_dat(next((case / "postProcessing/solid/baseTt").glob("*/surfaceFieldValue.dat")))
    n = min(len(base), len(log["bulk_t"]), len(log["dpdx_pa_per_m"]))
    t = np.array([b[0] for b in base[:n]])
    r_th = np.array([lp.r_th_pins(cell, base[i][1], log["bulk_t"][i]) for i in range(n)])
    dp = np.array([cell.dp_core_pa(g) for g in log["dpdx_pa_per_m"][:n]])
    probe = parse_dat(next((case / "postProcessing/wake/fluid").glob("*/U")))
    tp = np.array([p[0] for p in probe])
    uy = np.array([[p[1 + 3 * k + 1] for k in range(3)] for p in probe])
    half = t >= t[-1] / 2
    hp = tp >= tp[-1] / 2
    return {"end_time_s": float(t[-1]), "n_steps": int(n), "r_th_mean_last_half": float(r_th[half].mean()),
            "dp_mean_last_half": float(dp[half].mean()), "r_th_range_last_half": float(np.ptp(r_th[half])),
            "dp_range_last_half": float(np.ptp(dp[half])), "uy_wake_max_abs_last_half": float(np.abs(uy[hp]).max()),
            "uy_wake_max_abs_first_step": float(np.abs(uy[0]).max()), "ubar_m_s": s["case"]["ubar_m_s"],
            "series": {"t": t[:: max(1, n // 400)].tolist(), "r_th": r_th[:: max(1, n // 400)].tolist(),
                       "dp": dp[:: max(1, n // 400)].tolist()}}


# --------------------------------------------------------------------------------------------------------------------
# figures and summary


def plot_front(data: dict, out: Path) -> Path:
    from tools.report import draw_requirements
    level1 = json.loads(LEVEL1.read_text())
    budget = level1["dp_budget_pa"] / 1e3
    thr, obj = level1["r_th_threshold_k_per_w"] * 1e3, level1["r_th_objective_k_per_w"] * 1e3
    fig, ax = ps.figure(5.6)
    x_lo, x_hi, y_lo, y_hi = 0.9, 6.0, 11.0, 23.0
    draw_requirements(ax, x_lo, y_lo, budget, thr, obj, y_hi, labels_right=True)
    f1, f2 = data["straight_front_3d_reference"], data["straight_front_planview"]
    ax.plot(np.array(f1["dp_pa"]) / 1e3, np.array(f1["r_th_k_per_w"]) * 1e3, "-", color=ps.BLUE, lw=1.6, zorder=2,
            label="Straight fins, 3D reference front")
    ax.plot(np.array(f2["dp_pa"]) / 1e3, np.array(f2["r_th_k_per_w"]) * 1e3, "--", color=ps.BLUE, lw=1.0, alpha=0.6, zorder=2,
            label="Straight fins, plan-view front")
    for sid, e in data["straight"].items():
        x, y = e["level3"]["dp_pa"] / 1e3, e["level3"]["r_th_k_per_w"] * 1e3
        ax.plot([x], [y], "s", color=ps.BLUE, ms=7, mec="white", mew=0.6, zorder=4)
        if "planview" in e:
            ax.plot([e["planview"]["dp_pa"] / 1e3], [e["planview"]["r_th_k_per_w"] * 1e3], "s", mfc="none", mec=ps.BLUE,
                    ms=7, mew=1.3, zorder=4)
            ax.annotate("", xy=(x, y), xytext=(e["planview"]["dp_pa"] / 1e3, e["planview"]["r_th_k_per_w"] * 1e3),
                        arrowprops=dict(arrowstyle="->", color=ps.BLUE, lw=0.8, alpha=0.7))
    pins = data["pins"]
    for e in pins:
        pvx, pvy = e["planview"]["dp_pa"] / 1e3, e["planview"]["r_th_k_per_w"] * 1e3
        ax.plot([pvx], [pvy], "o", mfc="none", mec=ps.ACCENT, ms=7, mew=1.3, zorder=4)
        if "level3" in e:
            x, y = e["level3"]["dp_pa"] / 1e3, e["level3"]["r_th_k_per_w"] * 1e3
            ax.plot([x], [y], "o", color=ps.ACCENT, ms=7, mec="white", mew=0.6, zorder=5)
            ax.annotate("", xy=(x, y), xytext=(pvx, pvy), arrowprops=dict(arrowstyle="->", color=ps.ACCENT, lw=0.8, alpha=0.7))
            if e["is_pick"]:
                ax.scatter([x], [y], s=230, facecolors="none", edgecolors=ps.INK, linewidths=1.6, zorder=6)
    done = [(e["level3"]["dp_pa"] / 1e3, e["level3"]["r_th_k_per_w"] * 1e3) for e in pins if "level3" in e]
    if done:
        done = sorted(done[i] for i in pareto_indices([d[0] for d in done], [d[1] for d in done]))
        ax.plot([d[0] for d in done], [d[1] for d in done], "-", color=ps.ACCENT, lw=1.4, zorder=3,
                label="Circular pins, 3D front")
    pv = sorted((e["planview"]["dp_pa"] / 1e3, e["planview"]["r_th_k_per_w"] * 1e3) for e in pins)
    ax.plot([d[0] for d in pv], [d[1] for d in pv], "--", color=ps.ACCENT, lw=1.0, alpha=0.6, zorder=3,
            label="Circular pins, plan view")
    ax.plot([], [], "o", color=ps.INK, ms=6, label="3D (solid markers)")
    ax.plot([], [], "o", mfc="none", mec=ps.INK, ms=6, label="Plan view (hollow)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ps.plain_log_ticks(ax, x_lo, x_hi, "x", ticks=(1, 1.5, 2, 3, 5))
    ps.plain_log_ticks(ax, y_lo, y_hi, "y", ticks=(11, 12, 13, 14, 15, 16, 18, 20))
    ax.set_xlabel("Pressure drop ΔP (kPa)")
    ax.set_ylabel("Thermal resistance R_th (mK/W)")
    ax.legend(loc="lower left", fontsize=ps.SMALL, ncol=2)
    ax.set_title("Circular pins against straight fins, 3D and plan view")
    ps.save(fig, out)
    return out


def _read_raw(path: Path) -> np.ndarray:
    rows = [ln.split() for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    return np.array(rows, dtype=float)


def plot_fields(case: Path, out: Path) -> Path | None:
    """Speed on the mid-height plane (tiled over two periods, mirrored to a full pitch) and temperature on the y = 0 plane
    through pin A and the base, from the sampled surfaces of the production run."""
    from matplotlib.tri import Triangulation
    case = Path(case)
    s = json.loads((case / "setup.json").read_text())
    p = s["params_mm"]
    pp = case / "postProcessing"
    try:
        u_file = sorted(pp.glob("planes/*/U_midZ.raw"), key=lambda f: float(f.parent.name))[-1]
        ts_file = sorted(pp.glob("solidPlanes/*/T_lo.raw"), key=lambda f: float(f.parent.name))[-1]
    except IndexError:
        return None
    u = _read_raw(u_file)
    ts = _read_raw(ts_file)
    x, y = u[:, 0] * 1e3, u[:, 1] * 1e3
    speed = np.linalg.norm(u[:, 3:6], axis=1)
    lat, per = p["s_t"] / 2, 2 * p["s_l"]
    X, Y, S = [], [], []
    for k in range(2):
        for sign in (1, -1):
            X.append(x + k * per)
            Y.append(sign * y)
            S.append(speed)
    X, Y, S = map(np.concatenate, (X, Y, S))
    r = p["d"] / 2
    pins = [(cx + k * per, sy * cy) for k in range(-1, 3) for cx, cy in ((0, 0), (p["s_l"], lat)) for sy in (1, -1)]
    tri = Triangulation(X, Y)
    cx = X[tri.triangles].mean(1)
    cy = Y[tri.triangles].mean(1)
    mask = np.zeros(len(cx), bool)
    for px, py in pins:
        mask |= np.hypot(cx - px, cy - py) < r * 1.02
    edges = np.max(np.hypot(np.diff(X[tri.triangles][:, [0, 1, 2, 0]], axis=1), np.diff(Y[tri.triangles][:, [0, 1, 2, 0]], axis=1)), 1)
    mask |= edges > 4 * np.median(edges)
    tri.set_mask(mask)
    fig, axes = ps.figure(5.4, 2, 1, gridspec_kw={"height_ratios": [1.0, 1.0]})
    cs = axes[0].tricontourf(tri, S, levels=ps.field_levels(0, float(S.max())), cmap="viridis")
    ps.colorbar(fig, cs, axes[0], "Speed (m/s)", 0, float(S.max()), nbins=4, min_ticks=3, shrink=0.95)
    from matplotlib.patches import Circle
    for px, py in pins:
        axes[0].add_patch(Circle((px, py), r, color=ps.COPPER, zorder=3))
    axes[0].set_xlim(-p["s_l"] / 2, -p["s_l"] / 2 + 2 * per)
    axes[0].set_ylim(-lat, lat)
    axes[0].set_aspect("equal")
    axes[0].set_xlabel("x, flow direction (mm)")
    axes[0].set_ylabel("y (mm)")
    axes[0].set_title("Speed at mid-height (plan view of two periods, one lateral pitch)", fontsize=11.5, loc="left")
    ps.no_grid(axes[0])
    # temperature on the y = 0 plane: solid only (pin A and the base), relative to the period-mean bulk
    xs, zs, T = ts[:, 0] * 1e3, ts[:, 2] * 1e3, ts[:, 3]
    log = parse_solver_log((case / "log.solver").read_text())
    tb = log["bulk_t"][-1]
    # drawn with z across the page (equal scales): the chip interface on the left, the cover on the right
    tri2 = Triangulation(zs, xs)
    c2 = tri2.triangles
    zc, xc = zs[c2].mean(1), xs[c2].mean(1)
    in_pin = (zs[c2] > 2.0).any(1) & (np.abs(xs[c2]) > r).any(1)  # bridges between the base top and the pin side
    tri2.set_mask(in_pin | ((zc > 2.0) & (np.abs(xc) > 0.98 * r)))
    dT = T - tb
    cs2 = axes[1].tricontourf(tri2, dT, levels=ps.field_levels(float(dT.min()), float(dT.max())), cmap="inferno")
    ps.colorbar(fig, cs2, axes[1], "T − bulk T (K)", float(dT.min()), float(dT.max()), nbins=4, min_ticks=3, shrink=0.95)
    axes[1].set_xlim(0, 2.0 + p["h_fin"])
    axes[1].set_ylim(xs.min(), xs.max())
    axes[1].set_aspect("equal")
    axes[1].set_xlabel("z, up from the chip interface (mm); base 0–2 mm, pin 2–6 mm")
    axes[1].set_ylabel("x (mm)")
    axes[1].set_title("Temperature in the base and pin A on the plane y = 0", fontsize=11.5, loc="left")
    axes[1].set_facecolor("#eef2f7")
    ps.no_grid(axes[1])
    fig.suptitle("Circular pins, " + f"{p['d']:g} mm diameter, {p['s_t']:g} × {p['s_l']:g} mm pitch, 3D", fontsize=12.5)
    ps.save(fig, out)
    return out


def _pct(v):
    return "—" if v is None else f"{v:+.1f}%"


def write_summary(data: dict, out: Path, fields_png: bool) -> Path:
    a = ["# Level 3 stage two: circular pins against straight fins in 3D\n"]
    a.append(f"OpenFOAM v2406 `chtMultiRegionSimpleFoam` in Docker (`{data['image']}`), steady laminar, constant properties "
             "(PG25 at 40 °C, copper k = 390 W/m·K). Staggered circular pins as in the plan view (`level2/planforms.py`): one "
             "stagger period streamwise (periodic, mean velocity held), half a transverse pitch with symmetry planes, full pin "
             "height with a no-slip floor and an adiabatic no-slip cover on the pin tips, the 2 mm copper base under 1 kW "
             "uniform over the 60 × 60 mm footprint, conjugate conduction in pins and base. Flow 1.5 LPM over the 60 mm width "
             "(superficial velocity Q / (W H), as the plan view). R_th = (base-underside area average − period-mean bulk) / Q "
             "plus half the coolant rise, ΔP = dp/dx × 60 mm, both as stage one. Mesh: gmsh quadrilaterals in plan extruded "
             "to hexahedra, two conformal regions (`level3/pins.py`).\n")
    # P0
    if "p0" in data:
        p0 = data["p0"]
        a.append("## P0: plan-view cross-check (fluid only, slip floor and cover)\n")
        a.append(f"The pick `{data['pick']}` as a {lp.P0_SLAB_MM} mm slab one cell thick, symmetry (slip) floor and cover, no "
                 "base, uniform flux on the pin walls: the problem the plan-view solver solves. h = wall flux / (period-mean "
                 "wall T − period-mean velocity-weighted bulk T), the plan-view definition.\n")
        a.append("| Mesh | Cells | dp/dx (Pa/m) | h (W/m²K) | Final residual | Energy balance |\n|---|---|---|---|---|---|")
        for f, r in p0["runs"].items():
            tag = " (production)" if int(f) == PRODUCTION else " (2× finer in plan)"
            a.append(f"| factor {f}{tag} | {r['n_cells']:,} | {r['dpdx_pa_per_m']:.0f} | {r['h_w_m2k']:.0f} | "
                     f"{r['max_final_residual_used']:.1e} | {100 * r['energy_balance_error']:+.1e}% |")
        pvr = p0["planview"]
        a.append(f"| plan view (Taylor-Hood, {pvr['cells']:g} cells per gap) | — | {pvr['dpdx_pa_per_m']:.0f} | {pvr['h_w_m2k']:.0f} | — | — |\n")
        a.append("| Gate | Value | Limit | Result |\n|---|---|---|---|")
        a += [f"| {k} | {100 * v['value']:+.2f}% | ±{100 * v['limit']:g}% | {'pass' if v['pass'] else 'FAIL'} |"
              for k, v in p0["gates"].items()]
        a.append("")
    if "parallel_check" in data:
        pc = data["parallel_check"]
        a.append("## Serial against parallel (coarse conjugate pick)\n")
        a.append(f"Factor {COARSE} mesh ({pc['serial']['n_cells']:,} cells), {pc['serial']['iterations']} iterations: serial "
                 f"R_th {pc['serial']['r_th_k_per_w'] * 1e3:.4f} mK/W, ΔP {pc['serial']['dp_pa']:.2f} Pa; four processes "
                 f"(scotch) {pc['parallel_np4']['r_th_k_per_w'] * 1e3:.4f} mK/W, {pc['parallel_np4']['dp_pa']:.2f} Pa "
                 f"(difference {100 * pc['difference']['r_th_k_per_w']:+.1e}% and {100 * pc['difference']['dp_pa']:+.1e}%, "
                 f"limit {100 * GATES['parallel']:g}%: {'pass' if pc['pass'] else 'FAIL'}). Run times "
                 f"{pc['serial']['seconds']:.0f} s and {pc['parallel_np4']['seconds']:.0f} s.")
        pick3 = next((e.get("level3") for e in data["pins"] if e["is_pick"]), None)
        if pick3:
            a.append(f"The same coarse mesh (half the production resolution in every direction) differs from the production run "
                     f"by {100 * (pc['serial']['r_th_k_per_w'] / pick3['r_th_k_per_w'] - 1):+.2f}% on R_th and "
                     f"{100 * (pc['serial']['dp_pa'] / pick3['dp_pa'] - 1):+.2f}% on ΔP; the gate is production against 2× finer "
                     "(detached job below).")
        a.append("")
    # pins table
    a.append("## Pins in 3D against the plan view\n")
    a.append("| Design | Re_d,max | Plan view R_th (mK/W) | 3D R_th (mK/W) | 3D / plan | Plan view ΔP (kPa) | 3D ΔP (kPa) | 3D / plan | "
             "Cells | Residual | Flat (last 10%) | Energy balance | Run time (min) | Peak memory (GB) |")
    a.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for e in data["pins"]:
        pvr = e["planview"]
        name = f"`{e['design_id']}`" + (" (pick)" if e["is_pick"] else "")
        if "level3" not in e:
            a.append(f"| {name} | {pvr['reynolds_pin_max']:.1f} | {pvr['r_th_k_per_w'] * 1e3:.2f} | not run | | "
                     f"{pvr['dp_pa'] / 1e3:.3f} | | | | | | | | |")
            continue
        r, g = e["level3"], e["gates"]
        mem = f"{r['peak_memory_mb'] / 1024:.1f}" if r.get("peak_memory_mb") else "—"
        a.append(f"| {name} | {r['re_d_max']:.1f} | {pvr['r_th_k_per_w'] * 1e3:.2f} | {r['r_th_k_per_w'] * 1e3:.2f} | "
                 f"{e['ratio_3d_to_planview']['r_th']:.3f} | {pvr['dp_pa'] / 1e3:.3f} | {r['dp_pa'] / 1e3:.3f} | "
                 f"{e['ratio_3d_to_planview']['dp']:.3f} | {r['n_cells']:,} | {g['residual']['value']:.1e} "
                 f"({'pass' if g['residual']['pass'] else 'FAIL'}) | {100 * g['flat']['value']:.1e}% "
                 f"({'pass' if g['flat']['pass'] else 'FAIL'}) | {100 * r['energy_balance_error']:+.1e}% "
                 f"({'pass' if g['energy']['pass'] else 'FAIL'}) | {r['seconds'] / 60:.0f} | {mem} |")
    a.append("")
    for e in data["pins"]:
        if "heat_split" in e:
            hs = e["heat_split"]
            a.append(f"Heat path, `{e['design_id']}`: {100 * hs['floor_fraction']:.1f}% of the heat enters the coolant through "
                     f"the floor between the pins in 3D, against {100 * e['planview']['floor_fraction']:.1f}% in the plan-view "
                     "network, which gives the floor the pins' h (ASSUMPTIONS V4). From the wall heat flux on the interface "
                     "(`chtMultiRegionSimpleFoam -postProcess -func wallHeatFlux`), split by face orientation.\n")
    a.append("Gates per run: largest final initial residual (all velocity components, pressure and both energy equations, normalised) "
             f"< {GATES['residual']:g}; R_th and ΔP flat to {100 * GATES['flat']:g}% over the last 10% of iterations; heat "
             f"crossing the fluid-solid interface within {100 * GATES['energy']:g}% of the applied heat. For straight fins the "
             "plan view reads R_th 5.6% and ΔP 7.1% below Level 1 (no floor or cover); the 3D / plan ratios above are the "
             "same bias for pins.\n")
    # straight
    a.append("## Straight fins in 3D (stage-one model)\n")
    a.append("| Design | Role | 3D R_th (mK/W) | 3D ΔP (kPa) | Reference | 3D / reference R_th | 3D / reference ΔP |\n|---|---|---|---|---|---|---|")
    for sid, e in data["straight"].items():
        r = e["level3"]
        if "level1_cfd" in e:
            ref, rn = e["level1_cfd"], "Level 1 CFD"
        else:
            ref, rn = e["planview"], "plan view"
        a.append(f"| `{sid}` | {e['role']} | {r['r_th_k_per_w'] * 1e3:.3f} | {r['dp_pa'] / 1e3:.3f} | {rn} "
                 f"{ref['r_th_k_per_w'] * 1e3:.2f} mK/W, {ref['dp_pa'] / 1e3:.3f} kPa | {r['r_th_k_per_w'] / ref['r_th_k_per_w']:.3f} | "
                 f"{r['dp_pa'] / ref['dp_pa']:.3f} |")
    a.append("")
    sp = data["straight"].get(STRAIGHT_PICK.design_id)
    if sp:
        extra = ""
        if "mesh_change" in sp:
            extra = (f" Its 2× refinement changes R_th by {100 * sp['mesh_change']['r_th_k_per_w']:+.2f}% and ΔP by "
                     f"{100 * sp['mesh_change']['dp_pa']:+.2f}%.")
        if sp.get("vs_level1_front_pct") is not None:
            extra += (f" At its own ΔP the Level 1 CFD front reads {sp['level1_front_at_same_dp'] * 1e3:.3f} mK/W, so this 3D "
                      f"point sits {sp['vs_level1_front_pct']:+.2f}% from that front: 0.40 / 0.40 / 4.0 mm is not on the Level 1 "
                      "grid, and between its grid points (0.42 / 0.39 / 3.75 and 0.39 / 0.39 / 3.75 mm) the Level 1 front is a "
                      "chord that a better straight design lies below.")
        a.append("The 3D straight-fin reference used below is the Pareto front of the committed Level 1 CFD designs "
                 "(`examples/output/gpu_coldplate/results.json`) together with the two 3D straight points, interpolated in "
                 "log-log as the plan-view comparison does; between the 3D points it is the Level 1 CFD front, which stage one "
                 "showed the 3D model reproduces within 0.3%." + extra + "\n")
    # comparison
    a.append("## Pins against straight fins at equal ΔP\n")
    a.append("| Design | 3D ΔP (kPa) | 3D R_th (mK/W) | Straight, 3D reference at that ΔP (mK/W) | 3D difference | Against the Level 1 front alone | Plan-view difference |\n|---|---|---|---|---|---|---|")
    for e in data["pins"]:
        if "level3" not in e:
            continue
        a.append(f"| `{e['design_id']}` | {e['level3']['dp_pa'] / 1e3:.3f} | {e['level3']['r_th_k_per_w'] * 1e3:.2f} | "
                 f"{(e['straight_3d_front_at_same_dp'] or float('nan')) * 1e3:.2f} | {_pct(e['vs_straight_3d_pct'])} | "
                 f"{_pct(e['vs_level1_front_pct'])} | "
                 f"{_pct(e['planview_vs_straight_pct'])} |")
    if data["equal_dp_3d"]:
        a.append("\n| Equal ΔP (kPa) | Pins, 3D (mK/W) | Straight, 3D reference (mK/W) | Difference |\n|---|---|---|---|")
        for t in data["equal_dp_3d"]:
            a.append(f"| {t['dp_pa'] / 1e3:.2f} | {t['pins_3d'] * 1e3:.2f} | {t['straight_3d_reference'] * 1e3:.2f} | "
                     f"{_pct(t['difference_pct'])} |")
        a.append("\nPins between the 3D designs are interpolated in log-log along their 3D Pareto front (designs it "
                 "dominates are left out).")
        done = [e for e in data["pins"] if "level3" in e]
        pick = next((e for e in done if e["is_pick"]), None)
        pv_rng = [e["planview_vs_straight_pct"] for e in done]
        d3 = [e["vs_straight_3d_pct"] for e in done]
        rr = [e["ratio_3d_to_planview"]["r_th"] for e in done]
        pr = [e["ratio_3d_to_planview"]["dp"] for e in done]
        dominated = pick and any(e is not pick and e["level3"]["dp_pa"] <= pick["level3"]["dp_pa"]
                                 and e["level3"]["r_th_k_per_w"] <= pick["level3"]["r_th_k_per_w"] for e in done)
        a.append(f"\n**Answer.** The plan view puts these three pin designs {min(pv_rng):+.1f}% to {max(pv_rng):+.1f}% against "
                 f"straight fins at equal ΔP. In 3D they are {min(d3):+.1f}% to {max(d3):+.1f}%: the pin lead does not survive "
                 "as a lead. In 3D the pins' R_th rises by "
                 f"{100 * (min(rr) - 1):.0f}–{100 * (max(rr) - 1):.0f}% over the plan view (ΔP by {100 * (min(pr) - 1):.0f}–"
                 f"{100 * (max(pr) - 1):.0f}%), against about 6% (R_th) and 6–7% (ΔP) for straight fins. Most of the extra R_th is "
                 "at the floor: the plan-view network gives the floor between the pins the pins' h and so sends "
                 f"{100 * min(e['planview']['floor_fraction'] for e in done):.0f}–"
                 f"{100 * max(e['planview']['floor_fraction'] for e in done):.0f}% of the heat through it, while in 3D the "
                 f"floor carries {100 * min(e.get('heat_split', {}).get('floor_fraction', 1) for e in done):.1f}–"
                 f"{100 * max(e.get('heat_split', {}).get('floor_fraction', 0) for e in done):.1f}%; splitting the convective conductance by heat "
                 "share, about two thirds of what is lost is at the floor and one third at the pins, whose roots, where "
                 "most of their heat leaves, sit in the slow flow near the floor. "
                 + ("The cross-topology pick is dominated in 3D by its lower-ΔP neighbour. " if dominated else "")
                 + "Near 1.5 kPa, where the straight reference is anchored by two 3D points, pins are worse by about "
                 f"{min(d for d, e in zip(d3, done) if e['level3']['dp_pa'] < 1.7e3):.0f}% or more; near 2 kPa the best pin "
                 f"design is {d3[-1]:+.1f}% against a straight reference that is, there, a chord between two Level 1 grid "
                 "designs (the plan-view straight pick lies 3.3% below such a chord), so that difference is not resolved: "
                 "a tie within the uncertainty of the straight reference.")
    a.append("")
    a.append("![R_th against ΔP](pins_vs_straight.png)\n")
    a.append("R_th against channel-core ΔP at 1.5 LPM PG25 at 40 °C, 1 kW over 60 × 60 mm, copper, 2 mm base. Solid markers: "
             "3D conjugate (OpenFOAM); hollow markers: plan view (2D, tall-fin limit, no floor or cover); arrows go from each "
             "design's plan-view value to its 3D value. Orange circles: circular pins (ring: the cross-topology pick); blue "
             "squares: straight fins (Level 1 pick, from stage one, and the plan-view straight pick). Blue solid line: the "
             "3D straight-fin reference (Level 1 CFD designs and the two 3D straight points, Pareto front); dashed lines: "
             "plan-view fronts. Vertical line: 5 kPa budget; horizontal lines: 20 mK/W threshold and 15 mK/W objective.\n")
    if fields_png:
        a.append("![Field of the pick](fields_pick.png)\n")
        a.append("The pick in 3D. Top: speed on the mid-height plane (z = t_b + H / 2), the computed half pitch mirrored about "
                 "its symmetry planes and repeated over two stagger periods; pins in copper colour. Bottom: temperature above "
                 "the period-mean bulk temperature on the symmetry plane through pin A (y = 0), in the base and the pin "
                 "(the fluid on that plane is left blank).\n")
    # long jobs
    a.append("## Long jobs (detached)\n")
    lj = data.get("long_jobs") or {}
    if not lj:
        a.append("Not started.\n")
    for name, j in lj.items():
        line = f"- `{name}`: {j['what']}. State: **{j['state']}**."
        if "result" in j and name == "mesh_gate":
            r = j["result"]
            pick = next(e for e in data["pins"] if e["is_pick"]).get("level3")
            if pick:
                dr, dp = pick["r_th_k_per_w"] / r["r_th_k_per_w"] - 1, pick["dp_pa"] / r["dp_pa"] - 1
                ok = abs(dr) < GATES["mesh"] and abs(dp) < GATES["mesh"]
                line += (f" Fine mesh ({r['n_cells']:,} cells): R_th {r['r_th_k_per_w'] * 1e3:.3f} mK/W, ΔP "
                         f"{r['dp_pa'] / 1e3:.3f} kPa; production against fine {100 * dr:+.2f}% and {100 * dp:+.2f}% "
                         f"(limit ±2%: {'pass' if ok else 'FAIL'}).")
        if "result" in j and name == "transient":
            r = j["result"]
            line += (f" Ran to t = {r['end_time_s']:.3f} s; over the last half, R_th {r['r_th_mean_last_half'] * 1e3:.3f} mK/W "
                     f"(range {r['r_th_range_last_half'] * 1e3:.2e} mK/W), ΔP {r['dp_mean_last_half'] / 1e3:.3f} kPa (range "
                     f"{r['dp_range_last_half']:.2e} Pa); largest |u_y| on the wake centreline {r['uy_wake_max_abs_last_half']:.2e} m/s "
                     f"against {r['uy_wake_max_abs_first_step']:.2e} m/s just after the perturbation.")
        if "error" in j:
            line += f" Post-processing error: {j['error']}."
        if j.get("latest"):
            line += f" Latest: {j['latest']}."
        a.append(line)
    a.append("")
    out.write_text("\n".join(a))
    return out


def summarise() -> dict:
    data = collect()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "results.json").write_text(json.dumps(data, indent=2))
    plot_front(data, OUT / "pins_vs_straight.png")
    pick = next((e for e in data["pins"] if e["is_pick"] and "level3" in e), None)
    fields = False
    if pick:
        case = case_dir(pick["design_id"], "conj", PRODUCTION, pick["level3"]["n_procs"])
        fields = plot_fields(case, OUT / "fields_pick.png") is not None
    write_summary(data, OUT / "summary.md", fields)
    return data


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--p0", action="store_true")
    ap.add_argument("--parallel-check", action="store_true")
    ap.add_argument("--steady", action="store_true")
    ap.add_argument("--launch-long", action="store_true")
    ap.add_argument("--prepare-long", action="store_true", help="write the long-job cases without starting them")
    ap.add_argument("--summary-only", action="store_true")
    ap.add_argument("--finish", nargs="+", type=Path, help="analyse finished case directories (driver was not running)")
    a = ap.parse_args(argv)
    pv = json.loads(PLANVIEW.read_text())
    if a.p0:
        step_p0(pv)
    if a.parallel_check:
        step_parallel_check(pv)
    if a.steady:
        step_steady(pv)
    if a.finish:
        for c in a.finish:
            r = finish(c)
            print(f"{c}: R_th {r['r_th_k_per_w'] * 1e3:.4f} mK/W, ΔP {r['dp_pa']:.2f} Pa, {r['seconds']:.0f} s")
    if a.prepare_long:
        for name, spec in long_jobs(pv).items():
            print(prepare_long(name, spec), flush=True)
        return 0
    if a.launch_long:
        launch_long(pv)
        return 0
    data = summarise()
    if "p0" in data and not data["p0"]["pass"]:
        print("P0 GATE FAILED", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
