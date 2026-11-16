"""Entrance effects in the straight-fin study (issue #12): 1D developing-flow option against a resolved developing solve.

1D: `tools.screen_1d.analyze_1d(..., entrance_effects=True)` on every CFD-verified straight design of the Level 1 study
(read from examples/output/gpu_coldplate/results.json). Values are reported two ways: the 1D screen on its own, and the
Level 1 CFD value plus the 1D increment (CFD-anchored), which is what the table in summary.md uses.

Resolved check (tall-fin limit): the plan-view solver (`level2.planar`, the same Taylor-Hood Navier-Stokes + energy solve as
`level2.planview.run_developing`) on one straight half channel over the full 60 mm plate with a uniform inlet (velocity
and temperature both developing), and again with a parabolic inlet (thermal entrance only). Uniform flux on the fin wall.
The fins are of infinite height: no floor, no cover, so the channel is a parallel-plate channel of gap w heated on both
sides. Its length-mean h = q'' / (mean wall T - mean bulk T over 0..60 mm) and the pressure drop from the inlet plane to
x = 60 mm are compared with the developed-periodic values (`level2.planview.run_design`) and with the 1D correction
evaluated for the same parallel-plate channel (aspect ratio 0.01, as Muzychka & Yovanovich model parallel plates, with the
screen's developed Nu in the plate limit, 8.235 on 2w, and K(inf) for parallel plates). A plenum-to-plenum solve gives the
inlet and outlet losses, which the 1D option reports separately.

    uv run python -m level2.entrance              # everything, about 5 min; writes examples/output/entrance_effects/
    uv run python -m level2.entrance --plot-only  # redraw the figure and summary from results.json
"""
import argparse
import json
import math
import sys
import time
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np

from tools import plotstyle as ps
from tools.fluids import fluid_props
from tools.pareto import recommend_marginal_return_objective
from tools.screen_1d import (K_INF_PLATES, analyze_1d, effective_nu_factor, f_app_re_sqrt_a, f_re_sqrt_a, k_contraction, k_expansion,
                             nu_fd_sqrt_a, nu_h1)
from tools.types import Design, Requirements

from .planar import line_mean, solve_energy, solve_flow, window_means
from .planforms import mesh_planform, straight
from .planview import B_M, network, run_design, superficial_velocity

ROOT = Path(__file__).resolve().parents[1]
LEVEL1 = ROOT / "examples/output/gpu_coldplate/results.json"
OUT = ROOT / "examples/output/entrance_effects"
FLOW_LPM = 1.5
LENGTH_MM = 60.0
WINDOW_MM = 5.0  # reporting windows of the developing solve
CELLS, CELLS_FINE = 8.0, 16.0  # developing solve: cells across the half channel (wall cell w / (2 cells))
CELLS_PERIODIC = 10.0  # developed-periodic reference, as in the plan-view study
LEAD_IN_MM = 1.0  # slip lead-in upstream of the fin leading edge (uniform-inlet case)
PLENUM_IN_MM, PLENUM_OUT_MM = 2.0, 10.0  # plenum case: plenums as tall as the fins, before and after the fin array
PLATE_EPS = 0.01  # Muzychka & Yovanovich model a parallel-plate channel as a rectangular duct of aspect ratio 0.01
GATE_MESH, GATE_ENERGY, GATE_AGREEMENT = 0.02, 0.01, 0.10


def resolved_designs(rows: list[dict], recommended: str) -> list[str]:
    """Designs for the resolved check, taken from the Level 1 results (nothing hard-coded): the widest-channel (lowest ΔP)
    design, the recommended design and its neighbours in ΔP order, and the narrowest-channel (highest ΔP) design."""
    order = sorted(rows, key=lambda r: r["cfd_anchored"]["dp_fd"])
    ids = [r["design_id"] for r in order]
    i = ids.index(recommended)
    picks = [ids[0]] + ids[max(i - 1, 0):i + 2] + [ids[-1]]
    return [d for k, d in enumerate(picks) if d not in picks[:k]]


def requirements(level1: dict) -> Requirements:
    return Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, 40.0, level1["dp_budget_pa"], "PG25",
                        r_th_threshold_k_per_w=level1["r_th_threshold_k_per_w"],
                        r_th_objective_k_per_w=level1["r_th_objective_k_per_w"],
                        marginal_return_cutoff=level1["marginal_return_cutoff"],
                        min_dp_ratio=level1.get("min_dp_ratio", 1.10))


# --------------------------------------------------------------------------------------------------------------------
# 1D on the CFD-verified designs


def one_d_table(level1: dict, req: Requirements, fluid) -> list[dict]:
    rows = []
    for entry in level1["designs"]:
        g, cfd = entry["one_d"]["design"], entry["cfd"]
        r = analyze_1d(Design(g["w_ch_m"], g["t_fin_m"], g["h_fin_m"], g["flow_lpm"]), req, fluid, entrance_effects=True)
        fd_check = abs(r.r_th_fully_developed_k_per_w / entry["one_d"]["r_th_k_per_w"] - 1) + \
            abs(r.dp_fully_developed_pa / entry["one_d"]["dp_pa"] - 1)
        r_cfd, dp_cfd = cfd["r_th_k_per_w"], cfd["dp_pa"]
        d_th = r.r_th_thermal_entrance_k_per_w - r.r_th_fully_developed_k_per_w
        d_comb = r.r_th_k_per_w - r.r_th_fully_developed_k_per_w
        rows.append({
            "design_id": cfd["design_id"], "w_mm": g["w_ch_m"] * 1e3, "t_mm": g["t_fin_m"] * 1e3, "h_mm": g["h_fin_m"] * 1e3,
            "n_channels": r.n_channels, "reynolds": r.reynolds, "x_star": r.x_star, "l_plus_dh": r.l_plus_dh,
            "sigma": r.sigma, "k_contraction": r.k_contraction, "k_expansion": r.k_expansion,
            "nu_factor_thermal": r.nu_factor_thermal, "nu_factor_combined": r.nu_factor_combined,
            "dp_factor": r.dp_factor, "k_incremental": r.k_incremental, "dynamic_pressure_pa": r.dynamic_pressure_pa,
            "fd_matches_level1": fd_check < 1e-12,
            "one_d": {"r_th_fd": r.r_th_fully_developed_k_per_w, "r_th_thermal": r.r_th_thermal_entrance_k_per_w,
                      "r_th_combined": r.r_th_k_per_w, "dp_fd": r.dp_fully_developed_pa, "dp_core_developing": r.dp_pa,
                      "dp_hydrodynamic_entrance": r.dp_hydrodynamic_entrance_pa, "dp_inlet_outlet": r.dp_inlet_outlet_pa},
            "cfd_anchored": {"r_th_fd": r_cfd, "r_th_thermal": r_cfd + d_th, "r_th_combined": r_cfd + d_comb,
                             "dp_fd": dp_cfd, "dp_combined_core": dp_cfd + r.dp_hydrodynamic_entrance_pa,
                             "dp_inlet_outlet_separate": r.dp_inlet_outlet_pa},
        })
    return rows


def r_th_with_h_factor(row: dict, req: Requirements, fluid, factor: float) -> float:
    """1D screen R_th of a design with its developed h multiplied by `factor` (same network as tools.screen_1d)."""
    fd = analyze_1d(Design(row["w_mm"] * 1e-3, row["t_mm"] * 1e-3, row["h_mm"] * 1e-3, FLOW_LPM), req, fluid)
    w, t, h = row["w_mm"] * 1e-3, row["t_mm"] * 1e-3, row["h_mm"] * 1e-3
    htc = fd.nusselt * factor * fluid.k / (2 * w * h / (w + h))
    m = math.sqrt(2 * htc / (req.solid_k * t))
    eta = math.tanh(m * h) / (m * h)
    return fd.r_base + fd.r_caloric + 1 / (fd.n_channels * req.flow_length_m * htc * (w + 2 * eta * h))


def selection(rows: list[dict], req: Requirements, r_key: str, dp_key: str) -> dict:
    """The study's selection rule (SPEC Amendment 4) on the CFD-anchored values."""
    cands = [(r["cfd_anchored"][dp_key], r["cfd_anchored"][r_key]) for r in rows]
    rec = recommend_marginal_return_objective(cands, req.dp_budget_pa, req.r_th_threshold_k_per_w, req.r_th_objective_k_per_w,
                                              req.marginal_return_cutoff, req.min_dp_ratio)
    order_r = sorted(range(len(rows)), key=lambda i: cands[i][1])
    rid = lambda i: None if i is None else rows[i]["design_id"]
    obj = rec["objective_step"]
    return {"r_key": r_key, "dp_key": dp_key, "recommended": rid(rec["index"]), "rule": rec["rule"],
            "meets_objective": rec["meets_objective"], "front": [rid(i) for i in rec["front"]],
            "walk_stop": rid(rec["walk_stop"]),
            "objective_step": None if obj is None else {"from": rid(obj["from"]), "to": rid(obj["to"])},
            "steps": [{"from": rid(s["from"]), "to": rid(s["to"]),
                       "returns": [{"to": rid(e["to"]), "return": e["return"], "skipped": e["skipped"]} for e in s["returns"]]}
                      for s in rec["steps"]],
            "order_by_r_th": [rid(i) for i in order_r]}


# --------------------------------------------------------------------------------------------------------------------
# Resolved developing solve (tall-fin limit)


def solve_developing(pf, req: Requirements, fluid, cells: float, inlet: str) -> dict:
    """Half channel over 0..60 mm (plus 5 mm and a short buffer before the outlet); inlet 'uniform' or 'parabolic'.

    'uniform': uniform velocity and temperature at x = -LEAD_IN_MM, slip (symmetry) walls upstream of the fin leading
    edge at x = 0, so the velocity jump sits at the leading edge rather than at a corner of the inflow boundary (where the
    pressure on the inflow line is singular and the inlet-to-outlet pressure drop does not converge with the mesh).
    'parabolic': developed velocity at x = 0 (thermal entrance only), no lead-in.
    """
    n = int(round(LENGTH_MM / WINDOW_MM))
    pf5 = replace(pf, period_mm=WINDOW_MM)
    lead = LEAD_IN_MM if inlet == "uniform" else 0.0
    pm = mesh_planform(pf5, n_periods=n + 1, periodic=False, cells=cells, lead_in_mm=lead)
    us = superficial_velocity(req, FLOW_LPM, pf.params["h_fin"])
    re_b = fluid.rho * us * B_M / fluid.mu
    pr = fluid.cp * fluid.mu / fluid.k
    u_mean = pm.flow_rate / pm.inlet_width  # mean channel velocity in units of U_s
    half = pf.params["w"] / 2
    profile = (lambda x: u_mean + 0 * x[0]) if inlet == "uniform" else \
        (lambda x: 1.5 * u_mean * (1 - ((x[1] - half) / half) ** 2))  # wall at y = 0, symmetry at y = w/2
    wx = pm.mesh.p[0, pm.mesh.facets[:, pm.walls]].mean(axis=0)
    walls, sym = pm.walls[wx > 0], np.concatenate([pm.symmetry, pm.walls[wx < 0]])
    fixed = {0: [(pm.inlet, profile), (walls, None)], 1: [(pm.inlet, None), (walls, None), (sym, None)]}
    steps = [r for r in (10.0, 50.0) if r < 0.8 * re_b] + [re_b]
    t0 = time.perf_counter()
    flow = solve_flow(pm.mesh, steps, fixed)
    energy = solve_energy(flow, pr, pm.inlet, pm.heated, pm.outlet)
    whole = window_means(flow, energy, pm.heated, 0.0, LENGTH_MM)
    p0, p1 = line_mean(pm.mesh, flow.basis_p, flow.p, -lead), line_mean(pm.mesh, flow.basis_p, flow.p, LENGTH_MM)
    windows = []
    for k in range(n):
        lo, hi = pm.cuts[k], pm.cuts[k + 1]
        wm = window_means(flow, energy, pm.heated, lo, hi)
        dpw = line_mean(pm.mesh, flow.basis_p, flow.p, lo) - line_mean(pm.mesh, flow.basis_p, flow.p, hi)
        windows.append({"x_mm": [float(lo), float(hi)], "h_w_m2k": wm["nusselt_b"] * fluid.k / B_M,
                        "dpdx_pa_per_m": dpw / (hi - lo) * fluid.rho * us**2 / B_M})
    return {"inlet": inlet, "cells": cells, "reynolds_b": re_b, "n_unknowns": int(flow.basis_u.N + flow.basis_p.N + energy.basis.N),
            "n_triangles": int(pm.mesh.t.shape[1]), "newton_iterations": flow.newton_iterations,
            "energy_balance_error": energy.energy_balance_error,
            "h_eff_w_m2k": whole["nusselt_b"] * fluid.k / B_M,  # q'' / (mean wall T - mean bulk T) over 0..60 mm
            "dp_pa": (p0 - p1) * fluid.rho * us**2, "windows": windows, "seconds": time.perf_counter() - t0}


def mesh_core_with_plenums(w: float, t: float, cells: float):
    """Plan-view strip y in [0, p/2] (fin centre at y = 0, channel centre at y = p/2, both symmetry planes): inlet plenum,
    half a fin of thickness t from x = 0 to 60 mm, outlet plenum. Lengths in mm. Returns the mesh, facet sets and the
    x positions of the embedded pressure lines in the outlet plenum."""
    import gmsh

    from .planforms import _to_meshtri
    half_p = (w + t) / 2
    x_in, x_out = -PLENUM_IN_MM, LENGTH_MM + PLENUM_OUT_MM
    stations = [LENGTH_MM + d for d in (2.5, 5.0, 7.5)]
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("plenums")
        occ = gmsh.model.occ
        dom = occ.addRectangle(x_in, 0, 0, x_out - x_in, half_p)
        fin = occ.addRectangle(0, 0, 0, LENGTH_MM, t / 2)
        out, _ = occ.cut([(2, dom)], [(2, fin)])
        lines = [(1, occ.addLine(occ.addPoint(x, -1, 0), occ.addPoint(x, half_p + 1, 0))) for x in stations]
        occ.fragment(out, lines)
        occ.synchronize()
        on_surface = set()
        for _, sfc in gmsh.model.getEntities(2):
            on_surface.update(abs(c) for c in gmsh.model.getAdjacencies(2, sfc)[1])
        stray = [(1, c) for _, c in gmsh.model.getEntities(1) if c not in on_surface]
        if stray:
            gmsh.model.removeEntities(stray, recursive=True)
        eps = 1e-6
        fin_curves = []
        for _, c in gmsh.model.getEntities(1):
            bb = gmsh.model.getBoundingBox(1, c)
            on_fin = (abs(bb[0]) < eps and abs(bb[3]) < eps) or (abs(bb[1] - t / 2) < eps and abs(bb[4] - t / 2) < eps) or \
                (abs(bb[0] - LENGTH_MM) < eps and abs(bb[3] - LENGTH_MM) < eps)
            if on_fin:
                fin_curves.append(c)
        fd = gmsh.model.mesh.field
        dist = fd.add("Distance")
        fd.setNumbers(dist, "CurvesList", fin_curves)
        fd.setNumber(dist, "Sampling", 4000)
        thr = fd.add("Threshold")
        fd.setNumber(thr, "InField", dist)
        fd.setNumber(thr, "SizeMin", w / (2 * cells))
        fd.setNumber(thr, "SizeMax", w / cells)
        fd.setNumber(thr, "DistMin", 0.0)
        fd.setNumber(thr, "DistMax", w / 2)
        fd.setAsBackgroundMesh(thr)
        for opt in ("Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature", "Mesh.MeshSizeExtendFromBoundary"):
            gmsh.option.setNumber(opt, 0)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.model.mesh.generate(2)
        mesh = _to_meshtri()
    finally:
        gmsh.finalize()
    tol = 1e-6
    bnd = mesh.boundary_facets()
    mx, my = mesh.p[:, mesh.facets[:, bnd]].mean(axis=1)
    inlet, outlet = bnd[np.abs(mx - x_in) < tol], bnd[np.abs(mx - x_out) < tol]
    sym = bnd[(np.abs(my - half_p) < tol) | ((np.abs(my) < tol) & ((mx < 0) | (mx > LENGTH_MM)))]
    walls = np.setdiff1d(bnd, np.concatenate([inlet, outlet, sym]))
    return mesh, inlet, outlet, sym, walls, stations


def solve_core_with_plenums(pf, req: Requirements, fluid, cells: float) -> dict:
    """Pressure drop from the inlet plenum to the outlet plenum through the 60 mm fin array (tall-fin limit, flow only).

    Uniform inflow at the superficial velocity 2 mm upstream of the square fin leading edge; pressure differences to
    lines 2.5, 5 and 7.5 mm after the array and to the traction-free outlet 10 mm after it.
    """
    w, t = pf.params["w"], pf.params["t"]
    mesh, inlet, outlet, sym, walls, stations = mesh_core_with_plenums(w, t, cells)
    us = superficial_velocity(req, FLOW_LPM, pf.params["h_fin"])
    re_b = fluid.rho * us * B_M / fluid.mu
    fixed = {0: [(inlet, lambda x: 1.0 + 0 * x[0]), (walls, None)], 1: [(inlet, None), (walls, None), (sym, None)]}
    steps = [r for r in (10.0, 50.0) if r < 0.8 * re_b] + [re_b]
    t0 = time.perf_counter()
    flow = solve_flow(mesh, steps, fixed)
    p_in = line_mean(mesh, flow.basis_p, flow.p, -PLENUM_IN_MM)
    scale = fluid.rho * us**2
    dp = {f"{x - LENGTH_MM:.1f}": (p_in - line_mean(mesh, flow.basis_p, flow.p, x)) * scale for x in stations}
    dp[f"{PLENUM_OUT_MM:.1f}"] = (p_in - line_mean(mesh, flow.basis_p, flow.p, LENGTH_MM + PLENUM_OUT_MM)) * scale
    return {"cells": cells, "n_triangles": int(mesh.t.shape[1]), "n_unknowns": int(flow.basis_u.N + flow.basis_p.N),
            "dp_pa_to_mm_after_array": dp, "dp_pa": dp[f"{PLENUM_OUT_MM:.1f}"], "seconds": time.perf_counter() - t0}


def one_d_parallel_plates(w_mm: float, u_channel: float, dp_fd_pa: float, fluid) -> dict:
    """The 1D option's factors for a parallel-plate channel of gap w (both walls heated) at channel velocity u.

    Re-anchored model (the screen): Part II entrance asymptotes blended with the screen's developed Nu in the plate limit,
    Nu_H1(0) = 8.235 on Dh = 2w, and an incremental pressure drop K(inf) rho u^2 / 2. The earlier formulation (Part II
    blended with its own Eq. 31, Part I apparent friction) is kept for comparison only."""
    w = w_mm * 1e-3
    sqrt_a = math.sqrt(w * w / PLATE_EPS)
    pr = fluid.cp * fluid.mu / fluid.k
    l_plus = LENGTH_MM * 1e-3 * fluid.mu / (fluid.rho * u_channel * sqrt_a**2)
    nu_anchor = nu_h1(0.0) * sqrt_a / (2 * w)  # h = 8.235 k / (2w), expressed on sqrt(A)
    q = 0.5 * fluid.rho * u_channel**2
    return {"nu_fd_anchor_sqrt_a": nu_anchor, "nu_fd_eq31_over_anchor": nu_fd_sqrt_a(PLATE_EPS) / nu_anchor,
            "nu_factor_thermal": effective_nu_factor(PLATE_EPS, l_plus / pr, nu_fd=nu_anchor),
            "nu_factor_combined": effective_nu_factor(PLATE_EPS, l_plus / pr, pr, nu_fd=nu_anchor),
            "k_incremental": K_INF_PLATES, "dp_increase": K_INF_PLATES * q / dp_fd_pa,
            "earlier_nu_factor_thermal": effective_nu_factor(PLATE_EPS, l_plus / pr),
            "earlier_nu_factor_combined": effective_nu_factor(PLATE_EPS, l_plus / pr, pr),
            "earlier_dp_increase": f_app_re_sqrt_a(PLATE_EPS, l_plus) / f_re_sqrt_a(PLATE_EPS) - 1,
            "l_plus_sqrt_a": l_plus, "x_star_2w": LENGTH_MM * 1e-3 * fluid.mu / (fluid.rho * u_channel * 2 * w * 2 * w) / pr}


def resolved_case(row: dict, req: Requirements, fluid, mesh_check: bool) -> dict:
    pf = straight(row["w_mm"], row["t_mm"], row["h_mm"])
    per = run_design(pf, req, fluid, FLOW_LPM, cells=CELLS_PERIODIC)
    u_ch = superficial_velocity(req, FLOW_LPM, row["h_mm"]) / pf.porosity
    q = 0.5 * fluid.rho * u_ch**2
    out = {"design_id": row["design_id"], "periodic": {"dpdx_pa_per_m": per.dpdx_pa_per_m, "h_w_m2k": per.h_w_m2k,
                                                         "dp_pa": per.dp_pa, "r_th_k_per_w": per.r_th_k_per_w,
                                                         "energy_balance_error": per.energy_balance_error},
           "u_channel_m_s": u_ch, "reynolds_2w": fluid.rho * u_ch * 2 * row["w_mm"] * 1e-3 / fluid.mu,
           "one_d_plates": one_d_parallel_plates(row["w_mm"], u_ch, per.dp_pa, fluid)}
    f1 = out["one_d_plates"]
    levels = (CELLS, CELLS_FINE) if mesh_check else (CELLS,)
    for inlet in ("uniform", "parabolic"):
        sols = [solve_developing(pf, req, fluid, c, inlet) for c in levels]
        out[inlet] = sols[0]  # production mesh; the refinement only checks it
        if mesh_check:
            out[inlet]["mesh_change"] = {
                "h_eff": sols[1]["h_eff_w_m2k"] / sols[0]["h_eff_w_m2k"] - 1, "dp": sols[1]["dp_pa"] / sols[0]["dp_pa"] - 1,
                "fine_cells": CELLS_FINE, "fine_h_eff_w_m2k": sols[1]["h_eff_w_m2k"], "fine_dp_pa": sols[1]["dp_pa"],
                "h_increase_fine": sols[1]["h_eff_w_m2k"] / per.h_w_m2k - 1, "dp_increase_fine": sols[1]["dp_pa"] / per.dp_pa - 1,
                "fine_energy_balance_error": sols[1]["energy_balance_error"]}
    plen = [solve_core_with_plenums(pf, req, fluid, c) for c in levels]
    out["plenums"] = {**plen[0], "excess_over_developed_q": (plen[0]["dp_pa"] - per.dp_pa) / q,
                      "inlet_outlet_q": (plen[0]["dp_pa"] - out["uniform"]["dp_pa"]) / q,
                      "developing_excess_q": (out["uniform"]["dp_pa"] - per.dp_pa) / q,
                      "one_d_developing_excess_q": f1["k_incremental"],
                      "earlier_one_d_developing_excess_q": f1["earlier_dp_increase"] * per.dp_pa / q,
                      "one_d_kc_plus_ke": k_contraction(pf.porosity, 0.0) + k_expansion(pf.porosity, 0.0),
                      "sigma": pf.porosity, "dynamic_pressure_pa": q}
    if mesh_check:
        out["plenums"]["mesh_change"] = {"dp": plen[1]["dp_pa"] / plen[0]["dp_pa"] - 1, "fine_dp_pa": plen[1]["dp_pa"]}
    # plan-view network with developed, resolved and 1D-corrected values
    net = lambda dp_factor, h_factor: network(pf, per.dpdx_pa_per_m * dp_factor, per.h_w_m2k * h_factor, req, fluid, FLOW_LPM)
    fd = net(1.0, 1.0)
    dev = network(pf, out["uniform"]["dp_pa"] / (LENGTH_MM * 1e-3), out["uniform"]["h_eff_w_m2k"], req, fluid, FLOW_LPM)
    th = network(pf, per.dpdx_pa_per_m, out["parabolic"]["h_eff_w_m2k"], req, fluid, FLOW_LPM)
    one = net(1 + f1["dp_increase"], f1["nu_factor_combined"])
    one_th = net(1.0, f1["nu_factor_thermal"])
    early = net(1 + f1["earlier_dp_increase"], f1["earlier_nu_factor_combined"])
    out["network"] = {k: {"r_th_k_per_w": v["r_th_k_per_w"], "dp_pa": v["dp_pa"]} for k, v in
                      (("fully_developed", fd), ("resolved_combined", dev), ("resolved_thermal", th), ("one_d_combined", one),
                       ("one_d_thermal", one_th), ("earlier_one_d_combined", early))}
    res_dp = out["uniform"]["dp_pa"] / per.dp_pa - 1
    d_r_res = dev["r_th_k_per_w"] - fd["r_th_k_per_w"]
    d_r_1d, d_r_early = one["r_th_k_per_w"] - fd["r_th_k_per_w"], early["r_th_k_per_w"] - fd["r_th_k_per_w"]
    d_r_res_th, d_r_1d_th = th["r_th_k_per_w"] - fd["r_th_k_per_w"], one_th["r_th_k_per_w"] - fd["r_th_k_per_w"]
    out["comparison"] = {
        "h_increase_resolved": out["uniform"]["h_eff_w_m2k"] / per.h_w_m2k - 1, "h_increase_1d": f1["nu_factor_combined"] - 1,
        "h_increase_resolved_thermal": out["parabolic"]["h_eff_w_m2k"] / per.h_w_m2k - 1,
        "h_increase_1d_thermal": f1["nu_factor_thermal"] - 1,
        "dp_increase_resolved": res_dp, "dp_increase_1d": f1["dp_increase"],
        "r_th_change_resolved": d_r_res, "r_th_change_1d": d_r_1d,
        "r_th_change_resolved_thermal": d_r_res_th, "r_th_change_1d_thermal": d_r_1d_th,
        "r_th_contribution_disagreement": d_r_1d / d_r_res - 1,
        "r_th_thermal_contribution_disagreement": d_r_1d_th / d_r_res_th - 1,
        "dp_contribution_disagreement": f1["dp_increase"] / res_dp - 1,
        "earlier_h_increase_1d": f1["earlier_nu_factor_combined"] - 1, "earlier_dp_increase_1d": f1["earlier_dp_increase"],
        "earlier_r_th_contribution_disagreement": d_r_early / d_r_res - 1,
        "earlier_dp_contribution_disagreement": f1["earlier_dp_increase"] / res_dp - 1}
    return out


# --------------------------------------------------------------------------------------------------------------------
# Figure and summary


def plot_front(rows: list[dict], req: Requirements, recommended: str, out: Path) -> Path:
    from tools.report import draw_requirements
    fig, ax = ps.figure(5.4)
    budget_kpa = req.dp_budget_pa / 1e3
    fd = [(r["cfd_anchored"]["dp_fd"] / 1e3, r["cfd_anchored"]["r_th_fd"] * 1e3) for r in rows]
    comb = [(r["cfd_anchored"]["dp_combined_core"] / 1e3, r["cfd_anchored"]["r_th_combined"] * 1e3) for r in rows]
    y_lo = 0.7 * min(p[1] for p in comb + fd)
    y_hi = 1.5 * req.r_th_threshold_k_per_w * 1e3
    x_lo, x_hi = 0.3, 1.5 * budget_kpa
    draw_requirements(ax, x_lo, y_lo, budget_kpa, req.r_th_threshold_k_per_w * 1e3, req.r_th_objective_k_per_w * 1e3, y_hi)
    order = np.argsort([p[0] for p in fd])
    ax.plot([fd[i][0] for i in order], [fd[i][1] for i in order], color=ps.BLUE, lw=1.2, zorder=3)
    ax.scatter(*zip(*fd), s=50, color=ps.BLUE, edgecolors="white", linewidths=0.8, zorder=4,
               label="Fully developed (Level 1 CFD)")
    ax.plot([comb[i][0] for i in order], [comb[i][1] for i in order], color=ps.COPPER, lw=1.2, zorder=3)
    ax.scatter(*zip(*comb), s=50, color=ps.COPPER, edgecolors="white", linewidths=0.8, zorder=4,
               label="With developing flow, channel core")
    for a, b in zip(fd, comb):
        ax.annotate("", xy=b, xytext=a, arrowprops={"arrowstyle": "->", "color": ps.MUTED, "lw": 0.9}, zorder=2)
    i = next(k for k, r in enumerate(rows) if r["design_id"] == recommended)
    for p in (fd[i], comb[i]):
        ax.scatter([p[0]], [p[1]], s=220, facecolors="none", edgecolors=ps.ACCENT, linewidths=2.0, zorder=5)
    ax.scatter([], [], s=120, facecolors="none", edgecolors=ps.ACCENT, linewidths=2.0, label="Recommended design")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ps.plain_log_ticks(ax, x_lo, x_hi, "x")
    ps.plain_log_ticks(ax, y_lo, y_hi, "y")
    ax.set_xlabel("ΔP (kPa)")
    ax.set_ylabel("R_th (mK/W)")
    ax.set_title("Straight-fin front with and without entrance effects", loc="left")
    ax.legend(loc="upper right", bbox_to_anchor=(0.97, 0.88))
    ps.save(fig, out)
    return out


def geo(r: dict) -> str:
    return f"{r['w_mm']:.2f} / {r['t_mm']:.2f} / {r['h_mm']:.2f}"


def write_summary(data: dict, out: Path) -> Path:
    rows, res, req = data["one_d"], data["resolved"], data["requirements"]
    rec_id = data["level1_recommended"]
    rec = next(r for r in rows if r["design_id"] == rec_id)
    rres = next(r for r in res if r["design_id"] == rec_id)
    obj, thr = req["r_th_objective_k_per_w"], req["r_th_threshold_k_per_w"]
    ca, alt = rec["cfd_anchored"], rres["cfd_anchored_with_resolved_factors"]
    margin = lambda r: 100 * (1 - r / obj)
    sel_fd, sel_core = data["selection_fully_developed"], data["selection_core"]
    newly = [r for r in rows if r["cfd_anchored"]["r_th_fd"] > thr >= r["cfd_anchored"]["r_th_combined"]]
    by_id = {r["design_id"]: r for r in rows}
    walk = lambda s: " → ".join([geo(by_id[s["steps"][0]["from"]])] + [geo(by_id[x["to"]]) for x in s["steps"] if x["to"]])
    g = data["gates"]
    c = rres["comparison"]
    gate_r = next(v for k, v in g.items() if "R_th" in k)
    gate_dp = next(v for k, v in g.items() if "ΔP" in k)
    passed = gate_r["pass"] and gate_dp["pass"]
    k_inf = rec["k_incremental"]
    L = []
    a = L.append
    a("# Entrance effects in the straight-fin study (issue #12)")
    a("")
    a(f"**Answer.** For the recommended design ({rec['w_mm']:.2f} mm channels, {rec['t_mm']:.2f} mm fins, {rec['h_mm']:.2f} mm "
      f"tall, {rec['n_channels']} channels, Re {rec['reynolds']:.0f}, x* = L/(Dh·Re·Pr) = {rec['x_star']:.3f}), developing flow "
      f"lowers R_th from {ca['r_th_fd']*1e3:.2f} to {ca['r_th_combined']*1e3:.2f} mK/W "
      f"({100*(ca['r_th_combined']/ca['r_th_fd']-1):+.1f}%) and raises the channel-core ΔP from {ca['dp_fd']/1e3:.3f} to "
      f"{ca['dp_combined_core']/1e3:.3f} kPa ({100*(ca['dp_combined_core']/ca['dp_fd']-1):+.1f}%). "
      f"**The 15 mK/W objective is met with a margin of {margin(ca['r_th_combined']):.1f}%** "
      f"({margin(ca['r_th_fd']):.1f}% fully developed). Values are the Level 1 CFD result plus the 1D increment. "
      f"With the resolved tall-fin ratios in place of the 1D ones the design reads {alt['r_th_combined']*1e3:.2f} mK/W at "
      f"{alt['dp_combined_core']/1e3:.3f} kPa; plenum to plenum the resolved solve adds a further "
      f"{rres['plenums']['inlet_outlet_q']:.2f}·ρU²/2 ({rres['plenums']['inlet_outlet_q']*rec['dynamic_pressure_pa']:.0f} Pa) "
      "for the inlet and outlet.")
    a("")
    a(f"**Gate.** In the tall-fin (parallel-plate) limit the re-anchored 1D correction and the resolved solve agree on the "
      f"entrance contribution to within {gate_r['value']} on R_th and {gate_dp['value']} on ΔP for the recommended design "
      f"(limit ±10%): **{'pass' if passed else 'FAIL'}**. All five resolved designs are in the table below. The first version of "
      f"this check (earlier Level 1 pick, 0.42 / 0.44 / 4.0 mm) failed at −36% on R_th and +41% on ΔP; the earlier "
      f"formulation, applied to the current pick, still gives {100*c['earlier_r_th_contribution_disagreement']:+.0f}% and "
      f"{100*c['earlier_dp_contribution_disagreement']:+.0f}%. "
      + "Across the five designs the re-anchored correction differs from the resolved contribution by " + ", ".join(
          f"{100*x['comparison']['r_th_contribution_disagreement']:+.1f}%" for x in res) + " on R_th and " + ", ".join(
          f"{100*x['comparison']['dp_contribution_disagreement']:+.1f}%" for x in res) + " on ΔP (widest to narrowest); "
      "the widest channel is the only one outside ±10% on R_th, where x* is smallest and the entrance dominates. "
      "What changed, approved by Riccardo before the rerun and not fitted to the resolved solves:")
    a("- *Thermal:* the published entrance asymptotes are now blended with the screen's own developed Nusselt number (Shah & "
      "London Nu_H1) instead of the source's developed-limit formula, which at small aspect ratio reads "
      f"{100*(rres['one_d_plates']['nu_fd_eq31_over_anchor']-1):.0f}% above the exact parallel-plate value and so hid part of "
      "the entrance gain. The long-channel limit is now exactly the fully developed model.")
    a(f"- *Hydrodynamic:* fully developed friction plus an incremental pressure drop K(∞)·ρU²/2 with K(∞) = {k_inf:.3f} "
      "(parallel plates, from a published table; see Sources), in place of the apparent-friction blend, whose long-channel "
      f"tail was worth about {rres['plenums']['earlier_one_d_developing_excess_q']:.2f}·ρU²/2.")
    a(f"- *Inlet and outlet:* the Kays coefficients are reported as a separate term and no longer added to any combined "
      f"figure (reason below). For the recommended design they would add {rec['one_d']['dp_inlet_outlet']:.0f} Pa.")
    a("")
    a(f"**Front and recommendation.** The order of the {len(rows)} CFD-verified designs by R_th is "
      f"{'unchanged' if data['ranking_unchanged'] else 'changed'}; wide channels gain most (x* is smallest there). " + (
          f"{', '.join(geo(r) + ' mm' for r in newly)} now meet{'s' if len(newly) == 1 else ''} the 20 mK/W threshold "
          f"({', '.join(f'{r['cfd_anchored']['r_th_combined']*1e3:.2f}' for r in newly)} mK/W). " if newly else "") +
      f"The selection rule on the corrected values walks {walk(sel_core)} mm and picks "
      f"{'the same design' if sel_core['recommended'] == rec_id else geo(by_id[sel_core['recommended']]) + ' mm'} "
      f"(fully developed: {walk(sel_fd)} mm). "
      + ("The recommendation would not change." if sel_core["recommended"] == rec_id else change_text(sel_core, sel_fd, rec_id, by_id, req)))
    a("")
    a("Geometry is channel width / fin thickness / fin height in mm.")
    a("")
    a("## Sources and formulas")
    a("See `docs/literature.md` and `ASSUMPTIONS.md` rows P2 and P9–P12.")
    a("- **Thermal and simultaneous entrance, uniform wall flux:** Muzychka & Yovanovich (2002), Part II, Eqs. 50–52 with their "
      "Table 2 (C3 = 0.501, f(Pr) = 0.886/[1 + (1.909 Pr^(1/6))^(9/2)]^(2/9), m = 2.27 + 1.65 Pr^(1/3)), local form "
      "(C2 = C4 = 1), z* = z/(√A·Re_√A·Pr): Nu_G = [(C3 (f·Re_√A / z*)^(1/3))⁵ + Nu_fd⁵]^(1/5), "
      "Nu = [Nu_G^m + (f(Pr)/√z*)^m]^(1/m), with Nu_fd the screen's Shah & London Nu_H1 converted to √A (re-anchored). Valid "
      "for 0.1 < Pr < ∞; stated accuracy ±12% (Graetz) and ±15% (combined entry). R_th uses the mean of (wall − bulk) "
      "temperature over the length, so the effective h is the harmonic length-mean of the local h.")
    a(f"- **Hydrodynamic entrance:** ΔP = ΔP_fd + K(∞)·ρU²/2, U the channel mean velocity, K(∞) = {k_inf:.3f}: Ray, Durst & "
      "Ray (2017), arXiv:1707.04947, Table 2, parallel-plate channel, Kn = 10⁻⁴ (continuum), first-order slip (C2 = 0): "
      "K1 = 0.6802, the high-Re asymptote of their Navier–Stokes K_fd from a uniform inlet. Their full correlation adds K0/Re "
      "(K0 = 80.3 on Dh = 2 × gap), which comes from axial diffusion upstream of a uniform profile imposed at the inlet plane; "
      "in a plate fed from a plenum that profile is not imposed, and the resolved solve (with a slip lead-in ahead of the fin) "
      "gives a Re-independent value close to K1. The same K(∞) is used for the real duct (aspect ratio about 0.1); no open "
      "tabulation for rectangular ducts was obtained.")
    a("- **Inlet and outlet losses (reported separately):** laminar abrupt-contraction and abrupt-expansion coefficients of "
      "Kays (1950) for flat and square ducts, as polynomial fits in σ interpolated linearly in aspect ratio by dos Reis & "
      "Altemani (COBEM 1999, Eq. 5, Table 1); ΔP_io = (Kc + Ke)·ρU²/2, σ = w/(w + t), plenum as tall as the fins.")
    a("")
    a("## R_th and ΔP per design")
    a("")
    a("| Design | x* | R_th fully developed | thermal entrance only | combined entry | ΔP fully developed | ΔP with developing flow (core) | inlet/outlet, separate |")
    a("|---|---|---|---|---|---|---|---|")
    for r in rows:
        cc = r["cfd_anchored"]
        tag = " **(recommended)**" if r["design_id"] == rec_id else ""
        a(f"| {geo(r)}{tag} | {r['x_star']:.3f} | {cc['r_th_fd']*1e3:.2f} | {cc['r_th_thermal']*1e3:.2f} | "
          f"{cc['r_th_combined']*1e3:.2f} | {cc['dp_fd']:.0f} | {cc['dp_combined_core']:.0f} | {cc['dp_inlet_outlet_separate']:.0f} |")
    a("")
    a("R_th in mK/W, ΔP in Pa. Thermal entrance only: Graetz model (developed velocity). Combined entry: velocity and "
      "temperature developing together from a uniform inlet; this is the corrected R_th. The inlet/outlet column is "
      "(Kc + Ke)·ρU²/2 and is not included in any other column. Coefficients and factors:")
    a("")
    a("| Design | Re | σ | Kc | Ke | h factor, thermal | h factor, combined | ρU²/2 (Pa) | ΔP factor | L/(Dh·Re) |")
    a("|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        a(f"| {geo(r)} | {r['reynolds']:.0f} | {r['sigma']:.3f} | {r['k_contraction']:.3f} | {r['k_expansion']:.3f} | "
          f"{r['nu_factor_thermal']:.3f} | {r['nu_factor_combined']:.3f} | {r['dynamic_pressure_pa']:.1f} | {r['dp_factor']:.4f} | "
          f"{r['l_plus_dh']:.3f} |")
    a("")
    a("![Straight-fin front with and without entrance effects](front_entrance.png)")
    a("")
    off = [geo(r) for r in rows if r["cfd_anchored"]["dp_fd"] < 300.0]
    a(f"*R_th against channel-core ΔP for the {len(rows)} CFD-verified straight-fin designs of the Level 1 study at 1.5 LPM "
      "PG25: fully developed (Level 1 CFD, blue) and with developing flow (copper, arrows join each design); corrections are "
      "the re-anchored 1D increments added to the CFD values; inlet and outlet losses are not included. Budget 5 kPa "
      f"(channel core), threshold 20 mK/W, objective 15 mK/W. Circled: the Level 1 recommended design ({geo(rec)} mm). "
      "Region of interest as in the Level 1 report (from 0.3 kPa)" +
      (f": {', '.join(off)} mm (below 0.3 kPa) {'is' if len(off) == 1 else 'are'} off the plot and in the tables.*" if off else ".*"))
    a("")
    a("## Resolved check (tall-fin limit)")
    a("")
    a("Three solves per design with the plan-view solver (`level2/planar.py`, the solver behind "
      "`level2.planview.run_developing`), uniform flux on the fin walls, inlet temperature uniform:")
    a(f"1. **Combined entry:** one straight half channel over the full 60 mm, uniform velocity {LEAD_IN_MM:g} mm upstream "
      "of the fin leading edge with slip walls there (a uniform profile imposed right at the wall corner makes the inlet "
      "pressure singular, and the pressure drop then grows with mesh refinement; the slip lead-in removes that).")
    a("2. **Thermal entrance only:** the same channel with a parabolic inlet profile at x = 0.")
    a(f"3. **Plenum to plenum (ΔP only):** a {PLENUM_IN_MM:g} mm inlet plenum, half a fin of thickness t with square ends "
      f"over 60 mm, a {PLENUM_OUT_MM:g} mm outlet plenum; plenums as tall as the fins (σ = w/(w + t)).")
    a("")
    a("Length-mean h = q''/(mean wall T − mean bulk T over 0–60 mm); ΔP from the inlet line to x = 60 mm. Developed "
      "reference: the streamwise-periodic cell. The 1D columns evaluate the re-anchored correction for the same channel "
      "(parallel plates of gap w, aspect ratio 0.01, developed Nu 8.235 on 2w); the earlier columns are the first formulation.")
    a("")
    a("**Tall-fin limit.** The plan view has fins of infinite height: no floor, no cover, both channel walls heated. So it "
      "checks development between two fins, not the real duct with an aspect ratio of about 0.1 and an adiabatic cover "
      "(ASSUMPTIONS V1). Corner regions of the real duct develop differently and are not checked here.")
    a("")
    a("| Design | Re (2w) | x* (2w) | h increase: resolved / 1D / earlier 1D | ΔP increase: resolved / 1D / earlier 1D | "
      "ΔR_th resolved / 1D (mK/W) | 1D vs resolved, R_th contribution | 1D vs resolved, ΔP contribution | energy balance |")
    a("|---|---|---|---|---|---|---|---|---|")
    for r in res:
        cc = r["comparison"]
        tag = " **(recommended)**" if r["design_id"] == rec_id else ""
        a(f"| {geo(by_id[r['design_id']])}{tag} | {r['reynolds_2w']:.0f} | {r['one_d_plates']['x_star_2w']:.3f} | "
          f"{100*cc['h_increase_resolved']:+.1f}% / {100*cc['h_increase_1d']:+.1f}% / {100*cc['earlier_h_increase_1d']:+.1f}% | "
          f"{100*cc['dp_increase_resolved']:+.2f}% / {100*cc['dp_increase_1d']:+.2f}% / {100*cc['earlier_dp_increase_1d']:+.2f}% | "
          f"{1e3*cc['r_th_change_resolved']:+.3f} / {1e3*cc['r_th_change_1d']:+.3f} | "
          f"{100*cc['r_th_contribution_disagreement']:+.1f}% | {100*cc['dp_contribution_disagreement']:+.1f}% | "
          f"{100*r['uniform']['energy_balance_error']:+.3f}% |")
    a("")
    a("ΔR_th: change of the plan-view network R_th from the developed value. Thermal entrance only (parabolic inlet):")
    a("")
    a("| Design | h increase, resolved | h increase, 1D (Graetz) | ΔR_th resolved (mK/W) | ΔR_th 1D (mK/W) | energy balance |")
    a("|---|---|---|---|---|---|")
    for r in res:
        cc = r["comparison"]
        a(f"| {geo(by_id[r['design_id']])} | {100*cc['h_increase_resolved_thermal']:+.1f}% | {100*cc['h_increase_1d_thermal']:+.1f}% | "
          f"{1e3*cc['r_th_change_resolved_thermal']:+.3f} | {1e3*cc['r_th_change_1d_thermal']:+.3f} | "
          f"{100*r['parabolic']['energy_balance_error']:+.3f}% |")
    a("")
    a("Plenum to plenum, in units of the channel dynamic pressure ρU²/2:")
    a("")
    a("| Design | σ | ρU²/2 (Pa) | total excess over developed, resolved | of which developing channel (K, resolved) | "
      "of which inlet + outlet, resolved | K(∞) used in 1D | Kc + Ke (Kays, flat duct) |")
    a("|---|---|---|---|---|---|---|---|")
    for r in res:
        p = r["plenums"]
        a(f"| {geo(by_id[r['design_id']])} | {p['sigma']:.3f} | {p['dynamic_pressure_pa']:.1f} | {p['excess_over_developed_q']:.2f} | "
          f"{p['developing_excess_q']:.2f} | {p['inlet_outlet_q']:+.2f} | {p['one_d_developing_excess_q']:.2f} | "
          f"{p['one_d_kc_plus_ke']:.2f} |")
    a("")
    mc, mcp, mpl = rres["uniform"]["mesh_change"], rres["parabolic"]["mesh_change"], rres["plenums"]["mesh_change"]
    a(f"Mesh check, recommended design, {CELLS:g} → {CELLS_FINE:g} cells across the half channel: combined entry h "
      f"{100*mc['h_eff']:+.3f}%, ΔP {100*mc['dp']:+.4f}% (entrance increments {100*c['h_increase_resolved']:.2f}% → "
      f"{100*mc['h_increase_fine']:.2f}% on h and {100*c['dp_increase_resolved']:.3f}% → "
      f"{100*mc['dp_increase_fine']:.3f}% on ΔP); thermal entrance h {100*mcp['h_eff']:+.4f}%; plenum ΔP {100*mpl['dp']:+.4f}%.")
    a("")
    a("**Gates.** " + "; ".join(f"{k}: {v['value']} (limit {v['limit']}), {'pass' if v['pass'] else '**FAIL**'}"
                                for k, v in g.items()) + ".")
    a("")
    ks = [x["plenums"]["developing_excess_q"] for x in res]
    io = [x["plenums"]["inlet_outlet_q"] for x in res]
    kk = [x["plenums"]["one_d_kc_plus_ke"] for x in res]
    offs = [b - a_ for a_, b in zip(io, kk)]
    a("## Why inlet and outlet losses are kept separate")
    a(f"Plenum to plenum, the resolved inlet-plus-outlet loss is {min(io):+.2f} to {max(io):+.2f}·ρU²/2 over the five designs, "
      f"against Kays Kc + Ke of {min(kk):.2f} to {max(kk):.2f}. The gap is nearly constant ({min(offs):.2f} to {max(offs):.2f}) "
      f"and close to the resolved entrance excess of the channel itself ({min(ks):.2f} to {max(ks):.2f}). That pattern fits "
      "laminar contraction coefficients that already contain part of the entrance pressure defect, so adding them to K(∞) "
      "would count it twice. Kays's paper was not obtained, so this could not be confirmed from the source. The term is "
      "also outside the channel-core ΔP definition of the study. It is therefore listed on its own and added to nothing.")
    a("")
    a("## Not covered")
    a("- Flow maldistribution between channels, manifold, port and fitting losses (ASSUMPTIONS P5, P6).")
    a("- Axial conduction in the copper base, which moves heat toward the inlet where h is higher. With it the effective h lies "
      "between the harmonic length-mean used here and the arithmetic mean, so the R_th reductions here are the smaller "
      "estimate.")
    a("- Temperature-dependent viscosity (P3), and the 3D development of the real duct (floor, cover, corners); the K(∞) "
      "for the real duct is the parallel-plate value.")
    a("")
    a("Reproduce: `uv run python -m level2.entrance` (about 5 minutes); `--plot-only` redraws the figure and this file from "
      "`results.json`.")
    out.write_text("\n".join(L) + "\n")
    return out


def change_text(sel_core: dict, sel_fd: dict, rec_id: str, by_id: dict, req: dict) -> str:
    """Why the selection rule lands elsewhere on the corrected values: the return of the step it no longer takes."""
    new_id = sel_core["recommended"]
    ret = lambda sel, a, b: next((e["return"] for s in sel["steps"] if s["from"] == a for e in s["returns"] if e["to"] == b),
                                 None)
    r_new, r_fd = ret(sel_core, new_id, rec_id), ret(sel_fd, new_id, rec_id)
    cc = by_id[new_id]["cfd_anchored"]
    txt = (f"**The pick would move** to {geo(by_id[new_id])} mm ({cc['r_th_combined']*1e3:.2f} mK/W at "
           f"{cc['dp_combined_core']/1e3:.3f} kPa with entrance effects, {cc['r_th_fd']*1e3:.2f} mK/W at {cc['dp_fd']/1e3:.3f} kPa "
           "fully developed). ")
    if r_new is not None and r_fd is not None:
        txt += (f"The step from it to the current pick returns {r_new:.3f} with entrance effects against {r_fd:.3f} fully "
                f"developed, just under the {req['marginal_return_cutoff']:.2f} cut-off: the corrections lower R_th by "
                "slightly more on the wider design. ")
    return txt + ("This is the same one-step move the cut-off sensitivity in the Level 1 report shows near 0.25; the study keeps "
                  "fully developed values by default, so it is reported here and not applied.")


def gates_for(res: list[dict], rec_id: str) -> dict:
    r = next(x for x in res if x["design_id"] == rec_id)
    c, mc = r["comparison"], r["uniform"]["mesh_change"]
    mcp = r["parabolic"]["mesh_change"]
    worst_mesh = max(abs(mc["dp"]), abs(mc["h_eff"]), abs(mcp["dp"]), abs(mcp["h_eff"]), abs(r["plenums"]["mesh_change"]["dp"]))
    worst_energy = max([abs(x[k]["energy_balance_error"]) for x in res for k in ("uniform", "parabolic")] +
                       [abs(r[k]["mesh_change"]["fine_energy_balance_error"]) for k in ("uniform", "parabolic")])
    gates = {
        "mesh change at 2x refinement, recommended design": {"value": f"{100*worst_mesh:.2f}%", "limit": "2%",
                                                             "pass": worst_mesh < GATE_MESH},
        "energy balance, all resolved solves": {"value": f"{100*worst_energy:.3f}%", "limit": "1%", "pass": worst_energy < GATE_ENERGY},
        "1D vs resolved, entrance contribution to R_th (recommended)": {
            "value": f"{100*c['r_th_contribution_disagreement']:+.1f}%", "limit": "±10%",
            "pass": abs(c["r_th_contribution_disagreement"]) <= GATE_AGREEMENT},
        "1D vs resolved, entrance contribution to ΔP (recommended)": {
            "value": f"{100*c['dp_contribution_disagreement']:+.1f}%", "limit": "±10%",
            "pass": abs(c["dp_contribution_disagreement"]) <= GATE_AGREEMENT}}
    return gates


def run(plot_only: bool = False, level1_path: Path = LEVEL1, out_dir: Path = OUT) -> int:
    level1 = json.loads(Path(level1_path).read_text())
    req = requirements(level1)
    fluid = fluid_props(req.coolant, req.inlet_temp_c)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "results.json"
    if plot_only:
        data = json.loads(path.read_text())
    else:
        t0 = time.perf_counter()
        rows = one_d_table(level1, req, fluid)
        res = []
        for did in resolved_designs(rows, level1["recommended"]):
            row = next(r for r in rows if r["design_id"] == did)
            res.append(resolved_case(row, req, fluid, mesh_check=did == level1["recommended"]))
            c = res[-1]["comparison"]
            ca = row["cfd_anchored"]
            res[-1]["cfd_anchored_with_resolved_factors"] = {  # the tall-fin resolved ratios applied to the real duct instead
                "r_th_combined": ca["r_th_fd"] + r_th_with_h_factor(row, req, fluid, 1 + c["h_increase_resolved"])
                - row["one_d"]["r_th_fd"],
                "dp_combined_core": ca["dp_fd"] + res[-1]["plenums"]["developing_excess_q"] * row["dynamic_pressure_pa"]}
            print(f"{did}: h {100*c['h_increase_resolved']:+.2f}% (1D {100*c['h_increase_1d']:+.2f}%), "
                  f"dP {100*c['dp_increase_resolved']:+.2f}% (1D {100*c['dp_increase_1d']:+.2f}%), "
                  f"R_th contribution {100*c['r_th_contribution_disagreement']:+.1f}%, "
                  f"dP contribution {100*c['dp_contribution_disagreement']:+.1f}%, "
                  f"energy {100*res[-1]['uniform']['energy_balance_error']:+.3f}%", flush=True)
        sel_fd = selection(rows, req, "r_th_fd", "dp_fd")
        sel_core = selection(rows, req, "r_th_combined", "dp_combined_core")
        data = {"level1_recommended": level1["recommended"], "requirements": asdict(req), "fluid": asdict(fluid),
                "settings": {"flow_lpm": FLOW_LPM, "length_mm": LENGTH_MM, "window_mm": WINDOW_MM, "cells": CELLS,
                             "cells_fine": CELLS_FINE, "cells_periodic": CELLS_PERIODIC, "plate_aspect_ratio": PLATE_EPS,
                             "lead_in_mm": LEAD_IN_MM, "plenum_in_mm": PLENUM_IN_MM, "plenum_out_mm": PLENUM_OUT_MM,
                             "k_incremental": K_INF_PLATES},
                "one_d": rows, "resolved": res, "selection_fully_developed": sel_fd, "selection_core": sel_core,
                "ranking_unchanged": sel_fd["order_by_r_th"] == sel_core["order_by_r_th"],
                "gates": gates_for(res, level1["recommended"]), "seconds": time.perf_counter() - t0}
        path.write_text(json.dumps(data, indent=2))
    plot_front(data["one_d"], requirements(level1), level1["recommended"], out_dir / "front_entrance.png")
    write_summary(data, out_dir / "summary.md")
    for k, v in data["gates"].items():
        print(f"{k}: {v['value']} (limit {v['limit']}) {'pass' if v['pass'] else 'FAIL'}")
    print(f"wrote {out_dir}")
    return 0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--plot-only", action="store_true", help="redraw the figure and summary from results.json")
    ap.add_argument("--level1", type=Path, default=LEVEL1, help="Level 1 results.json to read the designs from")
    ap.add_argument("--out", type=Path, default=OUT, help="output folder")
    a = ap.parse_args()
    sys.exit(run(a.plot_only, a.level1, a.out))


if __name__ == "__main__":
    main()
