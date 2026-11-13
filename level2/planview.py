"""Level 2 plan view: fin planforms in the tall-fin limit, put into a plate-level resistance network.

2D solve (`level2.planar`, Taylor-Hood, steady laminar Navier-Stokes + energy) in the plane of the fin layer. Fins are
extruded shapes of infinite height: floor and cover drag and heat transfer are not in the 2D solve.
Nondimensional: lengths / b with b = 1 mm, velocities / U_s (superficial velocity), pressure / (rho U_s^2),
Re_b = rho U_s b / mu, theta = (T - T_ref) k_f / (q'' b) with unit flux on every fin wall, Pe = Re_b Pr.

Production values are developed-periodic, from a streamwise-periodic cell (Patankar, Liu & Sparrow 1977):
  dp/dx = beta rho U_s^2 / b                                          (beta: mean pressure gradient of the cell)
  h     = q'' / (mean fin-wall T over the period - mean bulk T over the period) = k_f / (b (theta_w - theta_b)),
          bulk = velocity-weighted mean, i.e. integral(u_x theta dA) / integral(u_x dA) over the period.
A developing strip (uniform inflow, lead-in, N periods) checks how fast these values are reached (`run_developing`).

Network to plate level (same definitions as Level 1, `tools/screen_1d.py`), footprint A_fp = W x L = 60 x 60 mm:
  U_s       = Q / (W H)                         flow per lateral unit: total flow over the 60 mm width and fin height H
  dP        = (dp/dx) L
  R_th      = R_base + 1 / (h (A_floor + eta A_fin)) + 1 / (2 m_dot cp)
  R_base    = t_base / (k_s A_fp)
  A_floor   = porosity x A_fp                    (plan fluid area; the same h is applied to it)
  A_fin     = (wetted fin perimeter per unit plan area) x A_fp x H
  eta       = tanh(mH) / (mH),  m = sqrt(h P_c / (k_s A_c))   (adiabatic tip; P_c, A_c: plan-section perimeter and area
                                                               of one fin; a continuous wall of thickness t gives
                                                               m = sqrt(2 h / (k_s t)))
Hydraulic diameter Dh = 4 A_fluid / P_wetted in plan (2w for straight fins); Re_Dh = rho (U_s / porosity) Dh / mu,
i.e. on the mean interstitial velocity. For pins Re_d,max = rho U_max d / mu with U_max = U_s S_T / (S_T - d) is also
reported. Entrance effects are excluded: the developed-periodic values are applied over the whole length.
"""
import json
import math
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
from matplotlib.tri import LinearTriInterpolator, Triangulation

from tools import plotstyle as ps
from tools.types import Fluid, Requirements

from .planar import (NewtonDivergence, line_mean, solve_energy, solve_energy_periodic, solve_flow, solve_flow_periodic,
                     window_means)
from .planforms import Planform, mesh_planform

B_M = 1e-3  # reference length, 1 mm


def superficial_velocity(req: Requirements, flow_lpm: float, h_fin_mm: float) -> float:
    return flow_lpm / 60e3 / (req.footprint_width_m * h_fin_mm * 1e-3)


def reynolds_b(req: Requirements, fluid: Fluid, flow_lpm: float, h_fin_mm: float) -> float:
    return fluid.rho * superficial_velocity(req, flow_lpm, h_fin_mm) * B_M / fluid.mu


def continuation(re: float) -> list[float]:
    """Reynolds steps for Newton; solve_flow_periodic also bisects a failed step."""
    return [r for r in (10.0, 50.0, 100.0, 200.0) if r < 0.8 * re] + [re]


def network(pf: Planform, dpdx_pa_m: float, h: float, req: Requirements, fluid: Fluid, flow_lpm: float) -> dict:
    """Plate-level dP and R_th from developed-periodic dp/dx and h (formulas in the module docstring)."""
    h_fin = pf.params["h_fin"] * 1e-3
    a_fp = req.footprint_width_m * req.flow_length_m
    a_floor = pf.porosity * a_fp
    a_fin = pf.perimeter_density_per_mm * 1e3 * a_fp * h_fin
    m = math.sqrt(h * pf.fin_section_perimeter_mm * 1e-3 / (req.solid_k * pf.fin_section_area_mm2 * 1e-6))
    eta = math.tanh(m * h_fin) / (m * h_fin)
    r_conv = 1 / (h * (a_floor + eta * a_fin))
    r_base = req.base_thickness_m / (req.solid_k * a_fp)
    mdot = fluid.rho * flow_lpm / 60e3
    r_cal = 1 / (2 * mdot * fluid.cp)
    return {"dp_pa": dpdx_pa_m * req.flow_length_m, "r_th_k_per_w": r_base + r_conv + r_cal, "r_base": r_base,
            "r_conv": r_conv, "r_caloric": r_cal, "fin_efficiency": eta, "m_per_m": m, "a_floor_m2": a_floor, "a_fin_m2": a_fin,
            "wall_flux_w_m2": req.heat_load_w / (a_floor + eta * a_fin)}


def design_id(pf: Planform) -> str:
    return pf.topology + "_" + "_".join(f"{k}{v * 1000:04.0f}" for k, v in pf.params.items() if k != "r_in")


@dataclass
class PlanviewResult:
    design_id: str
    topology: str
    params_mm: dict
    flow_lpm: float
    reynolds_b: float  # rho U_s (1 mm) / mu, the solver Reynolds number
    reynolds_dh: float  # rho (U_s / porosity) Dh / mu
    reynolds_pin_max: float | None
    dh_mm: float
    dpdx_pa_per_m: float
    h_w_m2k: float
    nu_dh: float  # h Dh / k_f
    dp_pa: float
    r_th_k_per_w: float
    r_base: float
    r_conv: float
    r_caloric: float
    fin_efficiency: float
    a_floor_m2: float
    a_fin_m2: float
    energy_balance_error: float
    undershoot_rel: float  # thermal undershoot below the inflow minimum, relative to theta_w - theta_b
    newton_iterations: int
    reynolds_path: list
    n_unknowns: int
    n_triangles: int
    cells: float
    n_periods: int
    windows: list = field(default_factory=list)  # per-period dp/dx and h
    seconds: float = 0.0
    geometry: dict = field(default_factory=dict)


def _flow_bcs(pm):
    return {0: [pm.walls], 1: [pm.walls] + ([pm.symmetry] if len(pm.symmetry) else [])}


def solve_periodic(pf: Planform, req: Requirements, fluid: Fluid, flow_lpm: float = 1.5, cells: float = 8.0, n_periods: int = 1):
    pm = mesh_planform(pf, n_periods=n_periods, periodic=True, cells=cells)
    re = reynolds_b(req, fluid, flow_lpm, pf.params["h_fin"])
    pr = fluid.cp * fluid.mu / fluid.k
    flow = solve_flow_periodic(pm.mesh, continuation(re), _flow_bcs(pm), pm.inlet, pm.outlet, pm.x_out - pm.x_in, pm.flow_rate)
    energy = solve_energy_periodic(flow, pr, pm.heated, pm.inlet, pm.outlet)
    windows = []
    for k in range(n_periods):
        lo, hi = pm.cuts[k], pm.cuts[k + 1]
        wm = window_means(flow, energy, pm.heated, lo, hi)
        dp_line = line_mean(pm.mesh, flow.basis_p, flow.p, lo) - line_mean(pm.mesh, flow.basis_p, flow.p, hi)
        windows.append({"beta": flow.pressure_gradient + dp_line / (hi - lo), **wm})
    return pm, flow, energy, windows


def _dimensional(pf, req, fluid, flow_lpm, beta, nu_b):
    us = superficial_velocity(req, flow_lpm, pf.params["h_fin"])
    return beta * fluid.rho * us**2 / B_M, nu_b * fluid.k / B_M


def run_design(pf: Planform, req: Requirements, fluid: Fluid, flow_lpm: float = 1.5, cells: float = 8.0, n_periods: int = 1,
               out_dir: Path | None = None, plot: bool = False) -> PlanviewResult:
    t0 = time.perf_counter()
    pm, flow, energy, windows = solve_periodic(pf, req, fluid, flow_lpm, cells, n_periods)
    beta = float(np.mean([w["beta"] for w in windows]))
    w_all = window_means(flow, energy, pm.heated, pm.cuts[0], pm.cuts[-1])
    dpdx, h = _dimensional(pf, req, fluid, flow_lpm, beta, w_all["nusselt_b"])
    net = network(pf, dpdx, h, req, fluid, flow_lpm)
    us = superficial_velocity(req, flow_lpm, pf.params["h_fin"])
    re_dh = fluid.rho * us / pf.porosity * pf.dh_mm * 1e-3 / fluid.mu
    re_pin = None
    if pf.topology.startswith("pin"):
        size = pf.params.get("d", pf.params.get("s"))
        extent = size if pf.topology == "pin_circle" else size * math.sqrt(2)
        s_t = pf.params["s_t"]
        re_pin = fluid.rho * us * s_t / (s_t - extent) * size * 1e-3 / fluid.mu
    for w in windows:
        w["dpdx_pa_per_m"], w["h_w_m2k"] = _dimensional(pf, req, fluid, flow_lpm, w["beta"], w["nusselt_b"])
    res = PlanviewResult(
        design_id=design_id(pf), topology=pf.topology, params_mm=dict(pf.params), flow_lpm=flow_lpm,
        reynolds_b=flow.reynolds, reynolds_dh=re_dh, reynolds_pin_max=re_pin, dh_mm=pf.dh_mm, dpdx_pa_per_m=dpdx, h_w_m2k=h,
        nu_dh=h * pf.dh_mm * 1e-3 / fluid.k, dp_pa=net["dp_pa"], r_th_k_per_w=net["r_th_k_per_w"], r_base=net["r_base"],
        r_conv=net["r_conv"], r_caloric=net["r_caloric"], fin_efficiency=net["fin_efficiency"], a_floor_m2=net["a_floor_m2"],
        a_fin_m2=net["a_fin_m2"], energy_balance_error=energy.energy_balance_error,
        undershoot_rel=energy.undershoot * w_all["nusselt_b"], newton_iterations=flow.newton_iterations,
        reynolds_path=list(flow.reynolds_path), n_unknowns=int(flow.basis_u.N + flow.basis_p.N + 1 + energy.basis.N),
        n_triangles=int(pm.mesh.t.shape[1]), cells=cells, n_periods=n_periods, windows=windows,
        seconds=time.perf_counter() - t0, geometry=pf.summary())
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "results.json").write_text(json.dumps(asdict(res), indent=2))
        if plot:
            plot_fields(pf, pm, flow, energy, res, net, req, fluid, out_dir / "fields.png")
    return res


def run_developing(pf: Planform, req: Requirements, fluid: Fluid, n_periods: int, flow_lpm: float = 1.5, cells: float = 8.0,
                   lead_in_mm: float = 1.0) -> dict:
    """Uniform inflow (theta = 0) through a lead-in and n_periods analysed periods; per-period dp/dx and h.

    One further period of the same geometry plus a short empty buffer sits between the last analysed period and the
    traction-free outlet, so the outlet does not disturb the analysed periods.
    """
    t0 = time.perf_counter()
    pm = mesh_planform(pf, n_periods=n_periods + 1, periodic=False, cells=cells, lead_in_mm=lead_in_mm)
    re = reynolds_b(req, fluid, flow_lpm, pf.params["h_fin"])
    pr = fluid.cp * fluid.mu / fluid.k
    u_in = pm.flow_rate / pm.inlet_width
    bcs = _flow_bcs(pm)
    fixed = {0: [(pm.inlet, lambda x: u_in + 0 * x[0])] + [(f, None) for f in bcs[0]],
             1: [(pm.inlet, None)] + [(f, None) for f in bcs[1]]}
    flow = solve_flow(pm.mesh, continuation(re), fixed)
    energy = solve_energy(flow, pr, pm.inlet, pm.heated, pm.outlet)
    periods = []
    for k in range(n_periods):
        lo, hi = pm.cuts[k], pm.cuts[k + 1]
        wm = window_means(flow, energy, pm.heated, lo, hi)
        beta = (line_mean(pm.mesh, flow.basis_p, flow.p, lo) - line_mean(pm.mesh, flow.basis_p, flow.p, hi)) / (hi - lo)
        dpdx, h = _dimensional(pf, req, fluid, flow_lpm, beta, wm["nusselt_b"])
        periods.append({"period": k + 1, "x_start_mm": float(lo), "beta": float(beta), "nusselt_b": wm["nusselt_b"],
                        "dpdx_pa_per_m": dpdx, "h_w_m2k": h})
    last, prev = periods[-1], periods[-2]
    theta_min = float(energy.theta.min())
    return {"design_id": design_id(pf), "n_periods": n_periods, "cells": cells, "lead_in_mm": lead_in_mm,
            "length_mm": float(pm.cuts[-1]), "reynolds_b": flow.reynolds, "energy_balance_error": energy.energy_balance_error,
            "theta_min": theta_min, "theta_min_rel": theta_min * last["nusselt_b"], "periods": periods,
            "last_two_dpdx_change_pct": 100 * (last["dpdx_pa_per_m"] / prev["dpdx_pa_per_m"] - 1),
            "last_two_h_change_pct": 100 * (last["h_w_m2k"] / prev["h_w_m2k"] - 1),
            "n_unknowns": int(flow.basis_u.N + flow.basis_p.N + energy.basis.N), "seconds": time.perf_counter() - t0}


# --------------------------------------------------------------------------------------------------------------------
# Field figure


def _tiled(pf, pm, flow, energy, n_tile):
    """Vertex coordinates, triangles and fields tiled over n_tile periods and mirrored to show a full lateral repeat."""
    m = flow.mesh
    nv = m.p.shape[1]
    ux = flow.u[flow.basis_u.nodal_dofs[0]]
    uy = flow.u[flow.basis_u.nodal_dofs[1]]
    theta = energy.theta[energy.basis.nodal_dofs[0]]
    x_rel = m.p[0] - flow.x_master
    t_rel = theta - energy.gamma * x_rel  # periodic part, follows the cell when tiled
    length = flow.length
    copies = []
    mirrors = [(1, 0.0)]
    if pf.kind == "strip":
        mirrors.append((-1, 0.0))
    elif pf.kind == "half_channel":
        mirrors.append((-1, 2 * pf.geo["half_width"]))
    else:  # wavy channel: add the neighbouring channel so the fin between them is visible
        mirrors.append((1, pf.lateral_unit_mm))
    for k in range(n_tile):
        for sign, off in mirrors:
            copies.append((m.p[0] + k * length, sign * m.p[1] + off, ux, sign * uy, t_rel, sign))
    px = np.concatenate([c[0] for c in copies])
    py = np.concatenate([c[1] for c in copies])
    tris = np.concatenate([(m.t.T if c[5] > 0 else m.t.T[:, ::-1]) + i * nv for i, c in enumerate(copies)])
    fields = [np.concatenate([c[j] for c in copies]) for j in (2, 3, 4)]
    # merge the duplicated vertices on tile seams (values agree there up to round-off)
    key = np.round(np.column_stack([px, py]) / 1e-7).astype(np.int64)
    _, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = inv.ravel()
    return (px[first], py[first], inv[tris], *[f[first] for f in fields])


def describe(topology: str, params: dict) -> str:
    """Plain-words name of a planform design, e.g. 'Circular pins, 0.4 mm diameter, 0.8 × 0.9 mm pitch'."""
    p = {k: f"{v:g}" for k, v in params.items()}
    wave = f"{p.get('amp')} mm amplitude, {p.get('lam')} mm wavelength"
    return {
        "straight": f"Straight fins, {p.get('w')} mm channels, {p.get('t')} mm fins",
        "sinusoidal": f"Sinusoidal fins, {p.get('w')} mm channels, {p.get('t')} mm fins, {wave}",
        "zigzag": f"Zigzag fins, {p.get('w')} mm channels, {p.get('t')} mm fins, {wave}",
        "offset_strip": f"Offset strip fins, {p.get('l')} mm long, {p.get('t')} mm thick, {p.get('p_t')} mm lateral pitch",
        "pin_circle": f"Circular pins, {p.get('d')} mm diameter, {p.get('s_t')} × {p.get('s_l')} mm pitch",
        "pin_diamond": f"Diamond pins, {p.get('s')} mm side, {p.get('s_t')} × {p.get('s_l')} mm pitch",
    }[topology]


def plot_fields(pf, pm, flow, energy, res, net, req, fluid, out: Path) -> Path:
    us = superficial_velocity(req, res.flow_lpm, pf.params["h_fin"])
    target_len = max(4 * pf.period_mm if pf.topology == "straight" else 2 * pf.period_mm, 3.0)
    n_tile = max(1, int(math.ceil(target_len / pf.period_mm - 1e-9)))
    n_tile = min(n_tile, 6)
    px, py, tris, ux, uy, t_rel = _tiled(pf, pm, flow, energy, n_tile)
    tri = Triangulation(px, py, tris)
    speed = np.hypot(ux, uy) * us
    w_all = window_means(flow, energy, pm.heated, pm.cuts[0], pm.cuts[-1])
    bulk_offset = w_all["theta_bulk"] - energy.gamma * 0.5 * flow.length  # mean of theta - gamma x over bulk
    scale = net["wall_flux_w_m2"] * B_M / fluid.k  # K per unit theta
    temp = (t_rel - bulk_offset) * scale
    xr, yr = (px.min(), px.max()), (py.min(), py.max())
    aspect = (yr[1] - yr[0]) / (xr[1] - xr[0])
    width = ps.WIDTH_IN
    panel_h = max(1.4, min(3.0, (width - 1.6) * aspect))
    fig, axes = ps.figure(2 * panel_h + 1.5, 2, 1, sharex=True)
    cs = axes[0].tricontourf(tri, speed, levels=ps.field_levels(0, speed.max()), cmap="viridis")
    ps.colorbar(fig, cs, axes[0], "Speed (m/s)", 0, float(speed.max()), nbins=4, min_ticks=3, shrink=0.95)
    nx = 400
    ny = max(40, int(nx * aspect))
    gx, gy = np.meshgrid(np.linspace(*xr, nx), np.linspace(*yr, ny))
    fu, fv = LinearTriInterpolator(tri, ux)(gx, gy), LinearTriInterpolator(tri, uy)(gx, gy)
    kw = {}
    if pf.topology == "straight":  # parallel flow: a few lines seeded across the channel, away from the walls
        kw["start_points"] = np.column_stack([np.full(5, xr[0] + 1e-6), np.linspace(yr[0], yr[1], 7)[1:-1]])
    axes[0].streamplot(gx, gy, np.ma.filled(fu, np.nan), np.ma.filled(fv, np.nan), color="white", linewidth=0.7,
                       density=(0.8, max(0.5, 0.8 * aspect * 2.5)), arrowsize=0.8,
                       broken_streamlines=pf.kind == "strip", **kw)  # broken lines keep recirculation cells from filling in
    cs = axes[1].tricontourf(tri, temp, levels=ps.field_levels(temp.min(), temp.max()), cmap="inferno")
    ps.colorbar(fig, cs, axes[1], "T − bulk T (K)", float(temp.min()), float(temp.max()), nbins=4, min_ticks=3, shrink=0.95)
    for ax in axes:
        ax.set_facecolor(ps.FIN_GREY)  # fins (outside the fluid) in grey
        ax.set_aspect("equal")
        ax.set_xlim(*xr)
        ax.set_ylim(*yr)
        ax.set_ylabel("y (mm)")
        ps.no_grid(ax)
    axes[1].set_xlabel("x, flow direction (mm)")
    axes[0].set_title("Speed and streamlines", fontsize=11.5, loc="left")
    axes[1].set_title("Temperature above the local bulk temperature", fontsize=11.5, loc="left")
    fig.suptitle(describe(pf.topology, pf.params), fontsize=12.5)
    ps.save(fig, out)
    return out
