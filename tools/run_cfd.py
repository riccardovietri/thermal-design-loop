"""run_cfd: 2D conjugate, fully developed laminar solve on the unit-cell cross-section (scikit-fem, P2).

Momentum:  -lap(u) = 1 in the fluid, u = 0 on walls;  w = (G/mu) u,  G = -dp/dx set by the flow rate.
Energy:    div(k grad T) = rho cp w dT/dx in the fluid, 0 in the solid; uniform flux on y = 0;
           dT/dx = Q / (rho cp V_dot L) from the overall energy balance; all other boundaries adiabatic/symmetry.
The temperature level is fixed by setting the local bulk fluid temperature to the axial mean bulk temperature.
"""
import json
from pathlib import Path

import numpy as np
from skfem import Basis, BilinearForm, ElementTriP2, FacetBasis, Functional, LinearForm, MeshTri, asm, condense, solve
from skfem.helpers import dot, grad

from .types import CfdResult, Design, Fluid, Requirements, to_dict


@BilinearForm
def _laplace(u, v, _):
    return dot(grad(u), grad(v))


@BilinearForm
def _mass(u, v, _):
    return u * v


@LinearForm
def _unit(v, _):
    return 1.0 * v


@Functional
def _dtdy(w):
    return w["temp"].grad[1]


def run_cfd(design: Design, req: Requirements, fluid: Fluid, mesh_path: Path, out_dir: Path | None = None) -> CfdResult:
    data = np.load(mesh_path)
    m = MeshTri(np.ascontiguousarray(data["p"]), np.ascontiguousarray(data["t"]))
    fluid_el = np.nonzero(data["fluid"])[0]
    solid_el = np.nonzero(~data["fluid"])[0]
    elem = ElementTriP2()
    basis = Basis(m, elem)
    basis_f = Basis(m, elem, elements=fluid_el)
    basis_s = Basis(m, elem, elements=solid_el)

    t2, tb, h = design.t_fin_m / 2, req.base_thickness_m, design.h_fin_m
    half_pitch = t2 + design.w_ch_m / 2
    tol = 1e-6 * half_pitch
    n_ch = design.n_channels(req.footprint_width_m)
    length = req.flow_length_m

    # --- momentum ---
    walls = m.facets_satisfying(
        lambda x: ((np.abs(x[0] - t2) < tol) & (x[1] > tb - tol))
        | ((np.abs(x[1] - tb) < tol) & (x[0] > t2 - tol))
        | ((np.abs(x[1] - (tb + h)) < tol) & (x[0] > t2 - tol))
    )
    wall_dofs = basis.get_dofs(facets=walls).all()
    free = np.setdiff1d(basis.get_dofs(elements=fluid_el).all(), wall_dofs)
    fixed = np.setdiff1d(np.arange(basis.N), free)
    k_f, m_f, f_f = asm(_laplace, basis_f), asm(_mass, basis_f), asm(_unit, basis_f)
    u = solve(*condense(k_f, f_f, D=fixed))
    int_u = f_f @ u  # m^4, half channel
    area_half = design.w_ch_m / 2 * h
    vdot = design.flow_lpm / 60e3
    grad_p = fluid.mu * vdot / (2 * n_ch * int_u)  # Pa/m
    dp = grad_p * length
    w_mean = vdot / (2 * n_ch * area_half)
    dh = 2 * design.w_ch_m * h / (design.w_ch_m + h)
    reynolds = fluid.rho * w_mean * dh / fluid.mu
    f_re = dh**2 / (2 * int_u / area_half)

    # --- energy ---
    q_cell_per_len = req.heat_load_w / (2 * n_ch * length)  # W/m into the half cell
    flux = q_cell_per_len / half_pitch
    dtdx = req.heat_load_w / (fluid.rho * fluid.cp * vdot * length)
    bottom = m.facets_satisfying(lambda x: np.abs(x[1]) < tol)
    fb_bottom = FacetBasis(m, elem, facets=bottom)
    b_bottom = asm(_unit, fb_bottom)
    a = req.solid_k * asm(_laplace, basis_s) + fluid.k * k_f
    rhs = flux * b_bottom - fluid.rho * fluid.cp * dtdx * (grad_p / fluid.mu) * (m_f @ u)
    pin = np.array([basis.nodal_dofs[0, 0]])
    temp = solve(*condense(a, rhs, D=pin))
    residual = a @ temp - rhs
    energy_residual = abs(residual[pin[0]]) / q_cell_per_len
    linear_residual = np.linalg.norm(np.delete(residual, pin)) / np.linalg.norm(rhs)
    # Independent check: conducted flux recovered from the temperature gradient at the chip interface vs applied flux.
    recovered = -req.solid_k * _dtdy.assemble(fb_bottom, temp=fb_bottom.interpolate(temp))
    base_flux_error = recovered / q_cell_per_len - 1
    temp -= (u @ (m_f @ temp)) / int_u  # local bulk temperature -> 0

    coolant_rise = req.heat_load_w / (fluid.rho * fluid.cp * vdot)
    base_avg = (b_bottom @ temp) / half_pitch
    base_max = temp[basis.get_dofs(facets=bottom).all()].max()
    r_th = (base_avg + coolant_rise / 2) / req.heat_load_w
    r_th_max = (base_max + coolant_rise) / req.heat_load_w
    # Effective Nu: wetted-perimeter-average h referenced to (mean floor temperature - bulk), on Dh.
    floor = m.facets_satisfying(lambda x: (np.abs(x[1] - tb) < tol) & (x[0] > t2 - tol))
    b_floor = asm(_unit, FacetBasis(m, elem, facets=floor))
    floor_avg = (b_floor @ temp) / (design.w_ch_m / 2)
    nu_eff = q_cell_per_len / ((design.w_ch_m / 2 + h) * floor_avg) * dh / fluid.k

    result = CfdResult(
        design_id=design.design_id, n_channels=n_ch, reynolds=float(reynolds), dp_pa=float(dp),
        pump_power_w=float(dp * vdot), f_re_fanning=float(f_re), nu_eff=float(nu_eff),
        r_th_k_per_w=float(r_th), r_th_max_k_per_w=float(r_th_max),
        base_avg_temp_c=float(req.inlet_temp_c + coolant_rise / 2 + base_avg),
        base_max_temp_c=float(req.inlet_temp_c + coolant_rise + base_max),
        coolant_rise_k=float(coolant_rise), energy_residual_rel=float(energy_residual), n_dofs=int(basis.N),
        linear_residual_rel=float(linear_residual), base_flux_error_rel=float(base_flux_error),
    )
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "inputs.json").write_text(json.dumps(
            {"design": to_dict(design), "requirements": to_dict(req), "fluid": to_dict(fluid)}, indent=2))
        (out_dir / "results.json").write_text(json.dumps(to_dict(result), indent=2))
        np.savez(out_dir / "fields.npz", u=u * grad_p / fluid.mu, temp=temp, dofs_x=basis.doflocs[0], dofs_y=basis.doflocs[1])
    return result
