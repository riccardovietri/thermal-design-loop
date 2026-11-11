"""Level 2: 2D planar steady laminar Navier-Stokes + energy, nondimensional (scikit-fem, Taylor-Hood, P2 temperature).

Lengths / b (reference width), velocities / V (mean inlet velocity), pressure / (rho V^2), Re = rho V b / mu.
Temperature theta = (T - T_inlet) k_f / (q'' b) with uniform flux q'' on the heated wall, Pe = Re Pr.
"""
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from scipy.sparse import bmat
from skfem import Basis, BilinearForm, ElementTriP1, ElementTriP2, ElementVector, FacetBasis, LinearForm, MeshTri, asm, condense, solve
from skfem.helpers import ddot, div, dot, grad, mul


@BilinearForm
def _viscous(u, v, _):
    return ddot(grad(u), grad(v))


@BilinearForm
def _divergence(u, q, _):
    return -div(u) * q


@BilinearForm
def _advection_jac(u, v, w):
    z = w["wind"]
    return dot(mul(grad(u), z), v) + dot(mul(grad(z), u), v)


@LinearForm
def _advection_res(v, w):
    z = w["wind"]
    return dot(mul(grad(z), z), v)


@BilinearForm
def _laplace(u, v, _):
    return dot(grad(u), grad(v))


@BilinearForm
def _scalar_advection(u, v, w):
    return dot(w["wind"], grad(u)) * v


@LinearForm
def _unit(v, _):
    return 1.0 * v


@LinearForm
def _vorticity(v, w):
    z = w["wind"]
    return (z.grad[1][0] - z.grad[0][1]) * v


@dataclass
class FlowSolution:
    mesh: MeshTri
    basis_u: Basis
    basis_p: Basis
    u: np.ndarray
    p: np.ndarray
    reynolds: float
    newton_iterations: int
    residual_norm: float

    def speed_at_nodes(self) -> np.ndarray:
        return np.hypot(self.u[self.basis_u.nodal_dofs[0]], self.u[self.basis_u.nodal_dofs[1]])


def component_dofs(basis_u: Basis, facets: np.ndarray, component: int) -> np.ndarray:
    dofs = basis_u.get_dofs(facets=facets).all()
    comp = np.concatenate([basis_u.nodal_dofs[component], basis_u.facet_dofs[component]])
    return dofs[np.isin(dofs, comp)]


def solve_flow(mesh: MeshTri, reynolds_steps: list[float], fixed: dict[int, list[tuple[np.ndarray, Callable | None]]],
               tol: float = 1e-9, max_iter: int = 30) -> FlowSolution:
    """Newton solve with Reynolds continuation.

    fixed maps velocity component (0 or 1) -> list of (facets, fn) where fn(x) gives the boundary value (None = zero).
    Boundaries not fixed are traction-free (outlets) or symmetric in the unfixed component.
    """
    bu = Basis(mesh, ElementVector(ElementTriP2()), intorder=4)
    bp = bu.with_element(ElementTriP1())
    nu = bu.N
    up = np.zeros(nu + bp.N)
    d_list = []
    for comp, entries in fixed.items():
        for facets, fn in entries:
            dofs = component_dofs(bu, facets, comp)
            d_list.append(dofs)
            if fn is not None:
                up[dofs] = fn(bu.doflocs[:, dofs])
    d_u = np.unique(np.concatenate(d_list))
    a_visc, b = asm(_viscous, bu), asm(_divergence, bu, bp)
    its, norm = 0, np.inf
    for re in reynolds_steps:
        for its in range(1, max_iter + 1):
            wind = bu.interpolate(up[:nu])
            jac = bmat([[a_visc / re + asm(_advection_jac, bu, wind=wind), b.T], [b, None]], format="csr")
            res = np.concatenate([a_visc @ up[:nu] / re + asm(_advection_res, bu, wind=wind) + b.T @ up[nu:], b @ up[:nu]])
            delta = solve(*condense(jac, -res, D=d_u))
            up += delta
            norm = np.linalg.norm(delta) / max(np.linalg.norm(up), 1e-30)
            if norm < tol:
                break
        else:
            raise RuntimeError(f"Newton did not converge at Re = {re}: last relative step {norm:.2e}")
    return FlowSolution(mesh, bu, bp, up[:nu], up[nu:], reynolds_steps[-1], its, float(norm))


def stream_function(flow: FlowSolution, wall_values: list[tuple[np.ndarray, float]]) -> np.ndarray:
    """psi with -lap(psi) = vorticity; wall_values gives (facets, psi) on boundaries where psi is known."""
    b1 = flow.basis_p
    psi0, d = np.zeros(b1.N), []
    for facets, value in wall_values:
        dofs = b1.get_dofs(facets=facets).all()
        psi0[dofs] = value
        d.append(dofs)
    rhs = asm(_vorticity, b1, wind=flow.basis_u.interpolate(flow.u))
    return solve(*condense(asm(_laplace, b1), rhs, x=psi0, D=np.unique(np.concatenate(d))))


@dataclass
class EnergySolution:
    basis: Basis
    theta: np.ndarray
    peclet: float
    energy_balance_error: float


def solve_energy(flow: FlowSolution, prandtl: float, inlet: np.ndarray, heated: np.ndarray, outlet: np.ndarray) -> EnergySolution:
    """Advection-diffusion with theta = 0 at the inlet and unit nondimensional flux on the heated facets."""
    bt = flow.basis_u.with_element(ElementTriP2())
    pe = flow.reynolds * prandtl
    wind = flow.basis_u.interpolate(flow.u)
    a = asm(_scalar_advection, bt, wind=wind) + asm(_laplace, bt) / pe
    fb_heat = FacetBasis(flow.mesh, ElementTriP2(), facets=heated)
    rhs = asm(_unit, fb_heat) / pe
    theta = solve(*condense(a, rhs, D=bt.get_dofs(facets=inlet).all()))
    # energy balance: enthalpy leaving through the outlet vs heat added (both divided by rho cp V b (q'' b / k))
    fb_out = FacetBasis(flow.mesh, ElementVector(ElementTriP2()), facets=outlet, intorder=4)
    fbt_out = fb_out.with_element(ElementTriP2())

    @LinearForm
    def outflow(v, w):
        return dot(w["wind"], w.n) * v

    convected = asm(outflow, fbt_out, wind=fb_out.interpolate(flow.u)) @ theta
    added = rhs.sum()
    return EnergySolution(bt, theta, pe, float(convected / added - 1))


def boundary_profile(sol: EnergySolution, facets: np.ndarray, axis: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Coordinates along `axis` and theta at the degrees of freedom on a boundary, sorted."""
    dofs = sol.basis.get_dofs(facets=facets).all()
    x = sol.basis.doflocs[axis, dofs]
    order = np.argsort(x)
    return x[order], sol.theta[dofs][order]
