"""Level 2: 2D planar steady laminar Navier-Stokes + energy, nondimensional (scikit-fem, Taylor-Hood, P2 temperature).

Lengths / b (reference width), velocities / V (mean inlet velocity), pressure / (rho V^2), Re = rho V b / mu.
Temperature theta = (T - T_inlet) k_f / (q'' b) with uniform flux q'' on the heated wall, Pe = Re Pr.

Streamwise-periodic variants (`solve_flow_periodic`, `solve_energy_periodic`) solve one or more geometric periods of a
fully developed (periodic) flow, following Patankar, Liu & Sparrow (1977): p = -beta x + periodic part, with beta found
from an imposed flow rate, and theta = gamma x + periodic part with gamma fixed by the heat added per period.
"""
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from scipy.sparse import bmat, block_diag, csr_matrix, identity
from scipy.spatial import cKDTree
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


@LinearForm
def _x_component(v, _):
    return v[0]


@LinearForm
def _wind_x(v, w):
    return w["wind"][0] * v


@LinearForm
def _dvdx(v, _):
    return v.grad[0]


class NewtonDivergence(RuntimeError):
    """Newton failed at a Reynolds number; keeps the Reynolds number and the last relative step for the record."""

    def __init__(self, reynolds: float, last_step: float):
        super().__init__(f"Newton did not converge at Re = {reynolds}: last relative step {last_step:.2e}")
        self.reynolds, self.last_step = float(reynolds), float(last_step)


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
            raise NewtonDivergence(re, norm)
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


# ---------------------------------------------------------------------------------------------------------------------
# Streamwise-periodic (fully developed) flow and heat transfer


def periodic_reduction(basis: Basis, master: np.ndarray, slave: np.ndarray, shift: float, tol: float = 1e-7):
    """Identify each degree of freedom on the slave facets (x = x_master + shift) with its master twin.

    Returns (index map full -> reduced, prolongation matrix P with full = P @ reduced).
    Components of a vector basis are kept apart, so u_x pairs with u_x and u_y with u_y.
    """
    comp = np.zeros(basis.N)
    if basis.nodal_dofs.shape[0] > 1:
        for c in range(basis.nodal_dofs.shape[0]):
            comp[basis.nodal_dofs[c]] = c
            if basis.facet_dofs.size:
                comp[basis.facet_dofs[c]] = c
    m_dofs = basis.get_dofs(facets=master).all()
    s_dofs = basis.get_dofs(facets=slave).all()
    key = lambda d, dx: np.column_stack([basis.doflocs[0, d] + dx, basis.doflocs[1, d], 10.0 * comp[d]])
    dist, j = cKDTree(key(m_dofs, shift)).query(key(s_dofs, 0.0))
    if dist.max() > tol or len(m_dofs) != len(s_dofs):
        raise ValueError(f"periodic boundaries do not match (worst distance {dist.max():.2e}, "
                         f"{len(m_dofs)} master vs {len(s_dofs)} slave dofs)")
    index = np.arange(basis.N)
    index[s_dofs] = m_dofs[j]
    _, red = np.unique(index, return_inverse=True)
    p = csr_matrix((np.ones(basis.N), (np.arange(basis.N), red)), shape=(basis.N, red.max() + 1))
    return red, p


@dataclass
class PeriodicFlowSolution(FlowSolution):
    """p is the periodic part of the pressure; the full pressure is p - pressure_gradient * x."""
    pressure_gradient: float = 0.0  # beta = -d(mean p)/dx, nondimensional (rho V^2 / b)
    length: float = 0.0
    flow_rate: float = 0.0
    x_master: float = 0.0
    reynolds_path: tuple = ()


def solve_flow_periodic(mesh: MeshTri, reynolds_steps: list[float], walls: dict[int, list[np.ndarray]],
                        master: np.ndarray, slave: np.ndarray, length: float, flow_rate: float,
                        tol: float = 1e-9, max_iter: int = 25, max_halvings: int = 4) -> PeriodicFlowSolution:
    """Fully developed periodic flow through one or more periods, at an imposed flow rate (per unit depth).

    walls maps velocity component -> facet arrays where that component is zero (no-slip: both; symmetry: component 1).
    master/slave are the upstream/downstream periodic boundaries, slave = master shifted by `length` in x.
    The mean pressure gradient beta is an unknown with the constraint (integral of u_x over master) = flow_rate.
    Reynolds continuation as in `solve_flow`; a failed step is retried from the last converged state at the midpoint
    Reynolds number, up to `max_halvings` times, before NewtonDivergence is raised.
    """
    bu = Basis(mesh, ElementVector(ElementTriP2()), intorder=4)
    bp = bu.with_element(ElementTriP1())
    nu, npr = bu.N, bp.N
    red_u, pu = periodic_reduction(bu, master, slave, length)
    red_p, pp = periodic_reduction(bp, master, slave, length)
    pfull = block_diag((pu, pp, identity(1)), format="csr")
    nur, npr_r = pu.shape[1], pp.shape[1]
    d_full = np.unique(np.concatenate([component_dofs(bu, f, c) for c, fl in walls.items() for f in fl]))
    d_red = np.unique(np.concatenate([red_u[d_full], [nur + red_p[0]]]))
    a_visc, b = asm(_viscous, bu), asm(_divergence, bu, bp)
    f = asm(_x_component, bu)
    g = asm(_x_component, FacetBasis(mesh, ElementVector(ElementTriP2()), facets=master, intorder=4))
    f_col, g_row = csr_matrix(-f[:, None]), csr_matrix(g[None, :])
    x = np.zeros(nu + npr + 1)

    def newton(x, re):
        norm = np.inf
        for its in range(1, max_iter + 1):
            u, p, beta = x[:nu], x[nu:nu + npr], x[-1]
            wind = bu.interpolate(u)
            jac = bmat([[a_visc / re + asm(_advection_jac, bu, wind=wind), b.T, f_col], [b, None, None], [g_row, None, None]],
                       format="csr")
            res = np.concatenate([a_visc @ u / re + asm(_advection_res, bu, wind=wind) + b.T @ p - beta * f, b @ u,
                                  [g @ u - flow_rate]])
            jr, rr = (pfull.T @ jac @ pfull).tocsr(), pfull.T @ res
            delta = pfull @ solve(*condense(jr, -rr, D=d_red))
            x = x + delta
            norm = np.linalg.norm(delta[:nu]) / max(np.linalg.norm(x[:nu]), 1e-30)
            if norm < tol:
                return x, its, norm
        raise NewtonDivergence(re, norm)

    path, its, norm, re_done = [], 0, np.inf, 0.0
    for target in reynolds_steps:
        re_try, halvings = target, 0
        while True:
            try:
                x, its, norm = newton(x, re_try)
                path.append(float(re_try))
                re_done = re_try
                if re_try >= target:
                    break
                re_try = target
            except NewtonDivergence:
                halvings += 1
                if halvings > max_halvings:
                    raise
                re_try = 0.5 * (re_done + re_try)
    return PeriodicFlowSolution(mesh, bu, bp, x[:nu], x[nu:nu + npr], reynolds_steps[-1], its, float(norm),
                                pressure_gradient=float(x[-1]), length=length, flow_rate=flow_rate,
                                x_master=float(mesh.p[0, mesh.facets[:, master]].mean()), reynolds_path=tuple(path))


@dataclass
class PeriodicEnergySolution(EnergySolution):
    """theta = gamma (x - x_master) + periodic part; theta_periodic is the periodic part (zero at one point)."""
    gamma: float = 0.0
    theta_periodic: np.ndarray | None = None
    heated_length: float = 0.0
    undershoot: float = 0.0  # min theta in the domain minus min theta on the inflow boundary (<= 0; maximum principle)


def solve_energy_periodic(flow: PeriodicFlowSolution, prandtl: float, heated: np.ndarray, master: np.ndarray,
                          slave: np.ndarray) -> PeriodicEnergySolution:
    """Developed periodic temperature with unit nondimensional flux on the heated facets, other boundaries adiabatic.

    gamma = (heated length / Pe) / (flow rate x period length): the bulk rise per unit length from the energy balance.
    The global energy equation (test function 1) is the dropped row; its residual is reported as energy_balance_error.
    """
    bt = flow.basis_u.with_element(ElementTriP2())
    pe = flow.reynolds * prandtl
    wind = flow.basis_u.interpolate(flow.u)
    k_adv = asm(_scalar_advection, bt, wind=wind)
    a = k_adv + asm(_laplace, bt) / pe
    heat = asm(_unit, FacetBasis(flow.mesh, ElementTriP2(), facets=heated))
    s_heat = heat.sum()
    gamma = s_heat / (pe * flow.flow_rate * flow.length)
    rhs = heat / pe - gamma * asm(_wind_x, bt, wind=wind) - gamma / pe * asm(_dvdx, bt)
    red, pt = periodic_reduction(bt, master, slave, flow.length)
    ar, br = (pt.T @ a @ pt).tocsr(), pt.T @ rhs
    phi = pt @ solve(*condense(ar, br, D=np.array([0])))
    theta = gamma * (bt.doflocs[0] - flow.x_master) + phi
    added = s_heat / pe
    err = (k_adv @ theta).sum() / added - 1
    inflow = bt.get_dofs(facets=master).all()
    return PeriodicEnergySolution(bt, theta, pe, float(err), gamma=float(gamma), theta_periodic=phi, heated_length=float(s_heat),
                                  undershoot=float(theta.min() - theta[inflow].min()))


def line_mean(mesh: MeshTri, basis_p1: Basis, values: np.ndarray, x0: float, tol: float = 1e-7) -> float:
    """Length-weighted mean of a P1 field over the mesh edges lying on the line x = x0 (fluid part only)."""
    f = mesh.facets
    on = np.all(np.abs(mesh.p[0, f] - x0) < tol, axis=0)
    a, b = f[0, on], f[1, on]
    if not on.any():
        raise ValueError(f"no mesh edges on x = {x0}")
    lengths = np.abs(mesh.p[1, a] - mesh.p[1, b])
    dof = basis_p1.nodal_dofs[0]
    return float(np.sum(lengths * 0.5 * (values[dof[a]] + values[dof[b]])) / lengths.sum())


def window_means(flow: FlowSolution, energy: EnergySolution, heated: np.ndarray, x_lo: float, x_hi: float) -> dict:
    """Mean wall temperature and velocity-weighted bulk temperature over x_lo <= x <= x_hi (mesh must conform to both lines).

    Bulk: integral of u_x theta over the window area / integral of u_x (equal to the x-average of the cross-section
    bulk temperature when the flow rate is the same through every cross-section).
    """
    mesh = flow.mesh
    xc = mesh.p[0, mesh.t].mean(axis=0)
    elems = np.nonzero((xc > x_lo) & (xc < x_hi))[0]
    fmid = mesh.p[0, mesh.facets[:, heated]].mean(axis=0)
    walls = heated[(fmid > x_lo) & (fmid < x_hi)]
    bt = Basis(mesh, ElementTriP2(), elements=elems, intorder=4)
    bu = Basis(mesh, ElementVector(ElementTriP2()), elements=elems, intorder=4)
    wx = asm(_wind_x, bt, wind=bu.interpolate(flow.u))
    fw = asm(_unit, FacetBasis(mesh, ElementTriP2(), facets=walls))
    theta_b, theta_w = (wx @ energy.theta) / wx.sum(), (fw @ energy.theta) / fw.sum()
    return {"theta_wall": float(theta_w), "theta_bulk": float(theta_b), "wall_length": float(fw.sum()),
            "flow_rate": float(wx.sum() / (x_hi - x_lo)), "nusselt_b": float(1.0 / (theta_w - theta_b))}
