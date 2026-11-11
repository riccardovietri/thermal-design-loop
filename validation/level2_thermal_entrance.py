"""Level 2 validation: thermally developing laminar flow between parallel plates, uniform wall heat flux.

Checks the planar Navier-Stokes + energy solver against
  (a) plane Poiseuille pressure gradient, 12/Re per unit length (exact),
  (b) fully developed Nu = 8.235 on Dh = 2 x gap (exact, 140/17),
  (c) local Nu(x*) from an independent finite-difference marching solution of the Graetz problem,
  (d) the overall energy balance.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import solve_banded
from skfem import MeshTri

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from level2.planar import boundary_profile, solve_energy, solve_flow  # noqa: E402

RE, PR, LENGTH = 20.0, 11.5, 100.0  # Re on the gap; Pe = 230; x* = x / (4 Pe) reaches 0.109


def graetz_marching(pe: float, x_out: np.ndarray, ny: int = 801) -> np.ndarray:
    """Crank-Nicolson march of u(y) dtheta/dx = (1/Pe) d2theta/dy2 with unit wall flux; returns Nu on Dh at x_out."""
    y = np.linspace(0, 1, ny)
    dy = y[1] - y[0]
    u = 6 * y * (1 - y)
    u[0] = u[-1] = 1.0  # placeholder; boundary rows are replaced by the flux condition
    theta, x, out = np.zeros(ny), 0.0, []
    xs = np.unique(np.concatenate([np.geomspace(1e-5, x_out.max(), 4000), x_out]))
    for x_next in xs:
        dx = x_next - x
        r = dx / (2 * pe * dy**2) / u
        ab = np.zeros((3, ny))
        ab[0, 1:], ab[1, :], ab[2, :-1] = -r[:-1], 1 + 2 * r, -r[1:]
        rhs = theta.copy()
        rhs[1:-1] += r[1:-1] * (theta[2:] - 2 * theta[1:-1] + theta[:-2])
        # second-order one-sided flux condition: dtheta/dn = 1 at both walls
        ab[1, 0], ab[0, 1] = 1.0, -1.0
        ab[1, -1], ab[2, -2] = 1.0, -1.0
        rhs[0] = rhs[-1] = dy
        theta = solve_banded((1, 1), ab, rhs)
        # refine the wall value with a second-order extrapolation consistent with the flux
        wall = (4 * theta[1] - theta[2] + 2 * dy) / 3
        x = x_next
        if np.any(np.isclose(x_out, x)):
            out.append(2.0 / (wall - 2 * x / pe))
    return np.array(out)


def main() -> int:
    xg = np.concatenate([[0.0], np.geomspace(0.01, LENGTH, 400)])
    s = np.linspace(0, 1, 41)
    yg = 0.5 * (1 - np.cos(np.pi * s))  # clustered at both walls
    mesh = MeshTri.init_tensor(xg, yg)
    tol = 1e-9
    inlet = mesh.facets_satisfying(lambda x: x[0] < tol, boundaries_only=True)
    outlet = mesh.facets_satisfying(lambda x: x[0] > LENGTH - tol, boundaries_only=True)
    walls = mesh.facets_satisfying(lambda x: (x[1] < tol) | (x[1] > 1 - tol), boundaries_only=True)
    bottom = mesh.facets_satisfying(lambda x: x[1] < tol, boundaries_only=True)
    parabola = lambda x: 6 * x[1] * (1 - x[1])
    flow = solve_flow(mesh, [RE], {0: [(inlet, parabola), (walls, None)], 1: [(inlet, None), (walls, None)]})
    energy = solve_energy(flow, PR, inlet, walls, outlet)
    pe = energy.peclet

    p_in = flow.p[flow.basis_p.get_dofs(facets=inlet).all()].mean()
    p_out = flow.p[flow.basis_p.get_dofs(facets=outlet).all()].mean()
    dp_err = (p_in - p_out) / (12 / RE * LENGTH) - 1

    x, theta_w = boundary_profile(energy, bottom)
    x_star = x / (4 * pe)
    nu = 2.0 / (theta_w - 2 * x / pe)
    probes = np.array([0.002, 0.005, 0.01, 0.02, 0.05, 0.1])
    nu_fe = np.interp(probes, x_star[1:], nu[1:])
    nu_ref = graetz_marching(pe, probes * 4 * pe)
    nu_fd_err = nu_fe[-1] / 8.235 - 1

    print(f"Re = {RE:.0f}, Pr = {PR}, Pe = {pe:.0f}; {flow.basis_u.N + flow.basis_p.N + energy.basis.N} unknowns; Newton iterations {flow.newton_iterations}")
    print(f"(a) pressure drop vs plane Poiseuille: {100*dp_err:+.4f}%")
    print(f"(b) Nu at x* = {probes[-1]}: {nu_fe[-1]:.4f} vs fully developed 8.235 ({100*nu_fd_err:+.3f}%)")
    print("(c) local Nu vs independent marching solution:")
    print(f"{'x*':>8} {'Nu FE':>9} {'Nu ref':>9} {'diff %':>8}")
    worst = 0.0
    for xs, a, b in zip(probes, nu_fe, nu_ref):
        worst = max(worst, abs(a / b - 1))
        print(f"{xs:8.3f} {a:9.4f} {b:9.4f} {100*(a/b-1):8.3f}")
    print(f"(d) energy balance error: {100*energy.energy_balance_error:+.4f}%")
    ok = abs(dp_err) < 0.01 and abs(nu_fd_err) < 0.01 and worst < 0.02 and abs(energy.energy_balance_error) < 0.01
    Path(__file__).with_name("level2_thermal_entrance.json").write_text(json.dumps({
        "dp_error": dp_err, "nu_fully_developed_error": nu_fd_err, "nu_local_worst_error": worst,
        "energy_balance_error": energy.energy_balance_error,
        "table": [{"x_star": float(a), "nu_fe": float(b), "nu_ref": float(c)} for a, b, c in zip(probes, nu_fe, nu_ref)]}, indent=2))
    print("PASS" if ok else "FAIL", "(gates: 1% on pressure, fully developed Nu and energy; 2% on local Nu)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
