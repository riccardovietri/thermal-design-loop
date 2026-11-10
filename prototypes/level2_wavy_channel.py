"""PROTOTYPE (not part of V1): Level 2 = 2D planar steady Navier-Stokes, developing laminar flow.

Plan view of one fin gap (infinite fin height). Straight vs sine-wave walls at two Reynolds numbers.
Nondimensional: lengths / gap w, velocities / mean velocity U, pressure / (rho U^2), Re_w = rho U w / mu.
Taylor-Hood (P2-P1), Newton with Reynolds continuation. Isothermal: flow field and pressure drop only.
"""
import sys
from pathlib import Path

import gmsh
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.tri import Triangulation
from scipy.sparse import bmat
from skfem import Basis, BilinearForm, ElementTriP1, ElementTriP2, ElementVector, LinearForm, MeshTri, asm, condense, solve
from skfem.helpers import ddot, div, dot, grad, mul

LEAD_IN, N_WAVES = 2.0, 5


def centerline(x, amp, lam):
    s = np.clip((x - LEAD_IN) / (N_WAVES * lam), 0, 1)
    return amp * np.sin(2 * np.pi * N_WAVES * s)


def make_mesh(amp: float, lam: float, size: float) -> tuple[MeshTri, float]:
    length = 2 * LEAD_IN + N_WAVES * lam
    xs = np.linspace(0, length, int(length / size) + 1)
    yc = centerline(xs, amp, lam)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    geo = gmsh.model.geo
    bot = [geo.addPoint(x, y - 0.5, 0, size) for x, y in zip(xs, yc)]
    top = [geo.addPoint(x, y + 0.5, 0, size) for x, y in zip(xs, yc)]
    lines = [geo.addSpline(bot), geo.addLine(bot[-1], top[-1]), geo.addSpline(top[::-1]), geo.addLine(top[0], bot[0])]
    geo.addPlaneSurface([geo.addCurveLoop(lines)])
    geo.synchronize()
    gmsh.model.mesh.generate(2)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    _, _, conn = gmsh.model.mesh.getElements(2)
    gmsh.finalize()
    idx = np.zeros(int(tags.max()) + 1, dtype=np.int64)
    idx[tags.astype(np.int64)] = np.arange(len(tags))
    used = np.unique(conn[0].astype(np.int64))
    remap = np.full(int(tags.max()) + 1, -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    p = coords.reshape(-1, 3)[idx[used], :2].T
    t = remap[conn[0].astype(np.int64)].reshape(-1, 3).T
    return MeshTri(np.ascontiguousarray(p), np.ascontiguousarray(t)), length


@BilinearForm
def viscous(u, v, _):
    return ddot(grad(u), grad(v))


@BilinearForm
def divergence(u, q, _):
    return -div(u) * q


@BilinearForm
def advection_jac(u, v, w):
    z = w["wind"]
    return dot(mul(grad(u), z), v) + dot(mul(grad(z), u), v)


@LinearForm
def advection_res(v, w):
    z = w["wind"]
    return dot(mul(grad(z), z), v)


@BilinearForm
def laplace(u, v, _):
    return dot(grad(u), grad(v))


@LinearForm
def vorticity(v, w):
    z = w["wind"]
    return (z.grad[1][0] - z.grad[0][1]) * v


def solve_channel(amp: float, lam: float, reynolds_w: list[float], size: float = 0.07) -> list[dict]:
    mesh, length = make_mesh(amp, lam, size)
    bu = Basis(mesh, ElementVector(ElementTriP2()), intorder=4)
    bp = bu.with_element(ElementTriP1())
    tol = 1e-8
    inlet = mesh.facets_satisfying(lambda x: x[0] < tol, boundaries_only=True)
    outlet = mesh.facets_satisfying(lambda x: x[0] > length - tol, boundaries_only=True)
    walls = np.setdiff1d(mesh.boundary_facets(), np.concatenate([inlet, outlet]))
    d_u = np.concatenate([bu.get_dofs(facets=inlet).all(), bu.get_dofs(facets=walls).all()])
    nu_dofs = bu.N

    def inflow(x):
        eta = x[1] + 0.5
        return np.stack([6 * eta * (1 - eta), np.zeros_like(eta)])

    up = np.zeros(nu_dofs + bp.N)
    inlet_dofs = bu.get_dofs(facets=inlet).all()
    up[inlet_dofs] = bu.project(inflow)[inlet_dofs] if False else _nodal_inflow(bu, inlet_dofs, inflow)
    a_visc, b = asm(viscous, bu), asm(divergence, bu, bp)

    out = []
    first = True
    for re in reynolds_w:
        for it in range(25):
            wind = bu.interpolate(up[:nu_dofs])
            jac = bmat([[a_visc / re + asm(advection_jac, bu, wind=wind), b.T], [b, None]], format="csr")
            res = np.concatenate([
                a_visc @ up[:nu_dofs] / re + asm(advection_res, bu, wind=wind) + b.T @ up[nu_dofs:],
                b @ up[:nu_dofs],
            ])
            delta = solve(*condense(jac, -res, D=d_u))
            up += delta
            norm = np.linalg.norm(delta) / max(np.linalg.norm(up), 1e-30)
            if norm < 1e-9:
                break
        else:
            raise RuntimeError(f"Newton did not converge at Re_w={re} (last step {norm:.2e})")
        # pressure drop across the wavy section, compared with fully developed plane Poiseuille (12/Re per unit length)
        bp_in = bp.get_dofs(facets=inlet).all()
        p = up[nu_dofs:]
        dp = p[bp_in].mean() - p[bp.get_dofs(facets=outlet).all()].mean()
        dp_straight_1d = 12.0 / re * length
        # stream function: -lap(psi) = vorticity, psi = 0 / 1 on bottom / top wall
        b1 = bp
        wind = bu.interpolate(up[:nu_dofs])
        wall_dofs = b1.get_dofs(facets=walls).all()
        psi0 = np.zeros(b1.N)
        xw, yw = b1.doflocs[0, wall_dofs], b1.doflocs[1, wall_dofs]
        psi0[wall_dofs] = (yw > centerline(xw, amp, lam)).astype(float)
        psi = solve(*condense(asm(laplace, b1), asm(vorticity, b1, wind=wind), x=psi0, D=wall_dofs))
        n = mesh.p.shape[1]
        speed = np.hypot(up[bu.nodal_dofs[0]], up[bu.nodal_dofs[1]])
        out.append({"re_w": re, "re_dh": 2 * re, "newton_its": it + 1, "dp": float(dp), "dp_ratio_vs_1d": float(dp / dp_straight_1d),
                    "psi_min": float(psi.min()), "psi_max": float(psi.max()), "mesh": mesh, "speed": speed, "psi": psi[:n],
                    "length": length, "amp": amp, "lam": lam})
    return out


def _nodal_inflow(bu, dofs, fn):
    vals = fn(bu.doflocs[:, dofs])
    comp = np.isin(dofs, np.concatenate([bu.nodal_dofs[1], bu.facet_dofs[1]])).astype(int)
    return vals[comp, np.arange(len(dofs))]


def main() -> int:
    gap_mm = 0.3
    cases = [("Straight fins", 0.0, 2.0 / gap_mm), ("Sine-wave fins (amplitude 0.3 mm, wavelength 2 mm)", 1.0, 2.0 / gap_mm)]
    re_list = [10.0, 25.0, 50.0, 100.0]  # Re_w; Re on Dh = 2w is double
    show = [25.0, 100.0]
    fig, axes = plt.subplots(4, 1, figsize=(11, 9.5), constrained_layout=True)
    row = 0
    for title, amp, lam in cases:
        results = solve_channel(amp, lam, re_list)
        for r in results:
            print(f"{title[:14]:14s} Re_Dh={r['re_dh']:5.0f} newton={r['newton_its']:2d} dp/(rho U^2)={r['dp']:7.3f} "
                  f"dp / fully-developed straight 1D = {r['dp_ratio_vs_1d']:.3f}  psi range [{r['psi_min']:.4f}, {r['psi_max']:.4f}]")
        for r in [r for r in results if r["re_w"] in show]:
            ax = axes[row]
            row += 1
            m = r["mesh"]
            tri = Triangulation(m.p[0] * gap_mm, m.p[1] * gap_mm, m.t.T)
            cs = ax.tricontourf(tri, r["speed"], levels=np.linspace(0, max(1.6, r["speed"].max()), 30), cmap="viridis")
            ax.tricontour(tri, r["psi"], levels=np.linspace(0.05, 0.95, 10), colors="white", linewidths=0.5)
            ax.tricontour(tri, r["psi"], levels=[-0.02, -0.005, -0.0005, 1.0005, 1.005, 1.02], colors="#ff6b6b", linewidths=0.9)
            ax.set_aspect("equal")
            ax.set_xlim(LEAD_IN * gap_mm + 2.0, LEAD_IN * gap_mm + 8.0)
            ax.set_ylabel("mm")
            ax.set_title(f"{title}, Re = {r['re_dh']:.0f}:  ΔP = {r['dp_ratio_vs_1d']:.2f} × straight-channel 1D prediction", fontsize=10)
            fig.colorbar(cs, ax=ax, label="speed / mean", shrink=0.9, ticks=[0, 0.5, 1.0, 1.5, 2.0])
    axes[-1].set_xlabel("flow direction (mm); white = streamlines, red = recirculating streamlines")
    out = Path(__file__).parent / "level2_wavy_channel.png"
    fig.savefig(out, dpi=180)
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
