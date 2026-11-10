"""Toolchain validation: fully developed laminar flow + H1 heat transfer in a rectangular duct.

Geometry/mesh: gmsh. Solve: scikit-fem (P2 triangles). Reference: Shah & London (1978),
Fanning f*Re and Nu_H1 for rectangular ducts.

Nondimensional problem on the duct cross-section (width a, height b):
    -lap(u) = 1,           u = 0 on walls    ->  f*Re (Fanning) = Dh^2 / (2 u_mean)
     lap(theta) = u/u_mean, theta = 0 on walls ->  Nu_H1 = Dh^2 / (4 * -theta_bulk)
"""
import json
import sys
from pathlib import Path

import gmsh
import numpy as np
from skfem import Basis, ElementTriP2, MeshTri, asm, condense, solve
from skfem.helpers import dot, grad
from skfem import BilinearForm, LinearForm

# alpha = short/long side : (Fanning f*Re, Nu_H1), Shah & London 1978
REFERENCE = {1.0: (14.227, 3.608), 0.5: (15.548, 4.123), 0.25: (18.233, 5.331), 0.125: (20.585, 6.490)}


def rect_mesh(a: float, b: float, h: float) -> MeshTri:
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.occ.addRectangle(0, 0, 0, a, b)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.model.mesh.generate(2)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    _, _, conn = gmsh.model.mesh.getElements(2)
    gmsh.finalize()
    idx = np.zeros(int(tags.max()) + 1, dtype=np.int64)
    idx[tags.astype(np.int64)] = np.arange(len(tags))
    p = coords.reshape(-1, 3)[:, :2].T
    t = idx[conn[0].astype(np.int64)].reshape(-1, 3).T
    return MeshTri(np.ascontiguousarray(p), np.ascontiguousarray(t))


@BilinearForm
def laplace(u, v, _):
    return dot(grad(u), grad(v))


@BilinearForm
def mass(u, v, _):
    return u * v


@LinearForm
def unit(v, _):
    return 1.0 * v


def solve_duct(alpha: float, h: float) -> dict:
    a, b = 1.0, alpha
    mesh = rect_mesh(a, b, h)
    basis = Basis(mesh, ElementTriP2())
    K, M, f = asm(laplace, basis), asm(mass, basis), asm(unit, basis)
    D = basis.get_dofs()
    u = solve(*condense(K, f, D=D))
    area = a * b
    u_mean = (f @ u) / area
    dh = 4 * area / (2 * (a + b))
    f_re = dh**2 / (2 * u_mean)
    theta = solve(*condense(K, -(M @ u) / u_mean, D=D))
    theta_bulk = (u @ (M @ theta)) / (u_mean * area)
    nu = dh**2 / (4 * -theta_bulk)
    return {"alpha": alpha, "h": h, "n_dofs": int(basis.N), "f_re": float(f_re), "nu_h1": float(nu)}


def main() -> int:
    rows, worst = [], 0.0
    for alpha, (f_ref, nu_ref) in REFERENCE.items():
        for h in (alpha / 8, alpha / 16, alpha / 32):
            r = solve_duct(alpha, h)
            r["f_re_err_pct"] = 100 * (r["f_re"] / f_ref - 1)
            r["nu_err_pct"] = 100 * (r["nu_h1"] / nu_ref - 1)
            rows.append(r)
        worst = max(worst, abs(rows[-1]["f_re_err_pct"]), abs(rows[-1]["nu_err_pct"]))
    print(f"{'alpha':>6} {'h':>8} {'dofs':>7} {'f*Re':>8} {'err%':>7} {'Nu_H1':>7} {'err%':>7}")
    for r in rows:
        print(f"{r['alpha']:6.3f} {r['h']:8.4f} {r['n_dofs']:7d} {r['f_re']:8.4f} {r['f_re_err_pct']:7.3f} {r['nu_h1']:7.4f} {r['nu_err_pct']:7.3f}")
    Path(__file__).with_name("duct_fully_developed.json").write_text(json.dumps(rows, indent=2))
    print(f"worst finest-mesh error: {worst:.3f}% (gate: 1%)")
    return 0 if worst < 1.0 else 1


if __name__ == "__main__":
    sys.exit(main())
