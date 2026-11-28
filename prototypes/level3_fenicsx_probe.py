"""FEniCSx capability probe for the Level 3 toolchain decision: 3D Taylor-Hood Stokes with a direct solve (MUMPS)."""
import resource
import sys
import time

import basix.ufl
import dolfinx
import numpy as np
import ufl
from dolfinx import fem, mesh
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI

n = int(sys.argv[1])
t0 = time.perf_counter()
msh = mesh.create_box(MPI.COMM_WORLD, [np.zeros(3), np.array([1.0, 1.0, 4.0])], [n, n, 4 * n], mesh.CellType.tetrahedron)
P2 = basix.ufl.element("Lagrange", msh.basix_cell(), 2, shape=(3,))
P1 = basix.ufl.element("Lagrange", msh.basix_cell(), 1)
W = fem.functionspace(msh, basix.ufl.mixed_element([P2, P1]))
(u, p), (v, q) = ufl.TrialFunctions(W), ufl.TestFunctions(W)
f = fem.Constant(msh, np.array([0.0, 0.0, 1.0]))
a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx - p * ufl.div(v) * ufl.dx - q * ufl.div(u) * ufl.dx
L = ufl.inner(f, v) * ufl.dx
V, _ = W.sub(0).collapse()
fdim = msh.topology.dim - 1
facets = mesh.locate_entities_boundary(msh, fdim, lambda x: np.isclose(x[0], 0) | np.isclose(x[0], 1) | np.isclose(x[1], 0) | np.isclose(x[1], 1))
u0 = fem.Function(V)
bc = fem.dirichletbc(u0, fem.locate_dofs_topological((W.sub(0), V), fdim, facets), W.sub(0))
pin = fem.locate_dofs_geometrical((W.sub(1), W.sub(1).collapse()[0]), lambda x: np.isclose(x[0], 0) & np.isclose(x[1], 0) & np.isclose(x[2], 0))
p0 = fem.Function(W.sub(1).collapse()[0])
bcp = fem.dirichletbc(p0, pin, W.sub(1))
prob = LinearProblem(a, L, bcs=[bc, bcp], petsc_options={"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps"})
wh = prob.solve()
print(f"dolfinx {dolfinx.__version__}: n={n}, dofs={W.dofmap.index_map.size_global * W.dofmap.index_map_bs}, "
      f"time {time.perf_counter() - t0:.1f} s, peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30:.2f} GB")
