import numpy as np
import pytest
from skfem import MeshTri

from level2.planar import boundary_profile, solve_energy, solve_flow
from level2.slot_jet import SlotJetDesign, run_slot_jet
from tools.fluids import fluid_props
from tools.types import Requirements


def test_plane_poiseuille_pressure_and_energy_balance():
    length, re = 10.0, 20.0
    mesh = MeshTri.init_tensor(np.linspace(0, length, 41), np.linspace(0, 1, 13))
    tol = 1e-9
    inlet = mesh.facets_satisfying(lambda x: x[0] < tol, boundaries_only=True)
    outlet = mesh.facets_satisfying(lambda x: x[0] > length - tol, boundaries_only=True)
    walls = mesh.facets_satisfying(lambda x: (x[1] < tol) | (x[1] > 1 - tol), boundaries_only=True)
    flow = solve_flow(mesh, [re], {0: [(inlet, lambda x: 6 * x[1] * (1 - x[1])), (walls, None)], 1: [(inlet, None), (walls, None)]})
    dp = flow.p[flow.basis_p.get_dofs(facets=inlet).all()].mean() - flow.p[flow.basis_p.get_dofs(facets=outlet).all()].mean()
    assert dp == pytest.approx(12 / re * length, rel=1e-6)
    energy = solve_energy(flow, 5.0, inlet, walls, outlet)
    assert abs(energy.energy_balance_error) < 5e-3
    x, theta = boundary_profile(energy, walls)
    assert theta.min() >= -1e-6 and theta[-1] > theta[1]


def test_slot_jet_small_case(tmp_path):
    req = Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, 40.0, 30e3, "PG25")
    d = SlotJetDesign(0.3e-3, 0.15e-3, 2.0e-3, 1.5)
    r = run_slot_jet(d, req, fluid_props("PG25", 40.0), tmp_path, cells_across=8)
    assert r.n_slots == 30
    assert r.reynolds_slot == pytest.approx(1010.6 * (2.5e-5 / (30 * 0.3e-3 * 0.06)) * 0.3e-3 / 1.394e-3, rel=1e-3)
    assert abs(r.energy_balance_error) < 0.01
    assert r.h_stagnation_w_m2k > r.h_mean_w_m2k > 0
    assert r.dp_pa > 0
    assert (tmp_path / "fields.png").stat().st_size > 0
