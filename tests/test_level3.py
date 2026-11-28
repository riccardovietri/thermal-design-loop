"""Level 3 (3D conjugate, OpenFOAM) bookkeeping and geometry: hand calculations, no solver needed."""
import math

import pytest

from level3.foam import parse_check_mesh, parse_solver_log
from level3.geometry import (MESH_LEVELS, UnitCell, blockmesh_dict, duct_f_re_exact, grid, r_th_from_temperatures,
                             unit_cell)
from tools.types import Design, Fluid, Requirements

REQ = Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, 40.0, 5e3, "PG25")
FLUID = Fluid("PG25", 40.0, 1010.6, 3958.0, 0.481, 1.394e-3)
REC = Design(0.387e-3, 0.387e-3, 3.75e-3, 1.5)


def test_unit_cell_matches_level1_bookkeeping():
    c = unit_cell(REC, REQ, FLUID)
    assert c.n_channels == 77  # floor(60 / 0.774)
    assert c.half_fin_m == pytest.approx(0.1935e-3) and c.half_channel_m == pytest.approx(0.1935e-3)
    vdot = 1.5 / 60e3
    assert c.u_bar == pytest.approx(vdot / (77 * 0.387e-3 * 3.75e-3))  # 0.2237 m/s
    assert c.heated_flux_w_m2 == pytest.approx(1000 / (77 * 0.774e-3 * 0.060))  # over n pitches, ASSUMPTIONS M5
    assert c.coolant_rise_k == pytest.approx(1000 / (1010.6 * 3958.0 * vdot))  # about 15.0 K
    assert c.cell_heat_w(1e-3) == pytest.approx(c.heated_flux_w_m2 * 0.387e-3 * 1e-3)
    # the cell's share of the total heat equals its share of the flow
    assert c.cell_heat_w(REQ.flow_length_m) * 2 * c.n_channels == pytest.approx(1000.0)
    assert c.cell_flow_m3_s() * 2 * c.n_channels == pytest.approx(vdot)


def test_r_th_adds_half_the_coolant_rise_as_level1():
    c = unit_cell(REC, REQ, FLUID)
    r = r_th_from_temperatures(c, t_base_avg=360.0, t_bulk=351.0)
    assert r == pytest.approx((9.0 + c.coolant_rise_k / 2) / 1000.0)


def test_dp_over_flow_length():
    c = unit_cell(REC, REQ, FLUID)
    assert c.dp_core_pa(dpdx_pa_per_m=2.5e4) == pytest.approx(1500.0)


@pytest.mark.parametrize("alpha, exact", [(1.0, 14.227), (0.5, 15.548), (0.25, 18.233), (0.125, 20.585), (1e-4, 24.0)])
def test_duct_f_re_exact_series(alpha, exact):
    # Shah & London (1978) rectangular-duct values, Fanning f on the hydraulic diameter
    assert duct_f_re_exact(alpha) == pytest.approx(exact, rel=2e-4)


def test_grid_refinement_doubles_cells_and_keeps_interfaces_conformal():
    c = unit_cell(REC, REQ, FLUID)
    g0, g1 = grid(c, MESH_LEVELS[1]), grid(c, MESH_LEVELS[2])
    for k in ("n_fin", "n_channel", "n_base", "n_fluid_z", "n_x"):
        assert g1[k] == 2 * g0[k]
    assert g1["n_cells_fluid"] == 8 * g0["n_cells_fluid"] and g1["n_cells_solid"] == 8 * g0["n_cells_solid"]
    # fin and channel blocks share the fluid-height division; base blocks share the fin and channel widths
    assert g0["n_cells_solid"] == g0["n_x"] * (g0["n_fin"] * g0["n_base"] + g0["n_channel"] * g0["n_base"]
                                               + g0["n_fin"] * g0["n_fluid_z"])
    assert g0["n_cells_fluid"] == g0["n_x"] * g0["n_channel"] * g0["n_fluid_z"]


def test_blockmesh_dicts_name_the_coupled_patches():
    c = unit_cell(REC, REQ, FLUID)
    f, s = blockmesh_dict(c, "fluid", MESH_LEVELS[1]), blockmesh_dict(c, "solid", MESH_LEVELS[1])
    assert "fluid_to_solid" in f and "samplePatch solid_to_fluid" in f and "sampleRegion solid" in f
    assert "solid_to_fluid" in s and "samplePatch fluid_to_solid" in s and "sampleRegion fluid" in s
    assert f.count("hex") == 1 and s.count("hex") == 3
    for d in (f, s):
        assert "type cyclic" in d and "neighbourPatch back" in d and "neighbourPatch front" in d
    assert "bottom" in s and "cover" in f


def test_cell_length_is_whole_cells():
    c = unit_cell(REC, REQ, FLUID)
    g = grid(c, MESH_LEVELS[1])
    assert g["lx_m"] / g["n_x"] == pytest.approx(g["dx_m"]) and isinstance(c, UnitCell)
    assert math.isclose(g["lx_m"], c.periodic_length_m)


LOG = """Pressure gradient source: uncorrected Ubar = 0.2, pressure gradient = 999
Time = 1

Pressure gradient source: uncorrected Ubar = 0.21, pressure gradient = 178773.2
L3MON bulkT 313.15 flow 1.6e-07 sink 0.08
Pressure gradient source: uncorrected Ubar = 0.22, pressure gradient = 64428.4
Time = 2

Pressure gradient source: uncorrected Ubar = 0.226, pressure gradient = 22575.1
L3MON bulkT 313.02 flow 1.7e-07 sink 0.08
Pressure gradient source: uncorrected Ubar = 0.2246, pressure gradient = 48934.4
End
"""


def test_parse_solver_log_keeps_the_corrected_gradient_per_iteration():
    p = parse_solver_log(LOG)
    assert p["iterations"] == [1, 2]
    assert p["dpdx_pa_per_m"] == [64428.4, 48934.4]  # last print of each iteration; the one before Time = 1 is ignored
    assert p["bulk_t"] == [313.15, 313.02] and p["flow_m3_s"] == [1.6e-07, 1.7e-07]


def test_parse_check_mesh():
    text = ("    cells:            3456\n    Max aspect ratio = 16.3 OK.\n    Mesh non-orthogonality Max: 0 average: 0\n"
            "    Max skewness = 8.7e-14 OK.\nMesh OK.\n")
    m = parse_check_mesh(text)
    assert m == {"cells": 3456, "max_non_orthogonality": 0.0, "max_skewness": 8.7e-14, "max_aspect_ratio": 16.3, "mesh_ok": True}


def test_flatness_over_last_tenth():
    from level3.straight import flatness
    v = [10.0] * 90 + [10.02, 10.01] + [10.0] * 8
    assert flatness(v) == pytest.approx(0.002)
    assert flatness([5.0] * 100) == 0.0


def test_richardson_recovers_second_order_limit():
    from level3.straight import richardson
    r = richardson([2.0, 1.25, 1.0625])  # 1 + h^2 at h = 1, 1/2, 1/4
    assert r["order"] == pytest.approx(2.0) and r["extrapolated"] == pytest.approx(1.0)


def test_committed_level3_results_pass_their_gates():
    import json
    from pathlib import Path
    from validation.level3_straight import check
    assert check() == 0
    data = json.loads(Path("examples/output/level3_straight_validation/results.json").read_text())
    level1 = json.loads(Path("examples/output/gpu_coldplate/results.json").read_text())
    assert data["recommended"] == level1["recommended"]
