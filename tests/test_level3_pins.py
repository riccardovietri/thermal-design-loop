"""Level 3 stage two (staggered circular pins, OpenFOAM): bookkeeping and plan geometry from hand calculations, no solver."""
import math

import pytest

from level3.pins import (P0_SLAB_MM, PinCell, layer_heights, pin_cell, pin_centres, plan_mesh_stats, r_th_pins,
                         planview_h)
from tools.types import Fluid, Requirements

REQ = Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, 40.0, 5e3, "PG25")
FLUID = Fluid("PG25", 40.0, 1010.6, 3958.0, 0.481, 1.394e-3)
PICK = {"d": 0.4, "s_t": 0.8, "s_l": 0.9, "h_fin": 4.0}


def test_pick_bookkeeping_matches_the_plan_view_definitions():
    c = pin_cell(PICK, REQ, FLUID)
    assert isinstance(c, PinCell)
    assert c.design_id == "pin_circle_d0400_s_t0800_s_l0900_h_fin4000"  # level2.planview.design_id
    assert c.lateral_m == pytest.approx(0.4e-3) and c.period_m == pytest.approx(1.8e-3)  # half pitch, one stagger period
    us = 2.5e-5 / (0.060 * 0.004)  # superficial velocity: total flow over 60 mm width and fin height
    assert c.u_superficial == pytest.approx(us)  # 0.1042 m/s
    assert c.porosity == pytest.approx(1 - math.pi * 0.2**2 / 0.72)  # one pin (two halves) in 1.8 x 0.4 mm
    assert c.cell_flow_m3_s == pytest.approx(us * 0.004 * 0.4e-3)
    assert c.cell_flow_m3_s * (60 / 0.4) == pytest.approx(2.5e-5)  # 150 half-pitch strips carry the total flow
    assert c.heated_flux_w_m2 == pytest.approx(1000 / 0.0036)  # 1 kW over the whole 60 x 60 mm footprint
    assert c.cell_heat_w == pytest.approx(1000 / 0.0036 * 0.72e-6)  # 0.2 W
    assert c.cell_heat_w * (60 / 0.4) * (60 / 1.8) == pytest.approx(1000.0)
    assert c.coolant_rise_k == pytest.approx(1000 / (1010.6 * 3958.0 * 2.5e-5))
    # plan view: Re_d,max = rho U_s S_T / (S_T - d) d / mu
    assert c.re_d_max == pytest.approx(1010.6 * us * 2 * 0.4e-3 / 1.394e-3)  # 60.4
    assert c.re_d_max == pytest.approx(60.4, abs=0.1)


def test_full_pitch_cell_doubles_the_lateral_width_and_keeps_flux_and_velocity():
    h, f = pin_cell(PICK, REQ, FLUID), pin_cell(PICK, REQ, FLUID, full_pitch=True)
    assert f.lateral_m == pytest.approx(2 * h.lateral_m)
    assert f.porosity == pytest.approx(h.porosity)
    assert f.cell_heat_w == pytest.approx(2 * h.cell_heat_w) and f.cell_flow_m3_s == pytest.approx(2 * h.cell_flow_m3_s)
    assert f.u_superficial == pytest.approx(h.u_superficial)


def test_mean_velocity_target_uses_the_meshed_fluid_area():
    c = pin_cell(PICK, REQ, FLUID)
    exact_fluid = 0.72e-6 - math.pi * 0.2e-3**2
    # meanVelocityForce holds the volume average of u_x over the fluid: Q_cell * L / V_fluid = U_s / porosity
    assert c.ubar(exact_fluid) == pytest.approx(c.u_superficial / c.porosity)
    assert c.ubar(0.99 * exact_fluid) == pytest.approx(c.u_superficial / c.porosity / 0.99)


def test_pin_centres_symmetric_and_full_pitch():
    # symmetric cell x in [-S_L/2, 3 S_L/2], y in [0, S_T/2]: half pin A on y = 0 at x = 0, half pin B on y = S_T/2 at x = S_L
    assert pin_centres(PICK, full_pitch=False) == [(0.0, 0.0, "up"), (0.9, 0.4, "down")]
    # full pitch y in [-S_T/2, S_T/2]: whole pin A, two half pins B on the lateral boundaries
    assert pin_centres(PICK, full_pitch=True) == [(0.0, 0.0, "full"), (0.9, -0.4, "up"), (0.9, 0.4, "down")]


def test_r_th_adds_half_the_coolant_rise_as_stage_one():
    c = pin_cell(PICK, REQ, FLUID)
    assert r_th_pins(c, t_base_avg=330.0, t_bulk=322.0) == pytest.approx((8.0 + c.coolant_rise_k / 2) / 1000.0)


def test_planview_h_definition():
    # h = q'' / (period-mean wall T - period-mean velocity-weighted bulk T)
    assert planview_h(flux_w_m2=5e4, t_wall=320.0, t_bulk=315.0) == pytest.approx(1e4)


def test_layer_heights_are_graded_symmetric_and_sum_to_one():
    h = layer_heights(8, 4.0, symmetric=True)
    assert len(h) == 8 and h[-1] == pytest.approx(1.0)
    dz = [b - a for a, b in zip([0.0] + h[:-1], h)]
    assert dz[3] / dz[0] == pytest.approx(4.0) and dz[4] / dz[7] == pytest.approx(4.0)  # last / first in each half
    assert dz[0] == pytest.approx(dz[7])
    one = layer_heights(6, 0.25, symmetric=False)  # base: fine at the top (interface), coarse at the chip side... ratio last/first
    dzb = [b - a for a, b in zip([0.0] + one[:-1], one)]
    assert dzb[-1] / dzb[0] == pytest.approx(0.25)
    assert layer_heights(1, 1.0, symmetric=False) == [1.0]


def test_plan_mesh_is_periodic_and_measures_the_fluid_area():
    s = plan_mesh_stats(PICK, wall_mm=0.02, core_mm=0.05, full_pitch=False)
    exact_fluid = 0.72 - math.pi * 0.2**2
    assert s["fluid_area_mm2"] == pytest.approx(exact_fluid, rel=2e-3)  # polygonal arcs: slightly more fluid
    assert s["fluid_area_mm2"] + s["pin_area_mm2"] == pytest.approx(0.72, rel=1e-9)
    assert s["pin_perimeter_mm"] == pytest.approx(math.pi * 0.4, rel=2e-3)
    assert s["periodic_nodes_match"] and s["n_plan_cells"] > 100


def test_p0_slab_is_thin():
    assert 0 < P0_SLAB_MM <= 0.1


def test_designs_are_the_pick_and_its_front_neighbours():
    import json
    from pathlib import Path
    from level3.pins_study import pin_designs, straight_planview_pick
    pv = json.loads(Path("examples/output/level2_planview/results.json").read_text())
    ids = [r["design_id"] for r in pin_designs(pv)]
    assert ids == ["pin_circle_d0400_s_t0900_s_l0700_h_fin4000", "pin_circle_d0400_s_t0800_s_l0900_h_fin4000",
                   "pin_circle_d0400_s_t0800_s_l0700_h_fin4000"]
    assert straight_planview_pick(pv)["params_mm"] == {"w": 0.4, "t": 0.4, "h_fin": 4.0}


def test_loglog_interpolation_and_straight_front():
    import json
    from pathlib import Path
    from level3.pins_study import interp_loglog, straight_front
    assert interp_loglog(2.0, [1.0, 4.0], [10.0, 2.5]) == pytest.approx(5.0)  # power law R = 10 / dP
    assert interp_loglog(5.0, [1.0, 4.0], [10.0, 2.5]) is None  # no extrapolation
    fx, fy = straight_front(json.loads(Path("examples/output/gpu_coldplate/results.json").read_text()))
    assert all(a < b for a, b in zip(fx, fx[1:])) and all(a > b for a, b in zip(fy, fy[1:]))  # a Pareto front


def test_patch_layout_symmetric_and_full_pitch():
    from level3.pins import patch_layout
    h, f = pin_cell(PICK, REQ, FLUID), pin_cell(PICK, REQ, FLUID, full_pitch=True)
    sym = {n: i for n, i, _ in patch_layout(h, "conjugate")["fluid"]}
    assert "symmetryPlane" in sym["symLo"] and "cyclic" in sym["front"] and "separationVector (0.0018 0 0)" in sym["front"]
    full = {n: i for n, i, _ in patch_layout(f, "conjugate")["solid"]}
    assert "neighbourPatch sideHi" in full["sideLo"] and "(0 0.0008 0)" in full["sideLo"] and "symLo" not in full
    p0 = {n: i for n, i, _ in patch_layout(h, "p0")["fluid"]}
    assert "symmetryPlane" in p0["floor"] and "symmetryPlane" in p0["cover"] and "wall" in p0["pins"]


def test_conjugate_mesh_extrusion_counts(tmp_path):
    from level3.pins import write_mesh
    m = write_mesh(tmp_path / "mesh.msh", PICK, 1, "conjugate")
    assert m["patches"] == ["backRaw", "bottomRaw", "coverRaw", "frontRaw", "hiRaw", "loRaw", "pinTopRaw"]
    # every plan cell (fluid and pin) is extruded through the base and the fin layer
    n_plan_total = m["n_cells"] / (m["n_base_layers"] + m["n_fluid_layers"])
    assert n_plan_total == int(n_plan_total) and n_plan_total > m["n_plan_cells"]
    assert m["periodic_nodes_match"]


def test_committed_stage_two_results_pass_their_recorded_gates():
    from validation.level3_pins import check
    assert check() == 0
