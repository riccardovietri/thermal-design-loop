import math

import numpy as np
import pytest

from level2.planforms import GeometryRejected, mesh_planform, offset_strip, pin_circle, pin_diamond, sinusoidal, straight, zigzag
from level2.planview import describe, network, run_design, superficial_velocity
from tools.fluids import fluid_props
from tools.types import Requirements

REQ = Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, 40.0, 5e3, "PG25", r_th_threshold_k_per_w=20e-3,
                   r_th_objective_k_per_w=15e-3)


def test_straight_geometry():
    pf = straight(0.4, 0.3, 3.0)
    assert pf.perimeter_density_per_mm == pytest.approx(2 / 0.7)
    assert pf.porosity == pytest.approx(0.4 / 0.7)
    assert pf.dh_mm == pytest.approx(0.8)
    assert pf.fin_section_perimeter_mm / pf.fin_section_area_mm2 == pytest.approx(2 / 0.3)


def test_offset_strip_geometry_and_gaps():
    pf = offset_strip(l=1.0, t=0.4, p_t=1.2, g=0.0, h_fin=4.0)
    # per period 2(l+g) = 2 mm and lateral half pitch 0.6 mm: one strip's perimeter and area
    assert pf.wetted_perimeter_mm == pytest.approx(2 * (1.0 + 0.4))
    assert pf.fluid_area_mm2 == pytest.approx(2.0 * 0.6 - 1.0 * 0.4)
    assert pf.perimeter_density_per_mm == pytest.approx(2.8 / 1.2)
    assert pf.min_gap_mm == pytest.approx(0.2)  # p_t/2 - t between rows
    with pytest.raises(GeometryRejected, match="min_fluid_gap"):
        offset_strip(l=1.0, t=0.4, p_t=1.0, g=0.0, h_fin=4.0)  # 0.1 mm between rows


def test_pin_geometry_and_limits():
    pf = pin_circle(d=0.4, s_t=0.7, s_l=0.6, h_fin=4.0)
    assert pf.wetted_perimeter_mm == pytest.approx(math.pi * 0.4)
    assert pf.fluid_area_mm2 == pytest.approx(0.6 * 0.7 - math.pi * 0.04)
    assert pf.min_gap_mm == pytest.approx(min(0.7 - 0.4, math.hypot(0.6, 0.35) - 0.4, 1.2 - 0.4))
    assert pf.fin_section_perimeter_mm / pf.fin_section_area_mm2 == pytest.approx(4 / 0.4)
    with pytest.raises(GeometryRejected, match="fin_aspect"):
        pin_circle(d=0.3, s_t=0.7, s_l=0.6, h_fin=4.0)  # H/d = 13.3
    with pytest.raises(GeometryRejected, match="min_fluid_gap"):
        pin_circle(d=0.4, s_t=0.5, s_l=0.6, h_fin=4.0)  # 0.1 mm between pins of a row
    with pytest.raises(GeometryRejected, match="min_solid"):
        pin_circle(d=0.1, s_t=0.7, s_l=0.6, h_fin=1.0)


def test_diamond_gap_by_hand():
    s, s_t, s_l = 0.4, 0.9, 0.7
    pf = pin_diamond(s=s, s_t=s_t, s_l=s_l, h_fin=4.0)
    a = s / math.sqrt(2)
    # diagonal neighbour at (s_l, s_t/2): facing edges are parallel lines x + y = const, distance (s_l + s_t/2 - 2a)/sqrt 2
    diag = (s_l + s_t / 2 - 2 * a) / math.sqrt(2)
    assert pf.min_gap_mm == pytest.approx(min(s_t - 2 * a, diag), abs=1e-9)
    assert pf.wetted_perimeter_mm == pytest.approx(4 * s)


def test_zigzag_geometry_by_hand():
    w, t, amp, lam, r_in = 0.4, 0.4, 0.3, 2.0, 0.05
    pf = zigzag(w, t, amp, lam, 4.0)
    phi = math.atan(4 * amp / lam)
    r = r_in + w / 2
    seg = math.hypot(lam / 2, 2 * amp)
    lc = 2 * (seg - 2 * r * math.tan(phi)) + 2 * r * 2 * phi
    assert pf.lateral_unit_mm == pytest.approx((w + t) / math.cos(phi))
    assert pf.wetted_perimeter_mm == pytest.approx(2 * lc)
    assert pf.fluid_area_mm2 == pytest.approx(w * lc)


def test_sinusoid_fin_thickness_and_fold_rejection():
    w, t, amp, lam = 0.4, 0.4, 0.3, 2.0
    pf = sinusoidal(w, t, amp, lam, 4.0)
    k = 2 * math.pi / lam
    x = np.linspace(0, lam, 20001)
    lc = np.sum(np.hypot(np.diff(x), np.diff(amp * np.cos(k * x))))
    assert pf.wetted_perimeter_mm == pytest.approx(2 * lc, rel=1e-4)
    # the lateral pitch exceeds w + t (fins are thicker at the crests) but not the sharp-wave value
    assert w + t < pf.lateral_unit_mm < (w + t) / math.cos(math.atan(amp * k))
    with pytest.raises(GeometryRejected, match="fold"):
        sinusoidal(w, t, 0.6, 2.0, 4.0)  # crest radius 0.17 mm < w/2 + 0.05


def test_network_hand_calculation():
    fluid = fluid_props("PG25", 40.0)
    pf = straight(0.4, 0.4, 4.0)
    h, dpdx = 5000.0, 20000.0
    net = network(pf, dpdx, h, REQ, fluid, 1.5)
    a_fp = 0.06 * 0.06
    m = math.sqrt(2 * h / (390.0 * 0.4e-3))
    eta = math.tanh(m * 4e-3) / (m * 4e-3)
    a_floor, a_fin = 0.5 * a_fp, 2 / 0.8e-3 * a_fp * 4e-3
    r = 0.002 / (390 * a_fp) + 1 / (h * (a_floor + eta * a_fin)) + 1 / (2 * fluid.rho * 2.5e-5 * fluid.cp)
    assert net["fin_efficiency"] == pytest.approx(eta)
    assert net["r_th_k_per_w"] == pytest.approx(r)
    assert net["dp_pa"] == pytest.approx(dpdx * 0.06)


def test_coarse_straight_periodic_solve():
    fluid = fluid_props("PG25", 40.0)
    pf = straight(0.4, 0.4, 4.0)
    r = run_design(pf, REQ, fluid, 1.5, cells=4)
    u_ch = superficial_velocity(REQ, 1.5, 4.0) / pf.porosity
    assert r.dpdx_pa_per_m == pytest.approx(12 * fluid.mu * u_ch / 0.4e-3**2, rel=1e-6)
    assert r.h_w_m2k * 0.8e-3 / fluid.k == pytest.approx(140 / 17, rel=1e-3)
    assert abs(r.energy_balance_error) < 1e-6


def test_periodic_mesh_heated_length():
    pm = mesh_planform(pin_circle(d=0.4, s_t=0.7, s_l=0.6, h_fin=4.0), n_periods=2, cells=4)
    f = pm.mesh.facets[:, pm.heated]
    length = np.sum(np.linalg.norm(pm.mesh.p[:, f[0]] - pm.mesh.p[:, f[1]], axis=0))
    assert length == pytest.approx(2 * math.pi * 0.4, rel=2e-3)  # polygonal circle
    assert pm.cuts.tolist() == pytest.approx([0.0, 1.2, 2.4])


def test_describe_uses_plain_words():
    assert describe("pin_circle", {"d": 0.4, "s_t": 0.8, "s_l": 0.9, "h_fin": 4.0}) == "Circular pins, 0.4 mm diameter, 0.8 × 0.9 mm pitch"
    assert describe("straight", {"w": 0.4, "t": 0.4, "h_fin": 4.0}) == "Straight fins, 0.4 mm channels, 0.4 mm fins"
