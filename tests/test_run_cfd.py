import math
from dataclasses import replace

import pytest

from tools.fluids import fluid_props
from tools.mesh import mesh
from tools.run_cfd import run_cfd
from tools.types import Design, Requirements

REQ = Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, 40.0, 30e3, "PG25")


def shah_london_f_re(alpha: float) -> float:
    return 24 * (1 - 1.3553 * alpha + 1.9467 * alpha**2 - 1.7012 * alpha**3 + 0.9564 * alpha**4 - 0.2537 * alpha**5)


@pytest.fixture(scope="module")
def pg25():
    return fluid_props("PG25", 40.0)


def test_friction_matches_shah_london(tmp_path, pg25):
    d = Design(0.5e-3, 0.5e-3, 2.0e-3, 1.5)
    r = run_cfd(d, REQ, pg25, mesh(d, REQ.base_thickness_m, tmp_path).path)
    assert r.f_re_fanning == pytest.approx(shah_london_f_re(0.25), rel=2e-3)
    assert r.n_channels == 60
    assert r.energy_residual_rel < 1e-8


def test_dp_scales_linearly_with_flow(tmp_path, pg25):
    d = Design(0.5e-3, 0.5e-3, 2.0e-3, 1.0)
    path = mesh(d, REQ.base_thickness_m, tmp_path).path
    r1, r2 = run_cfd(d, REQ, pg25, path), run_cfd(replace(d, flow_lpm=2.0), REQ, pg25, path)
    assert r2.dp_pa == pytest.approx(2 * r1.dp_pa, rel=1e-9)


def test_tall_channel_isothermal_fin_approaches_parallel_plate_nu(tmp_path, pg25):
    """h/w = 40 with a near-isothermal solid: Nu on Dh tends to the parallel-plate H value 8.235."""
    d = Design(0.1e-3, 0.4e-3, 4.0e-3, 1.5)
    req = replace(REQ, solid_k=1e7)
    r = run_cfd(d, req, pg25, mesh(d, req.base_thickness_m, tmp_path, cells_across_half_channel=4).path)
    assert r.nu_eff == pytest.approx(8.235, rel=0.05)


def test_base_conduction_and_caloric_terms(tmp_path, pg25):
    """R_th must exceed the exact floor: base conduction + half the coolant rise."""
    d = Design(0.5e-3, 0.5e-3, 2.0e-3, 1.5)
    r = run_cfd(d, REQ, pg25, mesh(d, REQ.base_thickness_m, tmp_path).path)
    caloric = 1 / (2 * pg25.rho * pg25.cp * 1.5 / 60e3)
    r_base = REQ.base_thickness_m / (REQ.solid_k * 0.060 * 0.060)
    assert r.coolant_rise_k == pytest.approx(1000 * 2 * caloric, rel=1e-9)
    assert r.r_th_k_per_w > r_base + caloric
    assert r.r_th_max_k_per_w > r.r_th_k_per_w
    assert math.isfinite(r.nu_eff) and r.nu_eff > 0
