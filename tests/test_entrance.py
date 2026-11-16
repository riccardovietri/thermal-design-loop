"""Developing-flow and inlet/outlet option of the 1D screen against the sources' own tabulated and worked values.

Sources (full text read, see docs/literature.md):
  MY-I   Muzychka & Yovanovich (2002), Laminar flow friction and heat transfer in non-circular ducts and channels,
         Part I: hydrodynamic problem. Compact Heat Exchangers (Festschrift for R.K. Shah), Grenoble, pp. 123-130.
  MY-II  Same authors, Part II: thermal problem, pp. 131-139.
  RA     dos Reis & Altemani (1999), Pressure drop and heat transfer in airflow cooled finned heat sinks using
         correlations, COBEM 1999 (fits to the Kays 1950 laminar loss coefficients).
"""
import math
from dataclasses import fields

import numpy as np
import pytest

from tools.screen_1d import (analyze_1d, effective_nu_factor, f_app_re_sqrt_a, f_pr_uwf, f_re_sqrt_a, hydrodynamic_entrance_length,
                             k_contraction, k_expansion, nu_combined_sqrt_a, nu_fd_sqrt_a, nu_graetz_sqrt_a, screen_1d,
                             thermal_entrance_length)
from tools.types import Design, DesignSpace, Fluid, Requirements, Screen1DEntranceResult, Screen1DResult, to_dict

REQ = Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, 40.0, 5e3, "PG25")
FLUID = Fluid("PG25-like", 40.0, 1010.6, 3958.0, 0.481, 1.394e-3)  # SPEC values, Pr = 11.5
SPACE = DesignSpace((0.3e-3, 0.42e-3, 1.0e-3), (0.26e-3, 0.44e-3), (2.5e-3, 4e-3), (1.5,), 0.15e-3, 10.0)

# MY-I Table 1: aspect ratio -> (f Re_sqrtA exact (Shah & London), f Re_sqrtA from Eq. 32)
MY_I_TABLE_1 = {0.01: (119.56, 119.56), 0.05: (52.77, 52.77), 0.1: (36.82, 36.81), 0.2: (25.59, 25.57), 0.3: (20.78, 20.76),
                0.4: (18.12, 18.09), 0.5: (16.49, 16.46), 0.6: (15.47, 15.43), 0.7: (14.84, 14.79), 0.8: (14.47, 14.40),
                0.9: (14.28, 14.20), 1.0: (14.23, 14.13)}


@pytest.mark.parametrize("eps", MY_I_TABLE_1)
def test_fre_sqrt_a_reproduces_my_part_i_table_1(eps):
    exact, eq32 = MY_I_TABLE_1[eps]
    assert f_re_sqrt_a(eps) == pytest.approx(eq32, abs=0.006)  # table printed to two decimals
    assert abs(f_re_sqrt_a(eps) / exact - 1) < 0.0075  # single-term error, at most 0.7% (MY-I text)


def test_hydrodynamic_entrance_length_matches_my_part_i_worked_values():
    # MY-I after Eq. 44: "L+_h = 0.059 when eps = 1 and L+_h = 0.00083 when eps = 0.01"
    assert hydrodynamic_entrance_length(1.0) == pytest.approx(0.059, abs=0.0005)
    assert hydrodynamic_entrance_length(0.01) == pytest.approx(0.00083, abs=0.000005)


def test_apparent_friction_asymptotes():
    eps = 0.1
    assert f_app_re_sqrt_a(eps, 1e-8) * math.sqrt(1e-8) == pytest.approx(3.44, rel=1e-4)  # short duct, MY-I Eq. 41
    assert f_app_re_sqrt_a(eps, 1e6) == pytest.approx(f_re_sqrt_a(eps), rel=1e-6)  # long duct


def test_thermal_entrance_length_matches_my_part_ii_worked_values():
    # MY-II after Eq. 53: circular duct (eps = 1) "L*_t = 0.0358 and L*_t = 0.0366 ... for (T) and (H)"
    assert thermal_entrance_length(1.0, 3.66, 1.0, 0.501) == pytest.approx(0.0366, abs=0.00005)
    assert thermal_entrance_length(1.0, 3.01, 1.0, 0.409) == pytest.approx(0.0358, abs=0.00005)


def test_flat_plate_function_reproduces_its_asymptotes():
    # MY-II Eq. 43, uniform flux: Nu_z / Re_z^(1/2) -> 0.886 Pr^(1/2) (Pr -> 0) and 0.464 Pr^(1/3) (Pr -> infinity);
    # Eq. 48 is written as Nu / (Re Pr)^(1/2) = f(Pr)
    assert f_pr_uwf(1e-12) == pytest.approx(0.886, rel=1e-3)
    pr = 1e12
    assert f_pr_uwf(pr) * pr ** 0.5 / pr ** (1 / 3) == pytest.approx(0.464, rel=1e-3)


def test_leveque_constant_of_my_part_ii():
    # MY-II Eq. 39: local uniform-flux Leveque coefficient Gamma(2/3) / 18^(1/3) = 0.5167 (the model uses 0.501)
    assert math.gamma(2 / 3) / 18 ** (1 / 3) == pytest.approx(0.5167, abs=5e-5)


def test_nusselt_models_reach_the_developed_limit():
    eps, pr = 0.1, 11.5
    nu_fd = nu_fd_sqrt_a(eps)
    assert nu_fd == pytest.approx(3.66 * f_re_sqrt_a(eps) / (8 * math.sqrt(math.pi) * eps**0.1))
    assert nu_graetz_sqrt_a(eps, 10.0) == pytest.approx(nu_fd, rel=1e-6)
    assert nu_combined_sqrt_a(eps, 10.0, pr) == pytest.approx(nu_fd, rel=1e-4)
    # near the inlet, the combined (developing-velocity) value is above the Graetz value
    assert nu_combined_sqrt_a(eps, 1e-4, pr) > nu_graetz_sqrt_a(eps, 1e-4) > nu_fd


def test_effective_nu_factor_is_a_harmonic_mean_and_vanishes_for_long_channels():
    eps, pr = 0.1, 11.5
    l_star = 0.05
    z = l_star * np.linspace(0, 1, 200001)[1:]
    harmonic = 1 / np.mean(1 / np.array([nu_graetz_sqrt_a(eps, zi) for zi in z[::50]]))
    assert effective_nu_factor(eps, l_star) == pytest.approx(harmonic / nu_fd_sqrt_a(eps), rel=2e-3)
    factors = [effective_nu_factor(eps, l, pr) for l in (0.01, 0.1, 1.0, 10.0, 1000.0)]
    assert all(a > b for a, b in zip(factors, factors[1:]))
    assert factors[-1] == pytest.approx(1.0, abs=1e-3)


def test_loss_coefficients_reproduce_the_source_fits():
    # RA Eq. 5 with Table 1 (laminar): flat duct at alpha = 0, square duct at alpha = 1, linear in between
    s = 0.5
    assert k_contraction(s, 0.0) == pytest.approx(0.800 + 0.029 * s - 0.430 * s**2)
    assert k_contraction(s, 1.0) == pytest.approx(1.190 - 0.011 * s - 0.389 * s**2)
    assert k_expansion(s, 0.0) == pytest.approx(1 - 2.400 * s + s**2)
    assert k_expansion(s, 1.0) == pytest.approx(1 - 2.800 * s + s**2)
    assert k_contraction(s, 0.25) == pytest.approx(0.75 * k_contraction(s, 0) + 0.25 * k_contraction(s, 1))


def test_laminar_expansion_fit_is_the_momentum_balance_with_a_parabolic_profile():
    # Momentum balance across an abrupt expansion gives Ke = 1 - 2 beta sigma + sigma^2, beta = mean(u^2) / U^2.
    # For plane Poiseuille flow beta = 6/5, so the flat-duct laminar fit (1 - 2.4 sigma + sigma^2) assumes a developed
    # profile at the core exit; at sigma = 0 every Ke is 1 (all dynamic pressure lost).
    y = np.linspace(-1, 1, 200001)
    u = 1 - y**2
    beta = np.trapezoid(u**2, y) / 2 / (np.trapezoid(u, y) / 2) ** 2
    assert beta == pytest.approx(1.2, rel=1e-6)
    for s in (0.2, 0.5, 0.8):
        assert k_expansion(s, 0.0) == pytest.approx(1 - 2 * beta * s + s**2, rel=1e-6)
    assert k_expansion(0.0, 0.3) == pytest.approx(1.0)


def test_default_path_is_unchanged():
    d = Design(0.42e-3, 0.44e-3, 4e-3, 1.5)
    default = analyze_1d(d, REQ, FLUID, SPACE)
    explicit = analyze_1d(d, REQ, FLUID, SPACE, entrance_effects=False)
    assert type(default) is Screen1DResult and type(explicit) is Screen1DResult
    assert to_dict(default) == to_dict(explicit)
    assert [to_dict(r) for r in screen_1d(REQ, FLUID, SPACE)] == [to_dict(r) for r in screen_1d(REQ, FLUID, SPACE, False)]
    on = analyze_1d(d, REQ, FLUID, SPACE, entrance_effects=True)
    assert isinstance(on, Screen1DEntranceResult)
    # the fully developed values carried in the entrance result are exactly the default ones
    assert on.dp_fully_developed_pa == default.dp_pa
    assert on.r_th_fully_developed_k_per_w == default.r_th_k_per_w
    assert on.nusselt_fully_developed == default.nusselt
    assert {f.name for f in fields(Screen1DResult)} <= set(to_dict(on))


def test_entrance_option_directions_and_bookkeeping():
    d = Design(0.42e-3, 0.44e-3, 4e-3, 1.5)
    r = analyze_1d(d, REQ, FLUID, entrance_effects=True)
    assert r.dp_pa > r.dp_fully_developed_pa and r.f_app_factor > 1
    assert r.r_th_k_per_w < r.r_th_fully_developed_k_per_w and r.nu_factor_combined > 1
    assert r.r_th_thermal_entrance_k_per_w < r.r_th_fully_developed_k_per_w
    assert r.dp_with_inlet_outlet_pa == pytest.approx(r.dp_pa + r.dp_inlet_outlet_pa)
    assert r.dp_hydrodynamic_entrance_pa == pytest.approx(r.dp_pa - r.dp_fully_developed_pa)
    assert r.r_th_k_per_w == pytest.approx(r.r_base + r.r_conv + r.r_caloric)
    dh = 2 * 0.42e-3 * 4e-3 / (0.42e-3 + 4e-3)
    assert r.x_star == pytest.approx(0.060 / (dh * r.reynolds * FLUID.cp * FLUID.mu / FLUID.k))


def test_entrance_effects_vanish_for_a_long_plate():
    d = Design(0.42e-3, 0.44e-3, 4e-3, 1.5)
    long_req = Requirements(1000.0, 0.060, 60.0, 0.002, 390.0, 40.0, 5e6, "PG25")  # 60 m flow length
    r = analyze_1d(d, long_req, FLUID, entrance_effects=True)
    assert r.f_app_factor == pytest.approx(1.0, abs=1e-4)
    assert r.nu_factor_combined == pytest.approx(1.0, abs=1e-3)
    assert r.r_th_k_per_w == pytest.approx(r.r_th_fully_developed_k_per_w, rel=1e-3)
