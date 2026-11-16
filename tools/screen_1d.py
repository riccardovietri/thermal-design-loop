"""screen_1d: straight-fin cold plate, fully developed laminar flow, 1D resistance network.

f*Re and Nu_H1 (four heated walls) from the Shah & London (1978) rectangular-duct polynomials,
adiabatic-tip fin efficiency on the fin sides, eta = 1 on the channel floor, adiabatic cover.
R_th = R_base + R_conv + Q/(2 m_dot cp), referenced to inlet temperature and average base temperature.

Optional, off by default (`entrance_effects=True`): developing flow and inlet/outlet losses, from
  Muzychka & Yovanovich (2002), "Laminar flow friction and heat transfer in non-circular ducts and channels",
    Part I (hydrodynamic): f_app Re_sqrtA = [ (f Re_sqrtA)^2 + (3.44 / sqrt(L+))^2 ]^(1/2),  L+ = L / (sqrt(A) Re_sqrtA),
      f Re_sqrtA = 12 / (sqrt(eps) (1 + eps) [1 - 192 eps / pi^5 tanh(pi / (2 eps))]),  eps = aspect ratio <= 1.
    Part II (thermal, uniform wall flux): local Nu_sqrtA(z*) for the Graetz problem (Eq. 50) and the combined entry
      problem (Eq. 51), z* = z / (sqrt(A) Re_sqrtA Pr), with C1 = 3.66, C2 = 1 (local), C3 = 0.501, C4 = 1 (local),
      gamma = 1/10, f(Pr) = 0.886 / [1 + (1.909 Pr^(1/6))^(9/2)]^(2/9), m = 2.27 + 1.65 Pr^(1/3) (their Table 2, Eq. 52).
  Both are used as ratios to their own fully developed limits, applied to the Shah & London f Re and Nu_H1 above, so the
  correction vanishes as L grows. Under uniform flux R_th depends on the length-mean of (wall - bulk) temperature, so the
  effective h is the harmonic mean of the local h over the length.
  Inlet/outlet: Kays (1950) laminar abrupt-contraction and abrupt-expansion coefficients for flat and square ducts, as
  polynomial fits in sigma interpolated linearly in aspect ratio by dos Reis & Altemani (COBEM 1999, Eq. 5 and Table 1);
  dP_io = (Kc + Ke) rho U^2 / 2 with U the channel velocity (incompressible; the 1 - sigma^2 terms cancel).
"""
import argparse
import itertools
import json
import math

from scipy.integrate import quad

from .fluids import fluid_props
from .types import Design, DesignSpace, Fluid, Requirements, Screen1DEntranceResult, Screen1DResult, to_dict

RE_LAMINAR_LIMIT = 2000.0


def f_re_fanning(alpha: float) -> float:
    a = alpha
    return 24 * (1 - 1.3553 * a + 1.9467 * a**2 - 1.7012 * a**3 + 0.9564 * a**4 - 0.2537 * a**5)


def nu_h1(alpha: float) -> float:
    a = alpha
    return 8.235 * (1 - 2.0421 * a + 3.0853 * a**2 - 2.4765 * a**3 + 1.0578 * a**4 - 0.1861 * a**5)


def analyze_1d(design: Design, req: Requirements, fluid: Fluid, space: DesignSpace | None = None,
               entrance_effects: bool = False) -> Screen1DResult:
    w, t, h = design.w_ch_m, design.t_fin_m, design.h_fin_m
    n, length = design.n_channels(req.footprint_width_m), req.flow_length_m
    alpha = min(w, h) / max(w, h)
    dh = 2 * w * h / (w + h)
    vdot = design.flow_lpm / 60e3
    v_mean = vdot / (n * w * h)
    reynolds = fluid.rho * v_mean * dh / fluid.mu
    fre = f_re_fanning(alpha)
    dp = 2 * fre * fluid.mu * v_mean * length / dh**2
    nu = nu_h1(alpha)
    htc = nu * fluid.k / dh
    m = math.sqrt(2 * htc / (req.solid_k * t))
    eta = math.tanh(m * h) / (m * h)
    r_conv = 1 / (n * length * htc * (w + 2 * eta * h))
    r_base = req.base_thickness_m / (req.solid_k * n * (w + t) * length)
    r_cal = 1 / (2 * fluid.rho * fluid.cp * vdot)

    violations = []
    if dp > req.dp_budget_pa:
        violations.append("dp_over_budget")
    if reynolds > RE_LAMINAR_LIMIT:
        violations.append("not_laminar")
    if space is not None:
        if min(w, t) < space.min_feature_m - 1e-12:
            violations.append("below_min_feature")
        if h / t > space.max_aspect_ratio + 1e-9:
            violations.append("fin_aspect_ratio")
        if h / w > space.max_aspect_ratio + 1e-9:
            violations.append("channel_aspect_ratio")
    result = Screen1DResult(
        design=design, n_channels=n, reynolds=reynolds, dp_pa=dp, r_th_k_per_w=r_base + r_conv + r_cal,
        r_th_max_k_per_w=r_base + r_conv + 2 * r_cal, r_base=r_base, r_conv=r_conv, r_caloric=r_cal,
        fin_efficiency=eta, nusselt=nu, f_re_fanning=fre, feasible=not violations, violations=tuple(violations),
    )
    if not entrance_effects:
        return result
    return _with_entrance_effects(result, req, fluid)


def screen_1d(req: Requirements, fluid: Fluid, space: DesignSpace, entrance_effects: bool = False) -> list[Screen1DResult]:
    """Full-factorial sweep of the design space."""
    grid = itertools.product(space.w_ch_m, space.t_fin_m, space.h_fin_m, space.flow_lpm)
    return [analyze_1d(Design(*g), req, fluid, space, entrance_effects) for g in grid]


# ---------------------------------------------------------------------------------------------------------------------
# Developing flow and inlet/outlet losses (optional). Muzychka & Yovanovich (2002) Parts I and II; Kays (1950) via
# dos Reis & Altemani (1999). Lengths on sqrt(A); eps = channel aspect ratio <= 1.

MY_SHORT_DUCT = 3.44  # Part I, Eq. (41): mean (apparent) f Re_L = 3.44 / sqrt(L+) near the inlet of any duct
MY_UWF_C1, MY_UWF_C3, MY_GAMMA = 3.66, 0.501, 0.1  # Part II, Table 2: uniform wall flux; gamma upper bound (right angles)
_KC_LAMINAR = ((0.800, 0.029, -0.430), (1.190, -0.011, -0.389))  # (flat duct A0..A2, square duct B0..B2), Table 1
_KE_LAMINAR = ((1.000, -2.400, 1.000), (1.000, -2.800, 1.000))


def f_re_sqrt_a(eps: float) -> float:
    """Fully developed f Re on sqrt(A), first term of the rectangular-duct series (Part I, Eq. 32)."""
    return 12 / (math.sqrt(eps) * (1 + eps) * (1 - 192 * eps / math.pi**5 * math.tanh(math.pi / (2 * eps))))


def f_app_re_sqrt_a(eps: float, l_plus: float) -> float:
    """Apparent (length-mean, includes the momentum change) f Re on sqrt(A) from a uniform inlet (Part I, Eq. 43)."""
    return math.hypot(f_re_sqrt_a(eps), MY_SHORT_DUCT / math.sqrt(l_plus))


def hydrodynamic_entrance_length(eps: float) -> float:
    """L+_h on sqrt(A), where the short- and long-duct asymptotes meet (Part I, Eq. 44)."""
    return (MY_SHORT_DUCT / f_re_sqrt_a(eps)) ** 2


def nu_fd_sqrt_a(eps: float, c1: float = MY_UWF_C1, gamma: float = MY_GAMMA) -> float:
    """Fully developed Nu on sqrt(A) (Part II, Eq. 31)."""
    return c1 * f_re_sqrt_a(eps) / (8 * math.sqrt(math.pi) * eps**gamma)


def f_pr_uwf(pr: float) -> float:
    """Simultaneously developing flat-plate limit, uniform flux (Part II, Eq. 48)."""
    return 0.886 / (1 + (1.909 * pr ** (1 / 6)) ** 4.5) ** (2 / 9)


def nu_graetz_sqrt_a(eps: float, z_star: float, c2: float = 1.0, c3: float = MY_UWF_C3) -> float:
    """Thermally developing, hydrodynamically developed Nu on sqrt(A); c2 = 1 local, 3/2 mean (Part II, Eq. 50)."""
    leveque = c2 * c3 * (f_re_sqrt_a(eps) / z_star) ** (1 / 3)
    return (leveque**5 + nu_fd_sqrt_a(eps) ** 5) ** 0.2


def nu_combined_sqrt_a(eps: float, z_star: float, pr: float, c2: float = 1.0, c4: float = 1.0) -> float:
    """Simultaneously developing Nu on sqrt(A); c2, c4 = 1 local, 3/2 and 2 mean (Part II, Eqs. 51-52)."""
    m = 2.27 + 1.65 * pr ** (1 / 3)
    flat_plate = c4 * f_pr_uwf(pr) / math.sqrt(z_star)
    return (nu_graetz_sqrt_a(eps, z_star, c2) ** m + flat_plate**m) ** (1 / m)


def thermal_entrance_length(eps: float, c1: float, c2: float, c3: float, gamma: float = MY_GAMMA) -> float:
    """z*_t on sqrt(A), where the Graetz asymptotes meet (Part II, Eq. 53)."""
    s = 1 - 192 * eps / math.pi**5 * math.tanh(math.pi / (2 * eps))
    return 19.80 * (c2 * c3 / c1) ** 3 * eps ** (1 + 3 * gamma) * (1 + eps) ** 2 * s**2


def effective_nu_factor(eps: float, l_star: float, pr: float | None = None) -> float:
    """Harmonic length-mean of the local Nu over 0 < z* < l_star, divided by the developed Nu.

    pr None: thermal entrance only (Graetz); otherwise simultaneously developing at that Prandtl number.
    The substitution z* = l_star s^3 removes the z*^(1/3) behaviour of 1/Nu at the inlet.
    """
    nu = (lambda z: nu_graetz_sqrt_a(eps, z)) if pr is None else (lambda z: nu_combined_sqrt_a(eps, z, pr))
    mean_inv, _ = quad(lambda s: 3 * s * s / nu(l_star * s**3), 0.0, 1.0, epsabs=0, epsrel=1e-10, limit=200)
    return 1 / (mean_inv * nu_fd_sqrt_a(eps))


def _poly(c, x):
    return c[0] + c[1] * x + c[2] * x * x


def k_contraction(sigma: float, alpha: float) -> float:
    """Laminar abrupt-contraction coefficient: Kays (1950) as fitted by dos Reis & Altemani (1999, Eq. 5)."""
    return (1 - alpha) * _poly(_KC_LAMINAR[0], sigma) + alpha * _poly(_KC_LAMINAR[1], sigma)


def k_expansion(sigma: float, alpha: float) -> float:
    """Laminar abrupt-expansion coefficient: Kays (1950) as fitted by dos Reis & Altemani (1999, Eq. 5)."""
    return (1 - alpha) * _poly(_KE_LAMINAR[0], sigma) + alpha * _poly(_KE_LAMINAR[1], sigma)


def _with_entrance_effects(fd: Screen1DResult, req: Requirements, fluid: Fluid) -> Screen1DEntranceResult:
    d = fd.design
    w, t, h = d.w_ch_m, d.t_fin_m, d.h_fin_m
    n, length = fd.n_channels, req.flow_length_m
    eps = min(w, h) / max(w, h)
    dh = 2 * w * h / (w + h)
    sqrt_a = math.sqrt(w * h)
    v_mean = d.flow_lpm / 60e3 / (n * w * h)
    pr = fluid.cp * fluid.mu / fluid.k
    l_plus = length * fluid.mu / (fluid.rho * v_mean * sqrt_a**2)  # L / (sqrt(A) Re_sqrtA)
    f_factor = f_app_re_sqrt_a(eps, l_plus) / f_re_sqrt_a(eps)
    dp = fd.dp_pa * f_factor
    sigma = w / (w + t)
    kc, ke = k_contraction(sigma, eps), k_expansion(sigma, eps)
    dp_io = (kc + ke) * 0.5 * fluid.rho * v_mean**2

    def r_conv(factor):
        htc = fd.nusselt * factor * fluid.k / dh
        m = math.sqrt(2 * htc / (req.solid_k * t))
        eta = math.tanh(m * h) / (m * h)
        return 1 / (n * length * htc * (w + 2 * eta * h)), eta

    nu_th, nu_comb = effective_nu_factor(eps, l_plus / pr), effective_nu_factor(eps, l_plus / pr, pr)
    r_th_only, _ = r_conv(nu_th)
    r_comb, eta = r_conv(nu_comb)
    violations = [v for v in fd.violations if v != "dp_over_budget"]
    if dp > req.dp_budget_pa:
        violations.insert(0, "dp_over_budget")
    base = {k: getattr(fd, k) for k in fd.__dataclass_fields__}
    base.update(dp_pa=dp, r_th_k_per_w=fd.r_base + r_comb + fd.r_caloric, r_th_max_k_per_w=fd.r_base + r_comb + 2 * fd.r_caloric,
                r_conv=r_comb, fin_efficiency=eta, nusselt=fd.nusselt * nu_comb, feasible=not violations,
                violations=tuple(violations))
    return Screen1DEntranceResult(
        **base, dp_fully_developed_pa=fd.dp_pa, dp_hydrodynamic_entrance_pa=dp - fd.dp_pa, dp_inlet_outlet_pa=dp_io,
        dp_with_inlet_outlet_pa=dp + dp_io, r_th_fully_developed_k_per_w=fd.r_th_k_per_w,
        r_th_thermal_entrance_k_per_w=fd.r_base + r_th_only + fd.r_caloric, nusselt_fully_developed=fd.nusselt,
        nu_factor_thermal=nu_th, nu_factor_combined=nu_comb, f_app_factor=f_factor, k_contraction=kc, k_expansion=ke,
        sigma=sigma, x_star=length / (dh * fd.reynolds * pr), l_plus_dh=length / (dh * fd.reynolds))


def main() -> None:
    ap = argparse.ArgumentParser(description="1D analysis of one straight-fin design (1 kW, 60x60 mm, 2 mm copper base).")
    for name in ("w_ch_mm", "t_fin_mm", "h_fin_mm", "flow_lpm"):
        ap.add_argument(f"--{name}", type=float, required=True)
    ap.add_argument("--coolant", default="PG25")
    ap.add_argument("--inlet_c", type=float, default=40.0)
    ap.add_argument("--entrance", action="store_true", help="add developing flow and inlet/outlet losses (off by default)")
    a = ap.parse_args()
    req = Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, a.inlet_c, 5e3, a.coolant, r_th_threshold_k_per_w=20e-3,
                       r_th_objective_k_per_w=15e-3)
    d = Design(a.w_ch_mm / 1e3, a.t_fin_mm / 1e3, a.h_fin_mm / 1e3, a.flow_lpm)
    print(json.dumps(to_dict(analyze_1d(d, req, fluid_props(a.coolant, a.inlet_c), entrance_effects=a.entrance)), indent=2))


if __name__ == "__main__":
    main()
