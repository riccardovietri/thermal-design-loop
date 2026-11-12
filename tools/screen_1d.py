"""screen_1d: straight-fin cold plate, fully developed laminar flow, 1D resistance network.

f*Re and Nu_H1 (four heated walls) from the Shah & London (1978) rectangular-duct polynomials,
adiabatic-tip fin efficiency on the fin sides, eta = 1 on the channel floor, adiabatic cover.
R_th = R_base + R_conv + Q/(2 m_dot cp), referenced to inlet temperature and average base temperature.
"""
import argparse
import itertools
import json
import math

from .fluids import fluid_props
from .types import Design, DesignSpace, Fluid, Requirements, Screen1DResult, to_dict

RE_LAMINAR_LIMIT = 2000.0


def f_re_fanning(alpha: float) -> float:
    a = alpha
    return 24 * (1 - 1.3553 * a + 1.9467 * a**2 - 1.7012 * a**3 + 0.9564 * a**4 - 0.2537 * a**5)


def nu_h1(alpha: float) -> float:
    a = alpha
    return 8.235 * (1 - 2.0421 * a + 3.0853 * a**2 - 2.4765 * a**3 + 1.0578 * a**4 - 0.1861 * a**5)


def analyze_1d(design: Design, req: Requirements, fluid: Fluid, space: DesignSpace | None = None) -> Screen1DResult:
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
    return Screen1DResult(
        design=design, n_channels=n, reynolds=reynolds, dp_pa=dp, r_th_k_per_w=r_base + r_conv + r_cal,
        r_th_max_k_per_w=r_base + r_conv + 2 * r_cal, r_base=r_base, r_conv=r_conv, r_caloric=r_cal,
        fin_efficiency=eta, nusselt=nu, f_re_fanning=fre, feasible=not violations, violations=tuple(violations),
    )


def screen_1d(req: Requirements, fluid: Fluid, space: DesignSpace) -> list[Screen1DResult]:
    """Full-factorial sweep of the design space."""
    grid = itertools.product(space.w_ch_m, space.t_fin_m, space.h_fin_m, space.flow_lpm)
    return [analyze_1d(Design(*g), req, fluid, space) for g in grid]


def main() -> None:
    ap = argparse.ArgumentParser(description="1D analysis of one straight-fin design (1 kW, 60x60 mm, 2 mm copper base).")
    for name in ("w_ch_mm", "t_fin_mm", "h_fin_mm", "flow_lpm"):
        ap.add_argument(f"--{name}", type=float, required=True)
    ap.add_argument("--coolant", default="PG25")
    ap.add_argument("--inlet_c", type=float, default=40.0)
    a = ap.parse_args()
    req = Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, a.inlet_c, 5e3, a.coolant, r_th_threshold_k_per_w=20e-3,
                       r_th_objective_k_per_w=15e-3)
    d = Design(a.w_ch_mm / 1e3, a.t_fin_mm / 1e3, a.h_fin_mm / 1e3, a.flow_lpm)
    print(json.dumps(to_dict(analyze_1d(d, req, fluid_props(a.coolant, a.inlet_c))), indent=2))


if __name__ == "__main__":
    main()
