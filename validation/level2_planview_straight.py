"""Level 2 plan-view validation on straight fins (gates G1-G3 of the plan-view study).

G1  developed-periodic straight channel: pressure gradient vs plane Poiseuille (Fanning f Re = 24 on Dh = 2w, i.e.
    dp/dx = 12 mu U / w^2) within 0.5%, and Nu on Dh = 2w vs 8.235 (both walls at uniform flux) within 1%.
G2  energy balance within 1%.
G3  Level 1 recommended straight design (examples/output/gpu_coldplate/results.json) through the plan-view network vs the
    Level 1 CFD: dP and R_th differences reported; stop if either exceeds 15%. Expected a few percent low on both,
    because the tall-fin model has no floor or cover drag and puts the side-wall h on the floor.
Informational: a developing straight channel over the real 60 mm length (uniform inflow) to show how far from developed
the straight-fin thermal field is over the plate (entrance effects are excluded from the study by design).

    uv run python validation/level2_planview_straight.py
"""
import json
import sys
import time
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from level2.planforms import straight  # noqa: E402
from level2.planview import run_design, run_developing, superficial_velocity  # noqa: E402
from tools.fluids import fluid_props  # noqa: E402
from tools.types import Requirements  # noqa: E402

REQ = Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, 40.0, 6e3, "PG25", r_th_target_k_per_w=15e-3)
FLOW_LPM, CELLS = 1.5, 10.0
G1_DP, G1_NU, G2_ENERGY, G3_LIMIT = 0.005, 0.01, 0.01, 0.15


def main() -> int:
    t0 = time.perf_counter()
    fluid = fluid_props(REQ.coolant, REQ.inlet_temp_c)
    level1 = json.loads(Path("examples/output/gpu_coldplate/results.json").read_text())
    rec = next(d for d in level1["designs"] if d["cfd"]["design_id"] == level1["recommended"])
    geo = rec["one_d"]["design"]
    w, t, h_fin = geo["w_ch_m"] * 1e3, geo["t_fin_m"] * 1e3, geo["h_fin_m"] * 1e3
    pf = straight(w, t, h_fin)
    r = run_design(pf, REQ, fluid, FLOW_LPM, cells=CELLS)

    u_ch = superficial_velocity(REQ, FLOW_LPM, h_fin) / pf.porosity  # mean channel velocity
    dpdx_exact = 12 * fluid.mu * u_ch / (w * 1e-3) ** 2
    nu_2w = r.h_w_m2k * 2 * w * 1e-3 / fluid.k
    g1_dp, g1_nu = r.dpdx_pa_per_m / dpdx_exact - 1, nu_2w / 8.235 - 1
    g3_dp = r.dp_pa / rec["cfd"]["dp_pa"] - 1
    g3_r = r.r_th_k_per_w / rec["cfd"]["r_th_k_per_w"] - 1

    print(f"Level 1 recommended design: w {w:.3f} mm, t {t:.3f} mm, H {h_fin:.1f} mm, {FLOW_LPM} LPM PG25")
    print(f"Re on Dh = 2w: {r.reynolds_dh:.1f} (mean channel velocity {u_ch:.4f} m/s); plan-view cell {r.n_unknowns} unknowns")
    print(f"G1 pressure gradient {r.dpdx_pa_per_m:.1f} Pa/m vs 12 mu U / w^2 = {dpdx_exact:.1f}: {100 * g1_dp:+.4f}% "
          f"(f Re = {24 * (1 + g1_dp):.4f})")
    print(f"G1 Nu on 2w {nu_2w:.4f} vs 8.235: {100 * g1_nu:+.4f}%")
    print(f"G2 energy balance {100 * r.energy_balance_error:+.2e}%")
    print(f"G3 dP   plan view {r.dp_pa:8.1f} Pa   vs Level 1 CFD {rec['cfd']['dp_pa']:8.1f} Pa   {100 * g3_dp:+.2f}%")
    print(f"G3 R_th plan view {1e3 * r.r_th_k_per_w:8.3f} mK/W vs Level 1 CFD {1e3 * rec['cfd']['r_th_k_per_w']:8.3f} mK/W "
          f"{100 * g3_r:+.2f}%  (base {1e3 * r.r_base:.3f}, convection {1e3 * r.r_conv:.3f}, caloric {1e3 * r.r_caloric:.3f}; "
          f"fin efficiency {r.fin_efficiency:.3f})")

    pf_long = replace(pf, period_mm=5.0)  # 5 mm reporting windows over the 60 mm plate
    dev = run_developing(pf_long, REQ, fluid, n_periods=12, cells=CELLS / 2, lead_in_mm=0.0)
    hs = [p["h_w_m2k"] for p in dev["periods"]]
    print("Informational, developing straight channel from a uniform inlet over 60 mm (5 mm windows): "
          "h / developed h = " + " ".join(f"{h / r.h_w_m2k:.2f}" for h in hs))
    print(f"  length-mean h is {100 * (sum(hs) / len(hs) / r.h_w_m2k - 1):+.1f}% above the developed value; "
          f"energy balance {100 * dev['energy_balance_error']:+.3f}%")

    ok = abs(g1_dp) < G1_DP and abs(g1_nu) < G1_NU and abs(r.energy_balance_error) < G2_ENERGY and \
        abs(dev["energy_balance_error"]) < G2_ENERGY and abs(g3_dp) < G3_LIMIT and abs(g3_r) < G3_LIMIT
    Path(__file__).with_name("level2_planview_straight.json").write_text(json.dumps({
        "design_mm": {"w": w, "t": t, "h_fin": h_fin}, "reynolds_dh": r.reynolds_dh,
        "g1_dpdx_error": g1_dp, "g1_nu_error": g1_nu, "g2_energy_balance_error": r.energy_balance_error,
        "g3_dp_pa": {"plan_view": r.dp_pa, "level1_cfd": rec["cfd"]["dp_pa"], "difference": g3_dp},
        "g3_r_th_k_per_w": {"plan_view": r.r_th_k_per_w, "level1_cfd": rec["cfd"]["r_th_k_per_w"], "difference": g3_r},
        "developing_60mm_h_ratio": [h / r.h_w_m2k for h in hs], "developing_energy_balance_error": dev["energy_balance_error"],
        "seconds": time.perf_counter() - t0}, indent=2))
    print(f"{'PASS' if ok else 'FAIL'} (G1: 0.5% on dp/dx, 1% on Nu; G2: 1%; G3: 15%) in {time.perf_counter() - t0:.0f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
