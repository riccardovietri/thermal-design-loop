"""report: Pareto plot, 1D-vs-CFD parity plot and report.md for a completed study."""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .types import CfdResult, Screen1DResult

ACCENT, MUTED = "#c2410c", "#9ca3af"


def parity_stats(one_d: list[float], cfd: list[float]) -> dict:
    diff = 100 * (np.asarray(one_d) / np.asarray(cfd) - 1)
    return {"mean_pct": float(diff.mean()), "mean_abs_pct": float(np.abs(diff).mean()), "max_abs_pct": float(np.abs(diff).max())}


def plot_pareto(sweep: list[Screen1DResult], front: list[Screen1DResult], cfd: list[CfdResult], recommended_id: str,
                dp_budget_pa: float, out: Path) -> Path:
    fig, ax = plt.subplots(figsize=(7.5, 5.2), constrained_layout=True)
    feas = [s for s in sweep if s.feasible]
    ax.scatter([s.dp_pa / 1e3 for s in feas], [s.r_th_k_per_w * 1e3 for s in feas], s=8, color=MUTED, alpha=0.5,
               label=f"1D screen, feasible ({len(feas)} of {len(sweep)})", linewidths=0)
    ax.plot([s.dp_pa / 1e3 for s in front], [s.r_th_k_per_w * 1e3 for s in front], color="#374151", lw=1.2, label="1D Pareto front")
    ax.scatter([c.dp_pa / 1e3 for c in cfd], [c.r_th_k_per_w * 1e3 for c in cfd], s=46, color="#1d4ed8", zorder=3,
               label=f"2D conjugate CFD ({len(cfd)} designs)")
    rec = next(c for c in cfd if c.design_id == recommended_id)
    ax.scatter([rec.dp_pa / 1e3], [rec.r_th_k_per_w * 1e3], s=190, facecolors="none", edgecolors=ACCENT, linewidths=2.2, zorder=4,
               label="Recommended")
    ax.axvline(dp_budget_pa / 1e3, color=ACCENT, ls="--", lw=1)
    ax.text(dp_budget_pa / 1e3 * 0.95, ax.get_ylim()[1] * 0.97, "ΔP budget", color=ACCENT, ha="right", va="top", fontsize=9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Channel pressure drop ΔP (kPa)")
    ax.set_ylabel("Thermal resistance R_th (mK/W)")
    ax.set_title("Thermal resistance vs pressure drop")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=9)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def plot_parity(pairs: list[tuple[Screen1DResult, CfdResult]], out: Path) -> Path:
    x = np.array([c.r_th_k_per_w * 1e3 for _, c in pairs])
    y = np.array([s.r_th_k_per_w * 1e3 for s, _ in pairs])
    eta = np.array([s.fin_efficiency for s, _ in pairs])
    stats = parity_stats(y, x)
    lo, hi = min(x.min(), y.min()) * 0.9, max(x.max(), y.max()) * 1.1
    fig, ax = plt.subplots(figsize=(6.2, 5.6), constrained_layout=True)
    line = np.array([lo, hi])
    ax.fill_between(line, line * 0.95, line * 1.05, color=MUTED, alpha=0.25, label="±5%")
    ax.plot(line, line, color="#374151", lw=1)
    sc = ax.scatter(x, y, c=eta, cmap="viridis", vmin=min(0.5, eta.min()), vmax=1.0, s=60, zorder=3, edgecolors="white", linewidths=0.6)
    fig.colorbar(sc, ax=ax, label="Fin efficiency (1D)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("R_th, 2D conjugate CFD (mK/W)")
    ax.set_ylabel("R_th, 1D screen (mK/W)")
    ax.set_title(f"1D vs CFD: mean |difference| {stats['mean_abs_pct']:.1f}%, max {stats['max_abs_pct']:.1f}%")
    ax.legend(frameon=False, loc="upper left")
    ax.grid(True, which="both", alpha=0.25)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def write_report(study: dict, out: Path) -> Path:
    req, fluid, pairs = study["requirements"], study["fluid"], study["pairs"]
    rec_s, rec_c = study["recommended"]
    best_s, best_c = study["lowest_r_th"]
    stats, conf = study["parity"], study["confidence"]
    d = rec_s.design
    lines = [
        f"# Cold-plate design study: {study['name']}",
        "",
        f"**Recommended design:** {d.w_ch_m*1e3:.2f} mm channels, {d.t_fin_m*1e3:.2f} mm fins, {d.h_fin_m*1e3:.1f} mm tall "
        f"({rec_c.n_channels} channels) at {d.flow_lpm:.1f} LPM. CFD gives R_th = {rec_c.r_th_k_per_w*1e3:.1f} mK/W at "
        f"ΔP = {rec_c.dp_pa/1e3:.1f} kPa, so the chip interface averages {rec_c.base_avg_temp_c:.1f} °C "
        f"(peak {rec_c.base_max_temp_c:.1f} °C) for {req.heat_load_w:.0f} W and {req.inlet_temp_c:.0f} °C inlet.",
        "",
        "## Why this design",
        f"- It is the knee of the CFD-verified front: the point farthest from the straight line between the two ends of the front in log-log space.",
        f"- The lowest-resistance design within budget ({best_s.design.w_ch_m*1e3:.2f} mm channels, {best_s.design.t_fin_m*1e3:.2f} mm fins, "
        f"{best_s.design.h_fin_m*1e3:.1f} mm tall, {best_s.design.flow_lpm:.1f} LPM) reaches {best_c.r_th_k_per_w*1e3:.1f} mK/W at "
        f"{best_c.dp_pa/1e3:.1f} kPa: {100*(1-best_c.r_th_k_per_w/rec_c.r_th_k_per_w):.0f}% lower resistance for "
        f"{best_c.dp_pa/rec_c.dp_pa:.1f}× the pressure drop.",
        f"- At the recommended point the resistance splits into base conduction {rec_s.r_base*1e3:.1f}, fin and convection {rec_s.r_conv*1e3:.1f}, "
        f"and coolant heating {rec_s.r_caloric*1e3:.1f} mK/W (1D breakdown).",
        "",
        "![Pareto](pareto.png)",
        "",
        "## 1D screen vs CFD",
        f"The 1D screen evaluated {study['n_sweep']} designs ({study['n_feasible']} feasible) in {study['t_screen_s']*1e3:.0f} ms. "
        f"{len(pairs)} designs spread along its Pareto front were solved with the 2D conjugate model in {study['t_cfd_s']:.1f} s.",
        "",
        f"- R_th: 1D differs from CFD by {stats['mean_abs_pct']:.1f}% on average, {stats['max_abs_pct']:.1f}% at worst (mean signed {stats['mean_pct']:+.1f}%).",
        f"- ΔP: identical to within {study['dp_max_abs_pct']:.2f}%. This is expected, not a result: both models assume fully developed flow and share the same f·Re.",
        "",
        "![Parity](parity.png)",
        "",
        "| Channel (mm) | Fin (mm) | Height (mm) | Flow (LPM) | Re | ΔP CFD (kPa) | R_th 1D (mK/W) | R_th CFD (mK/W) | 1D − CFD | Fin eff. |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for s, c in pairs:
        g = s.design
        mark = " **(recommended)**" if c.design_id == rec_c.design_id else ""
        lines.append(f"| {g.w_ch_m*1e3:.2f}{mark} | {g.t_fin_m*1e3:.2f} | {g.h_fin_m*1e3:.1f} | {g.flow_lpm:.1f} | {c.reynolds:.0f} | {c.dp_pa/1e3:.2f} | "
                     f"{s.r_th_k_per_w*1e3:.2f} | {c.r_th_k_per_w*1e3:.2f} | {100*(s.r_th_k_per_w/c.r_th_k_per_w-1):+.1f}% | {s.fin_efficiency:.2f} |")
    lines += ["", "## Confidence checks", "", "| Check | Value | Limit | Result |", "|---|---|---|---|"]
    lines += [f"| {c['check']} | {c['value']} | {c['limit']} | {'pass' if c['pass'] else 'FAIL'} |" for c in conf["checks"]]
    ms = conf["mesh_study"]
    lines += ["", f"Mesh study on `{ms['design_id']}`:", "", "| Cells across half channel | Unknowns | R_th (mK/W) | ΔP (kPa) |", "|---|---|---|---|"]
    lines += [f"| {r['cells_across_half_channel']} | {r['n_dofs']} | {r['r_th_k_per_w']*1e3:.4f} | {r['dp_pa']/1e3:.4f} |" for r in ms["levels"]]
    lines += [
        "", "![Fields](fields_recommended.png)", "",
        "## Inputs",
        f"- Heat load {req.heat_load_w:.0f} W, uniform over {req.footprint_width_m*1e3:.0f} × {req.flow_length_m*1e3:.0f} mm; "
        f"{req.base_thickness_m*1e3:.1f} mm base, k = {req.solid_k:.0f} W/m·K; ΔP budget {req.dp_budget_pa/1e3:.0f} kPa (channel core).",
        f"- {fluid.name} at {fluid.temp_c:.0f} °C (CoolProp): ρ = {fluid.rho:.1f} kg/m³, cp = {fluid.cp:.0f} J/kg·K, k = {fluid.k:.3f} W/m·K, μ = {fluid.mu*1e3:.3f} mPa·s.",
        "- R_th = (average chip-interface temperature − inlet temperature) / heat load. TIM and junction-to-case are excluded.",
        "",
        "## Design basis",
        "- Feature limits (≥ 0.15 mm, aspect ratio ≤ 10) sit inside published skived-copper practice: 0.1–0.8 mm fins with 0.1 mm minimum gap "
        "(Boyd fabrication guide); 100–200 µm fins and gaps, 1–2 mm tall, in direct-to-chip cold plates (ASME mini-channel study).",
        "- Base thickness is a fixed assumption. Published designs span 0.1 mm (manifold microcoolers) to 4.5 mm (a commercial 55 mm plate); no typical value was found.",
        "",
        "## What this model does not capture",
        "- Entrance effects and inlet/outlet losses (both models are fully developed).",
        "- Flow maldistribution between channels, manifolds and fittings.",
        "- Temperature-dependent properties, non-uniform heat flux and heat spreading from a die smaller than the plate.",
        "- Anything three-dimensional. Each CFD run is a 2D cross-section of one fin and channel.",
        "",
        f"Every CFD number traces to `runs/{study['name']}/<design_id>/` (inputs, mesh, results, fields).",
        "",
    ]
    out.write_text("\n".join(lines))
    return out
