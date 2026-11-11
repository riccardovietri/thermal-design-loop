"""report: Pareto plot, 1D-vs-CFD parity plot and report.md for a completed study."""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator

from .types import CfdResult, Screen1DResult

ACCENT, MUTED = "#c2410c", "#9ca3af"


def parity_stats(one_d: list[float], cfd: list[float]) -> dict:
    diff = 100 * (np.asarray(one_d) / np.asarray(cfd) - 1)
    return {"mean_pct": float(diff.mean()), "mean_abs_pct": float(np.abs(diff).mean()), "max_abs_pct": float(np.abs(diff).max())}


RULE_TEXT = {
    "meets_target_lowest_dp": "at least one CFD-verified design meets both the budget and the target, so the one with the lowest ΔP is recommended",
    "fallback_lowest_r_th": "no CFD-verified design meets the target within the budget, so the lowest R_th within budget is recommended instead",
}
KNEE_DEFINITION = ("the point on the CFD-verified front farthest from the straight line joining the front's two end points, "
                   "in log10(ΔP)–log10(R_th) coordinates; it depends on where the front's end points are")
X_MIN_KPA = 0.3
NICE_TICKS = (0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 70, 100, 150, 200)


def plain_log_ticks(ax, lo: float, hi: float, axis: str = "x") -> None:
    """Readable tick labels (0.3, 1, 2, ...) on a log axis instead of matplotlib's 3×10⁻¹ style."""
    a = ax.xaxis if axis == "x" else ax.yaxis
    ticks = [t for t in NICE_TICKS if lo <= t <= hi]
    a.set_major_locator(FixedLocator(ticks))
    a.set_major_formatter(lambda v, _: f"{v:g}")
    a.set_minor_locator(NullLocator())
    a.set_minor_formatter(NullFormatter())


def plot_pareto(sweep: list[Screen1DResult], front: list[Screen1DResult], cfd: list[CfdResult], recommended_id: str,
                dp_budget_pa: float, out: Path, r_th_target: float | None = None, knee_id: str | None = None) -> Path:
    """R_th vs ΔP restricted to the region of interest: 0.3 kPa to 1.5× budget, 0.7× lowest feasible R_th to 2× target."""
    fig, ax = plt.subplots(figsize=(7.5, 5.2), constrained_layout=True)
    feas = [s for s in sweep if s.feasible]
    budget_kpa = dp_budget_pa / 1e3
    in_budget = [c.r_th_k_per_w for c in cfd if c.dp_pa <= dp_budget_pa] + [s.r_th_k_per_w for s in feas]
    y_lo = 0.7 * min(in_budget) * 1e3
    y_hi = 2 * r_th_target * 1e3 if r_th_target is not None else 1.3 * max(c.r_th_k_per_w for c in cfd) * 1e3
    x_lo, x_hi = X_MIN_KPA, 1.5 * budget_kpa
    if r_th_target is not None:
        ax.fill_between([x_lo, budget_kpa], y_lo, r_th_target * 1e3, color="#16a34a", alpha=0.08, lw=0, zorder=0,
                        label="Meets budget and target")
    ax.scatter([s.dp_pa / 1e3 for s in feas], [s.r_th_k_per_w * 1e3 for s in feas], s=8, color=MUTED, alpha=0.5,
               label=f"1D screen, feasible ({len(feas)} of {len(sweep)})", linewidths=0)
    ax.plot([s.dp_pa / 1e3 for s in front], [s.r_th_k_per_w * 1e3 for s in front], color="#374151", lw=1.2, label="1D Pareto front")
    ax.scatter([c.dp_pa / 1e3 for c in cfd], [c.r_th_k_per_w * 1e3 for c in cfd], s=46, color="#1d4ed8", zorder=3,
               label=f"2D conjugate CFD ({len(cfd)} designs)")
    rec = next(c for c in cfd if c.design_id == recommended_id)
    ax.scatter([rec.dp_pa / 1e3], [rec.r_th_k_per_w * 1e3], s=190, facecolors="none", edgecolors=ACCENT, linewidths=2.2, zorder=4,
               label="Recommended")
    if knee_id is not None and knee_id != recommended_id:
        knee = next(c for c in cfd if c.design_id == knee_id)
        ax.scatter([knee.dp_pa / 1e3], [knee.r_th_k_per_w * 1e3], s=150, facecolors="none", edgecolors="#374151", linewidths=1.4,
                   ls="--", zorder=4, label="Knee of CFD front (secondary)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    plain_log_ticks(ax, x_lo, x_hi, "x")
    plain_log_ticks(ax, y_lo, y_hi, "y")
    ax.axvline(budget_kpa, color=ACCENT, ls="--", lw=1)
    ax.text(budget_kpa * 0.96, y_hi * 0.97, f"ΔP budget {budget_kpa:g} kPa", color=ACCENT, ha="right", va="top", fontsize=9)
    if r_th_target is not None:
        ax.axhline(r_th_target * 1e3, color="#15803d", ls="--", lw=1)
        ax.text(x_lo * 1.04, r_th_target * 1e3 * 1.02, f"R_th target {r_th_target*1e3:g} mK/W", color="#15803d", ha="left",
                va="bottom", fontsize=9)
    ax.set_xlabel("Channel pressure drop ΔP (kPa)")
    ax.set_ylabel("Thermal resistance R_th (mK/W)")
    ax.set_title("Thermal resistance vs pressure drop at rated flow (region of interest)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=9, loc="lower left")
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
    knee_s, knee_c = study["knee"]
    rule = study["rule"]
    target = req.r_th_target_k_per_w
    best_s, best_c = study["lowest_r_th"]
    stats, conf = study["parity"], study["confidence"]
    d = rec_s.design
    lines = [
        f"# Cold-plate design study: {study['name']}",
        "",
        f"**Recommended design:** {d.w_ch_m*1e3:.2f} mm channels, {d.t_fin_m*1e3:.2f} mm fins, {d.h_fin_m*1e3:.1f} mm tall "
        f"({rec_c.n_channels} channels) at {d.flow_lpm:.1f} LPM (rated flow). CFD gives R_th = {rec_c.r_th_k_per_w*1e3:.1f} mK/W at "
        f"ΔP = {rec_c.dp_pa/1e3:.1f} kPa, so the chip interface averages {rec_c.base_avg_temp_c:.1f} °C "
        f"(peak {rec_c.base_max_temp_c:.1f} °C) for {req.heat_load_w:.0f} W and {req.inlet_temp_c:.0f} °C inlet.",
        "",
        "## Why this design",
        "- Rule: among the CFD-verified designs that meet both the ΔP budget and the R_th target, recommend the one with the lowest ΔP. "
        "If none meets the target, recommend the lowest R_th within budget and say so.",
        f"- Budget: channel-core ΔP ≤ {req.dp_budget_pa/1e3:g} kPa at rated flow. Target: R_th ≤ "
        + (f"{target*1e3:g} mK/W." if target is not None else "none set."),
        f"- Rule that fired: `{rule}`: {RULE_TEXT[rule]}.",
        f"- Secondary: the knee of the CFD-verified front is {knee_s.design.w_ch_m*1e3:.2f} mm channels, {knee_s.design.t_fin_m*1e3:.2f} mm fins, "
        f"{knee_s.design.h_fin_m*1e3:.1f} mm tall: {knee_c.r_th_k_per_w*1e3:.1f} mK/W at {knee_c.dp_pa/1e3:.2f} kPa"
        + (" (the same design)." if knee_c.design_id == rec_c.design_id else ".")
        + f" The knee is {KNEE_DEFINITION}.",
        f"- The lowest-resistance design within budget ({best_s.design.w_ch_m*1e3:.2f} mm channels, {best_s.design.t_fin_m*1e3:.2f} mm fins, "
        f"{best_s.design.h_fin_m*1e3:.1f} mm tall) reaches {best_c.r_th_k_per_w*1e3:.1f} mK/W at {best_c.dp_pa/1e3:.2f} kPa"
        + (" (the recommended design)." if best_c.design_id == rec_c.design_id else
           f": {100*(1-best_c.r_th_k_per_w/rec_c.r_th_k_per_w):.0f}% lower resistance for {best_c.dp_pa/rec_c.dp_pa:.1f}× the pressure drop."),
        f"- At the recommended point the resistance splits into base conduction {rec_s.r_base*1e3:.1f}, fin and convection {rec_s.r_conv*1e3:.1f}, "
        f"and coolant heating {rec_s.r_caloric*1e3:.1f} mK/W (1D breakdown).",
        "",
        "![Pareto](pareto.png)",
        "",
        "The plot shows the region of interest only: ΔP from 0.3 kPa to 1.5× the budget, R_th from 0.7× the lowest feasible value "
        "to 2× the target. Designs outside it are screened but not drawn.",
        "",
        "## Flow sensitivity of the recommended design",
        "",
        "| Flow (LPM) | ΔP (kPa) | R_th (mK/W) | Average interface (°C) | Peak interface (°C) |",
        "|---|---|---|---|---|",
    ]
    for r in study["flow_sensitivity"]:
        mark = " (rated)" if abs(r["flow_lpm"] - d.flow_lpm) < 1e-9 else ""
        lines.append(f"| {r['flow_lpm']:.1f}{mark} | {r['dp_pa']/1e3:.2f} | {r['r_th_k_per_w']*1e3:.2f} | {r['base_avg_temp_c']:.1f} | "
                     f"{r['base_max_temp_c']:.1f} |")
    lines += [
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
        f"{req.base_thickness_m*1e3:.1f} mm base, k = {req.solid_k:.0f} W/m·K.",
        f"- Rated flow {d.flow_lpm:.1f} LPM ({d.flow_lpm/(req.heat_load_w/1e3):.1f} LPM per kW); the Pareto study is run at rated flow only. "
        f"ΔP budget {req.dp_budget_pa/1e3:g} kPa (channel core, at rated flow). R_th target "
        + (f"{target*1e3:g} mK/W." if target is not None else "none set."),
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
