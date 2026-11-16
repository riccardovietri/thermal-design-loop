"""report: Pareto plot, 1D-vs-CFD parity plot and report.md for a completed study."""
from pathlib import Path

import numpy as np

from . import plotstyle as ps
from .plotstyle import plain_log_ticks
from .types import CfdResult, Screen1DResult


def parity_stats(one_d: list[float], cfd: list[float]) -> dict:
    diff = 100 * (np.asarray(one_d) / np.asarray(cfd) - 1)
    return {"mean_pct": float(diff.mean()), "mean_abs_pct": float(np.abs(diff).mean()), "max_abs_pct": float(np.abs(diff).max())}


RULE_TEXT = {
    "marginal_return": "at least one CFD-verified design meets both the budget and the threshold, so the marginal-return walk along "
                       "their front picks the design",
    "marginal_return_then_objective": "the marginal-return walk stopped at a design that misses the objective, and a design "
                                      "further along the eligible front meets it, so the rule continued to the lowest-ΔP design "
                                      "that meets the objective",
    "fallback_lowest_r_th": "no CFD-verified design meets the threshold within the budget, so the lowest R_th within budget is "
                            "recommended instead; **the threshold is not met**",
}
RULE_DEFINITION = (
    "Candidates are the CFD-verified designs with ΔP ≤ budget and R_th ≤ threshold; take their non-dominated front sorted by ΔP. "
    "Start at the lowest-ΔP front point. From the current point i, compute the return e(i→j) = −ln(R_j / R_i) / ln(ΔP_j / ΔP_i) "
    "to every front point j further along. If any j has e ≥ the cut-off, move to the nearest such j and repeat; otherwise stop "
    "and recommend the current point. A cut-off of 0.25 means a doubling of ΔP must buy at least about 16% lower R_th. "
    "Guard (SPEC Amendment 4): a point j with ΔP_j / ΔP_i below the minimum ratio ({min_ratio:.2f}) is ignored when looking from i, "
    "so two nearly coincident designs cannot produce a large return. Objective clause (SPEC Amendment 4): if the walk stops at a "
    "design that misses the objective and a design further along the eligible front meets it, the rule continues to the "
    "lowest-ΔP design that meets it. If no design meets the threshold within budget, the lowest R_th within budget is "
    "recommended and the threshold is reported as not met.")
RATIO_DEFINITIONS = (
    "Capacity = threshold × heat load / R_th, the heat load the design could carry at the threshold rise "
    "(base-to-inlet rise of threshold × heat load, {rise:g} K here); capacity per kPa = capacity / ΔP; "
    "kW per pumping W = heat load / (ΔP × volumetric flow), using hydraulic power. Both ratios always favour the lowest ΔP "
    "inside the window, which is why they are reported and not used to choose.")
X_MIN_KPA = 0.3


def draw_requirements(ax, x_lo: float, y_lo: float, budget_kpa: float, threshold_mk: float | None, objective_mk: float | None,
                      y_hi: float, labels_right: bool = False) -> None:
    """Budget, threshold and objective lines with short labels; the threshold-and-budget window is lightly shaded and the
    objective-and-budget window a little darker. labels_right puts the R_th labels just left of the budget line."""
    lx, ha = (budget_kpa * 0.96, "right") if labels_right else (x_lo * 1.04, "left")
    if threshold_mk is not None:
        ax.fill_between([x_lo, budget_kpa], y_lo, threshold_mk, color=ps.THRESHOLD, alpha=0.06, lw=0, zorder=0)
        ax.axhline(threshold_mk, color=ps.THRESHOLD, ls="--", lw=1.1, zorder=1)
        ax.text(lx, threshold_mk * 1.015, f"Threshold {threshold_mk:g} mK/W", color=ps.THRESHOLD, ha=ha, va="bottom",
                fontsize=ps.SMALL)
    if objective_mk is not None:
        ax.fill_between([x_lo, budget_kpa], y_lo, objective_mk, color=ps.THRESHOLD, alpha=0.08, lw=0, zorder=0)
        ax.axhline(objective_mk, color=ps.OBJECTIVE, ls=":", lw=1.4, zorder=1)
        ax.text(lx, objective_mk * 1.015, f"Objective {objective_mk:g} mK/W", color=ps.OBJECTIVE, ha=ha, va="bottom",
                fontsize=ps.SMALL)
    ax.axvline(budget_kpa, color=ps.ACCENT, ls="--", lw=1.1, zorder=1)
    ax.text(budget_kpa * 0.96, y_hi * 0.975, f"ΔP budget {budget_kpa:g} kPa", color=ps.ACCENT, ha="right", va="top",
            fontsize=ps.SMALL)


def plot_pareto(sweep: list[Screen1DResult], front: list[Screen1DResult], cfd: list[CfdResult], recommended_id: str,
                dp_budget_pa: float, out: Path, r_th_threshold: float | None = None, r_th_objective: float | None = None) -> Path:
    """R_th vs ΔP restricted to the region of interest: 0.3 kPa to 1.5× budget, 0.7× lowest feasible R_th to 1.5× threshold."""
    fig, ax = ps.figure(5.4)
    feas = [s for s in sweep if s.feasible]
    budget_kpa = dp_budget_pa / 1e3
    in_budget = [c.r_th_k_per_w for c in cfd if c.dp_pa <= dp_budget_pa] + [s.r_th_k_per_w for s in feas]
    y_lo = 0.7 * min(in_budget) * 1e3
    y_hi = 1.5 * r_th_threshold * 1e3 if r_th_threshold is not None else 1.3 * max(c.r_th_k_per_w for c in cfd) * 1e3
    x_lo, x_hi = X_MIN_KPA, 1.5 * budget_kpa
    mk = lambda v: None if v is None else v * 1e3
    draw_requirements(ax, x_lo, y_lo, budget_kpa, mk(r_th_threshold), mk(r_th_objective), y_hi)
    ax.scatter([s.dp_pa / 1e3 for s in feas], [s.r_th_k_per_w * 1e3 for s in feas], s=9, color=ps.FIN_GREY, alpha=0.6,
               label=f"1D screen, feasible designs ({len(feas)} of {len(sweep)})", linewidths=0, zorder=2)
    ax.plot([s.dp_pa / 1e3 for s in front], [s.r_th_k_per_w * 1e3 for s in front], color=ps.INK, lw=1.4, label="1D Pareto front",
            zorder=3)
    ax.scatter([c.dp_pa / 1e3 for c in cfd], [c.r_th_k_per_w * 1e3 for c in cfd], s=50, color=ps.BLUE, zorder=4,
               edgecolors="white", linewidths=0.8, label=f"2D conjugate CFD ({len(cfd)} designs)")
    rec = next(c for c in cfd if c.design_id == recommended_id)
    ax.scatter([rec.dp_pa / 1e3], [rec.r_th_k_per_w * 1e3], s=220, facecolors="none", edgecolors=ps.ACCENT, linewidths=2.2, zorder=5,
               label="Recommended design")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    plain_log_ticks(ax, x_lo, x_hi, "x")
    plain_log_ticks(ax, y_lo, y_hi, "y")
    ax.set_xlabel("Pressure drop ΔP (kPa)")
    ax.set_ylabel("Thermal resistance R_th (mK/W)")
    ax.set_title("Straight fins: thermal resistance vs pressure drop")
    ax.legend(loc="lower left", borderaxespad=0.6)
    ps.save(fig, out)
    return out


def plot_parity(pairs: list[tuple[Screen1DResult, CfdResult]], out: Path) -> Path:
    x = np.array([c.r_th_k_per_w * 1e3 for _, c in pairs])
    y = np.array([s.r_th_k_per_w * 1e3 for s, _ in pairs])
    eta = np.array([s.fin_efficiency for s, _ in pairs])
    stats = parity_stats(y, x)
    lo, hi = min(x.min(), y.min()) * 0.9, max(x.max(), y.max()) * 1.1
    fig, ax = ps.figure(6.2, width=7.0)
    line = np.array([lo, hi])
    ax.fill_between(line, line * 0.95, line * 1.05, color=ps.FAINT, alpha=0.6, lw=0, label="±5% band")
    ax.plot(line, line, color=ps.INK, lw=1, label="1D = CFD")
    sc = ax.scatter(x, y, c=eta, cmap="viridis", vmin=min(0.5, eta.min()), vmax=1.0, s=70, zorder=3, edgecolors="white", linewidths=0.8)
    ps.colorbar(fig, sc, ax, "Fin efficiency, 1D", min(0.5, eta.min()), 1.0, shrink=0.85)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    plain_log_ticks(ax, lo, hi, "x")
    plain_log_ticks(ax, lo, hi, "y")
    ax.set_xlabel("R_th, 2D conjugate CFD (mK/W)")
    ax.set_ylabel("R_th, 1D screen (mK/W)")
    ax.set_title(f"1D screen vs CFD: mean difference {stats['mean_abs_pct']:.1f}%, max {stats['max_abs_pct']:.1f}%", loc="left")
    ax.legend(loc="upper left")
    ps.save(fig, out)
    return out


def _geom(s: Screen1DResult) -> str:
    g = s.design
    return f"{g.w_ch_m*1e3:.2f} / {g.t_fin_m*1e3:.2f} / {g.h_fin_m*1e3:.2f} mm"


def _ret_text(e: dict, cutoff: float, dp_from: float, dp_to: float) -> str:
    if e.get("skipped"):
        return f"*skipped* ({e['return']:.3f}, ΔP ratio {dp_to / dp_from:.2f})"
    return f"**{e['return']:.3f}**" if e["return"] >= cutoff else f"{e['return']:.3f}"


def clause_lines(sel: dict, name) -> list[str]:
    """Plain-words outcome of the objective clause; name(i) describes candidate i."""
    if sel.get("objective_step"):
        st = sel["objective_step"]
        return [f"Objective clause: the walk stopped at {name(st['from'])}, which misses the objective; the rule continued to "
                f"{name(st['to'])}, the lowest-ΔP design on the eligible front that meets it."]
    if sel.get("objective_unreachable"):
        return [f"Objective clause: the walk stopped at {name(sel['walk_stop'])}, which misses the objective, and no design "
                "further along the eligible front meets it, so the rule stays there. **The objective is not met.**"]
    if sel.get("walk_stop") is not None:
        return ["Objective clause: not needed; the walk stopped at a design that meets the objective."]
    return []


def plot_cutoff_sensitivity(rows: list[dict], out: Path, current_cutoff: float, r_th_threshold: float | None,
                            r_th_objective: float | None, dp_budget_pa: float, title: str) -> Path:
    """Recommended R_th (top) and ΔP (bottom) against the marginal-return cut-off. rows: [{"cutoff", "r_th_k_per_w", "dp_pa"}]."""
    c = [r["cutoff"] for r in rows]
    fig, (a1, a2) = ps.figure(5.6, 2, 1, sharex=True)
    a1.plot(c, [r["r_th_k_per_w"] * 1e3 for r in rows], color=ps.BLUE, marker="o", drawstyle="steps-mid", zorder=3)
    a2.plot(c, [r["dp_pa"] / 1e3 for r in rows], color=ps.BLUE, marker="o", drawstyle="steps-mid", zorder=3)
    y = [r["r_th_k_per_w"] * 1e3 for r in rows]
    lo, hi = min(y), max(y)
    for v, colour, ls, label in ((r_th_threshold, ps.THRESHOLD, "--", "Threshold"), (r_th_objective, ps.OBJECTIVE, ":", "Objective")):
        if v is not None:
            a1.axhline(v * 1e3, color=colour, ls=ls, lw=1.3, zorder=1)
            a1.text(c[-1], v * 1e3, f"{label} {v*1e3:g} mK/W", color=colour, ha="right", va="bottom", fontsize=ps.SMALL)
            lo, hi = min(lo, v * 1e3), max(hi, v * 1e3)
    pad = 0.08 * (hi - lo or 1.0)
    a1.set_ylim(lo - pad, hi + 3 * pad)
    dp = [r["dp_pa"] / 1e3 for r in rows]
    a2.set_ylim(0, max(max(dp), dp_budget_pa / 1e3) * 1.12)
    a2.axhline(dp_budget_pa / 1e3, color=ps.ACCENT, ls="--", lw=1.1, zorder=1)
    a2.text(c[-1], dp_budget_pa / 1e3, f"ΔP budget {dp_budget_pa/1e3:g} kPa", color=ps.ACCENT, ha="right", va="bottom",
            fontsize=ps.SMALL)
    for ax in (a1, a2):
        ax.axvline(current_cutoff, color=ps.MUTED, ls="-", lw=1.0, zorder=1)
    a1.text(current_cutoff, a1.get_ylim()[1], f" cut-off in use {current_cutoff:g}", color=ps.MUTED, ha="left", va="top",
            fontsize=ps.SMALL)
    a1.set_ylabel("Recommended R_th (mK/W)")
    a2.set_ylabel("Recommended ΔP (kPa)")
    a2.set_xlabel("Marginal-return cut-off")
    a2.set_xticks(c)
    a2.xaxis.set_major_formatter(lambda v, _: f"{v:.2f}")
    a2.tick_params(axis="x", labelsize=ps.SMALL)
    fig.suptitle(title)
    ps.save(fig, out)
    return out


def plot_local_returns(dp_front: list[float], r_front: list[float], out: Path, cutoff: float, min_dp_ratio: float,
                       recommended_dp: float, r_th_threshold: float | None, title: str) -> Path:
    """Return between consecutive points of a front against ΔP: one horizontal segment per pair, from ΔP_i to ΔP_j."""
    from .pareto import local_returns
    lr = local_returns(dp_front, r_front)
    fig, ax = ps.figure(4.8)
    for e in lr:
        x = [e["x_i"] / 1e3, e["x_j"] / 1e3]
        eligible = r_th_threshold is None or r_front[e["i"]] <= r_th_threshold
        colour = ps.BLUE if eligible else ps.FIN_GREY
        close = e["dp_ratio"] < min_dp_ratio
        ax.plot(x, [e["return"]] * 2, color=colour, lw=2.4, solid_capstyle="butt", zorder=3, ls=(0, (2, 1.5)) if close else "-")
        ax.scatter([np.sqrt(x[0] * x[1])], [e["return"]], s=34, zorder=4, edgecolors=colour, linewidths=1.2,
                   facecolors="white" if close else colour)
    ax.axhline(cutoff, color=ps.ACCENT, ls="--", lw=1.2, zorder=2)
    xs = [d / 1e3 for d in dp_front]
    ax.text(xs[-1], cutoff, f"cut-off {cutoff:g}", color=ps.ACCENT, ha="right", va="bottom", fontsize=ps.SMALL)
    ax.axvline(recommended_dp / 1e3, color=ps.INK, ls=":", lw=1.2, zorder=2)
    ax.text(recommended_dp / 1e3, 0.98, " recommended", color=ps.INK, ha="left", va="top", fontsize=ps.SMALL,
            transform=ax.get_xaxis_transform())
    ax.set_xscale("log")
    lo, hi = xs[0] * 0.9, xs[-1] * 1.1
    ax.set_xlim(lo, hi)
    plain_log_ticks(ax, lo, hi, "x")
    ymax = max([e["return"] for e in lr] + [cutoff])
    ax.set_ylim(min([0.0] + [e["return"] for e in lr]), ymax * 1.15)
    ax.set_xlabel("Pressure drop ΔP (kPa)")
    ax.set_ylabel("Return between neighbours")
    ax.set_title(title)
    ps.save(fig, out)
    return out


def range_text(rng) -> str:
    """'for cut-offs from a to b', or 'only at a (...)' when the stable run is a single swept value."""
    lo, hi = rng
    if abs(hi - lo) < 1e-9:
        return f"only at a cut-off of {lo:.2f}; the neighbouring swept values (±0.05) pick a different design"
    return f"for cut-offs from {lo:.2f} to {hi:.2f}"


def cutoff_section(study: dict) -> list[str]:
    """Report section: cut-off sensitivity table and figure, and the local-return figure."""
    rows, cutoff = study.get("cutoff_sensitivity"), study["requirements"].marginal_return_cutoff
    if not rows:
        return []
    rng = study["cutoff_stable_range"]
    pairs = study["pairs"]
    lines = ["## Sensitivity to the cut-off", "",
             "The rule re-run on the same CFD-verified designs for cut-offs from 0.05 to 0.60 (guard and objective clause on "
             "throughout). No new solves.", "",
             "| Cut-off | Recommended (channel / fin / height) | R_th (mK/W) | ΔP (kPa) | Rule | Walk stopped at |",
             "|---|---|---|---|---|---|"]
    for r in rows:
        s_, c_ = pairs[r["index"]]
        stop = _geom(pairs[r["walk_stop"]][0]) if r["walk_stop"] is not None and r["walk_stop"] != r["index"] else "same"
        mark = " (in use)" if abs(r["cutoff"] - cutoff) < 1e-9 else ""
        lines.append(f"| {r['cutoff']:.2f}{mark} | {_geom(s_)} | {c_.r_th_k_per_w*1e3:.2f} | {c_.dp_pa/1e3:.2f} | `{r['rule']}` | {stop} |")
    lines += ["", "![Cut-off sensitivity](cutoff_sensitivity.png)", "",
              "Recommended R_th (top) and ΔP (bottom) against the marginal-return cut-off; the vertical line is the cut-off in use. "
              + (f"The current recommendation is unchanged {range_text(rng)}." if rng else ""), ""]
    lf = study["local_front"]
    lines += ["## Local return along the front", "",
              "![Local marginal return](local_return.png)", "",
              f"Return e between neighbouring points of the CFD-verified front within budget ({len(lf)} points), drawn as a segment "
              "from the lower to the higher ΔP of each pair. Blue: pairs whose lower-ΔP point meets the threshold; grey: it does "
              f"not. Hollow markers and dashed segments: neighbours closer than the guard's ΔP ratio of {study['requirements'].min_dp_ratio:.2f}, "
              "which the rule ignores. The dashed line is the cut-off in use; the dotted line marks the recommended design. "
              "The rule compares the current point with every further point, not only its neighbour, so a single segment below "
              "the cut-off does not by itself stop the walk.", ""]
    return lines


def write_report(study: dict, out: Path) -> Path:
    req, fluid, pairs = study["requirements"], study["fluid"], study["pairs"]
    rec_s, rec_c = study["recommended"]
    rule, sel = study["rule"], study["selection"]
    threshold, objective, cutoff = req.r_th_threshold_k_per_w, req.r_th_objective_k_per_w, req.marginal_return_cutoff
    stats, conf, ratios = study["parity"], study["confidence"], study["ratios"]
    d = rec_s.design
    fmt_req = lambda v: f"{v*1e3:g} mK/W" if v is not None else "none set"
    lines = [
        f"# Cold-plate design study: {study['name']}",
        "",
        f"**Recommended design:** {d.w_ch_m*1e3:.2f} mm channels, {d.t_fin_m*1e3:.2f} mm fins, {d.h_fin_m*1e3:.2f} mm tall "
        f"({rec_c.n_channels} channels) at {d.flow_lpm:.1f} LPM (rated flow). CFD gives R_th = {rec_c.r_th_k_per_w*1e3:.1f} mK/W at "
        f"ΔP = {rec_c.dp_pa/1e3:.2f} kPa, so the chip interface averages {rec_c.base_avg_temp_c:.1f} °C "
        f"(peak {rec_c.base_max_temp_c:.1f} °C) for {req.heat_load_w:.0f} W and {req.inlet_temp_c:.0f} °C inlet. "
        f"It {'meets' if sel['meets_threshold'] else '**does not meet**'} the {fmt_req(threshold)} threshold and "
        f"{'meets' if sel['meets_objective'] else '**does not meet**'} the {fmt_req(objective)} objective.",
        "",
        "## Why this design",
        f"- Requirements: channel-core ΔP ≤ {req.dp_budget_pa/1e3:g} kPa at rated flow (budget); R_th ≤ {fmt_req(threshold)} "
        f"(threshold, must meet); R_th ≤ {fmt_req(objective)} (objective, want).",
        f"- Rule (marginal return, cut-off {cutoff:g}, minimum ΔP ratio {req.min_dp_ratio:.2f}): "
        + RULE_DEFINITION.format(min_ratio=req.min_dp_ratio),
        f"- Rule that fired: `{rule}`: {RULE_TEXT[rule]}.",
    ]
    if sel["steps"]:
        lines += ["", "Steps the rule took (returns from the current point to every further point on the eligible front; "
                  f"**bold** = at or above the cut-off {cutoff:g}; *skipped* = ΔP ratio below the "
                  f"guard's {req.min_dp_ratio:.2f}):", "",
                  "| Step | From (channel / fin / height) | R_th (mK/W) | ΔP (kPa) | Returns to further points | Moved to |",
                  "|---|---|---|---|---|---|"]
        for k, st in enumerate(sel["steps"], 1):
            fs, fc = pairs[st["from"]]
            rets = "; ".join(_ret_text(e, cutoff, pairs[st["from"]][1].dp_pa, pairs[e["to"]][1].dp_pa)
                             + f" to {_geom(pairs[e['to']][0])}" for e in st["returns"]) or "none (last front point)"
            moved = _geom(pairs[st["to"]][0]) if st["to"] is not None else "stop: recommended"
            lines.append(f"| {k} | {_geom(fs)} | {fc.r_th_k_per_w*1e3:.2f} | {fc.dp_pa/1e3:.2f} | {rets} | {moved} |")
    lines += [""] + clause_lines(sel, lambda i: f"{_geom(pairs[i][0])} ({pairs[i][1].r_th_k_per_w*1e3:.2f} mK/W at "
                                                f"{pairs[i][1].dp_pa/1e3:.2f} kPa)")
    lines += ["", "Reference designs (CFD-verified, within budget):", "",
              "| Reference | Channel / fin / height | R_th (mK/W) | ΔP (kPa) |", "|---|---|---|---|"]
    names = {"lowest_dp_meeting_threshold": f"Lowest ΔP meeting the threshold ({fmt_req(threshold)})",
             "lowest_dp_meeting_objective": f"Lowest ΔP meeting the objective ({fmt_req(objective)})",
             "lowest_r_th_within_budget": "Lowest R_th within budget"}
    for key, label in names.items():
        ref = study["references"][key]
        if ref is None:
            lines.append(f"| {label} | none | — | — |")
        else:
            mark = " (recommended)" if ref[1].design_id == rec_c.design_id else ""
            lines.append(f"| {label}{mark} | {_geom(ref[0])} | {ref[1].r_th_k_per_w*1e3:.2f} | {ref[1].dp_pa/1e3:.2f} |")
    lines += [
        "",
        f"- At the recommended point the resistance splits into base conduction {rec_s.r_base*1e3:.1f}, fin and convection {rec_s.r_conv*1e3:.1f}, "
        f"and coolant heating {rec_s.r_caloric*1e3:.1f} mK/W (1D breakdown).",
        "",
        "![Pareto](pareto.png)",
        "",
        f"Thermal resistance against channel pressure drop at the rated {d.flow_lpm:.1f} LPM. Grey dots: every feasible design from "
        "the 1D screen; dark line: its Pareto front; blue: the designs solved with 2D conjugate CFD; ring: the recommended design. "
        "The plot shows the region of interest only: ΔP from 0.3 kPa to 1.5× the budget, R_th from 0.7× the lowest feasible value "
        "to 1.5× the threshold. Designs outside it are screened but not drawn. Light shading: meets the budget and the threshold; "
        "darker: meets the budget and the objective.",
        "",
    ]
    lines += cutoff_section(study)
    lines += [
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
        f"The 1D screen evaluated {study['n_sweep']} designs ({study['n_feasible']} feasible) in {study['t_screen_s']*1e3:.0f} ms. " +
        (f"Every design on its Pareto front within budget ({len(pairs)}) was solved with the 2D conjugate model"
         if study.get("solved_whole_front") else f"{len(pairs)} designs spread along its Pareto front were solved with the 2D "
         "conjugate model") + f" in {study['t_cfd_s']:.1f} s.",
        "",
        f"- R_th: 1D differs from CFD by {stats['mean_abs_pct']:.1f}% on average, {stats['max_abs_pct']:.1f}% at worst (mean signed {stats['mean_pct']:+.1f}%).",
        f"- ΔP: identical to within {study['dp_max_abs_pct']:.2f}%. This is expected, not a result: both models assume fully developed flow and share the same f·Re.",
        "",
        "![Parity](parity.png)",
        "",
        "R_th from the 1D screen against the 2D conjugate CFD for each solved design, coloured by the 1D fin efficiency; "
        "the band is ±5%.",
        "",
        "Ratios (CFD): " + RATIO_DEFINITIONS.format(rise=(threshold or 0) * req.heat_load_w),
        "",
        "| Channel (mm) | Fin (mm) | Height (mm) | Flow (LPM) | Re | ΔP CFD (kPa) | R_th 1D (mK/W) | R_th CFD (mK/W) | 1D − CFD | Fin eff. "
        "| Capacity (kW) | kW per kPa | kW per pumping W |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for s, c in pairs:
        g = s.design
        mark = " **(recommended)**" if c.design_id == rec_c.design_id else ""
        lines.append(f"| {g.w_ch_m*1e3:.2f}{mark} | {g.t_fin_m*1e3:.2f} | {g.h_fin_m*1e3:.2f} | {g.flow_lpm:.1f} | {c.reynolds:.0f} | {c.dp_pa/1e3:.2f} | "
                     f"{s.r_th_k_per_w*1e3:.2f} | {c.r_th_k_per_w*1e3:.2f} | {100*(s.r_th_k_per_w/c.r_th_k_per_w-1):+.1f}% | {s.fin_efficiency:.2f} | "
                     + (f"{ratios[c.design_id]['capacity_kw']:.2f} | {ratios[c.design_id]['capacity_kw_per_kpa']:.2f} | "
                        f"{ratios[c.design_id]['kw_per_pumping_w']:.1f} |" if ratios.get(c.design_id) else "— | — | — |"))
    lines += ["", "## Confidence checks", "", "| Check | Value | Limit | Result |", "|---|---|---|---|"]
    lines += [f"| {c['check']} | {c['value']} | {c['limit']} | {'pass' if c['pass'] else 'FAIL'} |" for c in conf["checks"]]
    ms = conf["mesh_study"]
    lines += ["", f"Mesh study on `{ms['design_id']}`:", "", "| Cells across half channel | Unknowns | R_th (mK/W) | ΔP (kPa) |", "|---|---|---|---|"]
    lines += [f"| {r['cells_across_half_channel']} | {r['n_dofs']} | {r['r_th_k_per_w']*1e3:.4f} | {r['dp_pa']/1e3:.4f} |" for r in ms["levels"]]
    lines += [
        "", "![Fields](fields_recommended.png)", "",
        f"Recommended design at {d.flow_lpm:.1f} LPM {fluid.name}: axial velocity (left; copper in brown) and temperature at "
        f"mid-length (right) on the cross-section, mirrored to show three channels. Re = {rec_c.reynolds:.0f}, "
        f"ΔP = {rec_c.dp_pa/1e3:.2f} kPa, R_th = {rec_c.r_th_k_per_w*1e3:.1f} mK/W.", "",
        "## Inputs",
        f"- Heat load {req.heat_load_w:.0f} W, uniform over {req.footprint_width_m*1e3:.0f} × {req.flow_length_m*1e3:.0f} mm; "
        f"{req.base_thickness_m*1e3:.1f} mm base, k = {req.solid_k:.0f} W/m·K.",
        f"- Rated flow {d.flow_lpm:.1f} LPM ({d.flow_lpm/(req.heat_load_w/1e3):.1f} LPM per kW); the Pareto study is run at rated flow only. "
        f"ΔP budget {req.dp_budget_pa/1e3:g} kPa (channel core, at rated flow; one third of a 15 kPa whole-plate allocation). "
        f"R_th threshold {fmt_req(threshold)}, objective {fmt_req(objective)}, marginal-return cut-off {cutoff:g}.",
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
