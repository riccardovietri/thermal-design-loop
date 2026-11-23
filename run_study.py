"""run_study: requirements.yaml in, Pareto plot, parity plot and report.md out. No LLM involved.

    uv run python run_study.py examples/gpu_coldplate.yaml
"""
import argparse
import json
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import yaml

from tools.check_confidence import check_confidence, mesh_independence, run_failures
from tools.fluids import fluid_props
from tools.mesh import mesh
from tools.pareto import (capacity_ratios, cutoff_sensitivity, pareto_indices, recommend_marginal_return_objective,
                          spread_along_front, stable_cutoff_range)
from tools.postprocess import plot_fields
from tools.report import (parity_stats, plot_cutoff_sensitivity, plot_local_returns, plot_pareto, plot_parity,
                          write_report)
from tools.run_cfd import run_cfd
from tools.screen_1d import screen_1d
from tools.types import CfdResult, Design, DesignSpace, Requirements, Screen1DResult, to_dict


def _levels(spec) -> tuple[float, ...]:
    if isinstance(spec, list):
        return tuple(float(v) for v in spec)
    fn = np.geomspace if spec.get("spacing") == "geometric" else np.linspace
    return tuple(float(v) for v in fn(spec["min"], spec["max"], spec["levels"]))


def _mk(v):
    return None if v is None else v / 1e3


def load(path: Path) -> tuple[str, Requirements, DesignSpace, dict]:
    cfg = yaml.safe_load(Path(path).read_text())
    r, s, c = cfg["requirements"], cfg["design_space"], cfg["constraints"]
    req = Requirements(r["heat_load_w"], r["footprint_width_mm"] / 1e3, r["flow_length_mm"] / 1e3, r["base_thickness_mm"] / 1e3,
                       r["solid_k_w_mk"], r["inlet_temp_c"], r["dp_budget_kpa"] * 1e3, r["coolant"],
                       _mk(r.get("r_th_threshold_mk_w")), _mk(r.get("r_th_objective_mk_w")), float(r.get("marginal_return_cutoff", 0.25)),
                       float(r.get("min_dp_ratio", 1.10)))
    mm = lambda key: tuple(v / 1e3 for v in _levels(s[key]))
    space = DesignSpace(mm("w_ch_mm"), mm("t_fin_mm"), mm("h_fin_mm"), _levels(s["flow_lpm"]), c["min_feature_mm"] / 1e3, c["max_aspect_ratio"])
    return cfg["name"], req, space, cfg["study"]


def run_study(config: Path, out_dir: Path | None = None, runs_dir: Path | None = None) -> dict:
    name, req, space, opts = load(config)
    out_dir = Path(out_dir or Path("out") / name)
    runs_dir = Path(runs_dir or Path("runs") / name)
    out_dir.mkdir(parents=True, exist_ok=True)
    fluid = fluid_props(req.coolant, req.inlet_temp_c)

    t0 = time.perf_counter()
    sweep = screen_1d(req, fluid, space)
    t_screen = time.perf_counter() - t0
    feasible = [s for s in sweep if s.feasible]
    if not feasible:
        raise SystemExit("No feasible design in the 1D sweep; relax the constraints or widen the design space.")
    front = [feasible[i] for i in pareto_indices([s.dp_pa for s in feasible], [s.r_th_k_per_w for s in feasible])]
    whole_front = opts["n_cfd"] in ("all", None)
    shortlist = front if whole_front else [front[i] for i in spread_along_front([s.dp_pa for s in front], opts["n_cfd"])]

    t0 = time.perf_counter()
    pairs, qualities, excluded, solved = [], [], [], []
    for s in shortlist:
        run = runs_dir / s.design.design_id
        info = mesh(s.design, req.base_thickness_m, run, opts["mesh_cells_across_half_channel"])
        qualities.append(info.min_quality_sicn)
        c = run_cfd(s.design, req, fluid, info.path, run)
        solved.append(c)
        reasons = run_failures(c, info.min_quality_sicn)
        # a design whose own run fails a check is not a candidate for selection (issue #21)
        if reasons:
            excluded.append((s, c, reasons))
        else:
            pairs.append((s, c))
    t_cfd = time.perf_counter() - t0
    if not pairs:
        raise SystemExit("Every CFD design failed a per-run confidence check; no candidate for selection: "
                         + "; ".join(f"{c.design_id}: {', '.join(r)}" for _, c, r in excluded))

    cfd = [c for _, c in pairs]
    ok = [(s, c) for s, c in pairs if c.dp_pa <= req.dp_budget_pa]
    if not ok:
        raise SystemExit("No CFD-verified design meets the pressure-drop budget.")
    rec = recommend_marginal_return_objective([(c.dp_pa, c.r_th_k_per_w) for _, c in pairs], req.dp_budget_pa,
                                              req.r_th_threshold_k_per_w, req.r_th_objective_k_per_w,
                                              req.marginal_return_cutoff, req.min_dp_ratio)
    cut = selection_extras(pairs, req)
    recommended = pairs[rec["index"]]
    lowest = min(ok, key=lambda p: p[1].r_th_k_per_w)

    def lowest_dp_meeting(limit):
        meets = [p for p in ok if limit is not None and p[1].r_th_k_per_w <= limit]
        return min(meets, key=lambda p: p[1].dp_pa) if meets else None

    references = {"lowest_dp_meeting_threshold": lowest_dp_meeting(req.r_th_threshold_k_per_w),
                  "lowest_dp_meeting_objective": lowest_dp_meeting(req.r_th_objective_k_per_w),
                  "lowest_r_th_within_budget": lowest}

    cells = opts["mesh_cells_across_half_channel"]
    rec_design = recommended[0].design
    flows = sorted(set(opts.get("flow_sensitivity_lpm", [])) | {rec_design.flow_lpm})
    sens_dir = runs_dir / "_flow_sensitivity"
    sens_mesh = mesh(rec_design, req.base_thickness_m, sens_dir / "mesh", cells)
    sensitivity = []
    for q in flows:
        d = replace(rec_design, flow_lpm=q)
        c = run_cfd(d, req, fluid, sens_mesh.path, sens_dir / d.design_id)
        sensitivity.append({"flow_lpm": q, "design_id": d.design_id, "dp_pa": c.dp_pa, "r_th_k_per_w": c.r_th_k_per_w,
                            "base_avg_temp_c": c.base_avg_temp_c, "base_max_temp_c": c.base_max_temp_c, "reynolds": c.reynolds})
    mesh_study = mesh_independence(rec_design, req, fluid, runs_dir / "_mesh_study", (cells // 2, cells, cells * 2))
    confidence = check_confidence(solved, mesh_study, min(qualities))
    failed_checks = [c["check"] for c in confidence["checks"] if not c["pass"]]
    threshold = req.r_th_threshold_k_per_w
    ratios = {c.design_id: capacity_ratios(c.dp_pa, c.r_th_k_per_w, req.heat_load_w, s.design.flow_lpm / 6e4, threshold)
              if threshold is not None else None for s, c in pairs}

    study = {
        "name": name, "requirements": req, "fluid": fluid, "pairs": pairs, "recommended": recommended, "rule": rec["rule"],
        "selection": rec, "references": references, "ratios": ratios, "lowest_r_th": lowest, "flow_sensitivity": sensitivity,
        "solved_whole_front": whole_front, **cut, "failed_checks": failed_checks, "excluded": excluded,
        "n_sweep": len(sweep), "n_feasible": len(feasible), "t_screen_s": t_screen, "t_cfd_s": t_cfd, "confidence": confidence,
        "parity": parity_stats([s.r_th_k_per_w for s, _ in pairs], [c.r_th_k_per_w for _, c in pairs]),
        "dp_max_abs_pct": max(100 * abs(s.dp_pa / c.dp_pa - 1) for s, c in pairs),
    }
    plot_pareto(sweep, front, cfd, recommended[1].design_id, req.dp_budget_pa, out_dir / "pareto.png", req.r_th_threshold_k_per_w,
                req.r_th_objective_k_per_w)
    plot_parity(pairs, out_dir / "parity.png")
    plot_fields(runs_dir / recommended[1].design_id, out_dir / "fields_recommended.png")
    plot_selection(cut, pairs, recommended[1], req, out_dir)
    write_report(study, out_dir / "report.md")
    # Keys of a blocked study come first; a passing study writes exactly what it wrote before issue #21.
    blocked = {"status": "blocked", "failed_checks": failed_checks,
               "excluded_designs": [{"design_id": c.design_id, "reasons": r, "one_d": to_dict(s), "cfd": to_dict(c)}
                                    for s, c, r in excluded]} if failed_checks else {}
    (out_dir / "results.json").write_text(json.dumps({
        **blocked, "recommended": recommended[1].design_id, "rule": rec["rule"],
        "meets_threshold": rec["meets_threshold"], "meets_objective": rec["meets_objective"],
        "dp_budget_pa": req.dp_budget_pa, "r_th_threshold_k_per_w": threshold,
        "r_th_objective_k_per_w": req.r_th_objective_k_per_w, "marginal_return_cutoff": req.marginal_return_cutoff,
        "min_dp_ratio": req.min_dp_ratio,
        "selection_steps": [{"from": pairs[st["from"]][1].design_id, "to": pairs[st["to"]][1].design_id if st["to"] is not None else None,
                             "returns": [{"to": pairs[e["to"]][1].design_id, "return": e["return"], "skipped": e["skipped"]}
                                         for e in st["returns"]]}
                            for st in rec["steps"]],
        "walk_stop": pairs[rec["walk_stop"]][1].design_id if rec["walk_stop"] is not None else None,
        "objective_step": ({k: pairs[v][1].design_id for k, v in rec["objective_step"].items()}
                           if rec["objective_step"] else None),
        "objective_unreachable": rec["objective_unreachable"],
        "cutoff_sensitivity": [{**{k: r[k] for k in ("cutoff", "rule", "meets_objective")},
                                "recommended": pairs[r["index"]][1].design_id,
                                "walk_stop": pairs[r["walk_stop"]][1].design_id if r["walk_stop"] is not None else None,
                                "r_th_k_per_w": pairs[r["index"]][1].r_th_k_per_w, "dp_pa": pairs[r["index"]][1].dp_pa}
                               for r in cut["cutoff_sensitivity"]],
        "cutoff_stable_range": cut["cutoff_stable_range"],
        "references": {k: (v[1].design_id if v is not None else None) for k, v in references.items()},
        "parity": study["parity"], "confidence": confidence, "flow_sensitivity": sensitivity,
        "designs": [{"one_d": to_dict(s), "cfd": to_dict(c), "ratios": ratios[c.design_id]} for s, c in pairs],
    }, indent=2))
    return study


def selection_extras(pairs, req) -> dict:
    """Cut-off sensitivity of the rule over the CFD-verified designs, and the CFD front within budget for the local-return plot."""
    cands = [(c.dp_pa, c.r_th_k_per_w) for _, c in pairs]
    rows = cutoff_sensitivity(cands, req.dp_budget_pa, req.r_th_threshold_k_per_w, req.r_th_objective_k_per_w,
                              min_dp_ratio=req.min_dp_ratio)
    inb = [i for i, (dp, _) in enumerate(cands) if dp <= req.dp_budget_pa]
    front = [inb[k] for k in pareto_indices([cands[i][0] for i in inb], [cands[i][1] for i in inb])]
    return {"cutoff_sensitivity": rows, "cutoff_stable_range": stable_cutoff_range(rows, req.marginal_return_cutoff),
            "local_front": front}


def plot_selection(cut: dict, pairs, recommended, req, out_dir: Path) -> list[Path]:
    rows = [{"cutoff": r["cutoff"], "r_th_k_per_w": pairs[r["index"]][1].r_th_k_per_w, "dp_pa": pairs[r["index"]][1].dp_pa}
            for r in cut["cutoff_sensitivity"]]
    f = cut["local_front"]
    return [plot_cutoff_sensitivity(rows, out_dir / "cutoff_sensitivity.png", req.marginal_return_cutoff, req.r_th_threshold_k_per_w,
                                    req.r_th_objective_k_per_w, req.dp_budget_pa, "Straight fins: recommended design against cut-off"),
            plot_local_returns([pairs[i][1].dp_pa for i in f], [pairs[i][1].r_th_k_per_w for i in f], out_dir / "local_return.png",
                               req.marginal_return_cutoff, req.min_dp_ratio, recommended.dp_pa, req.r_th_threshold_k_per_w,
                               "Straight fins: return between neighbouring front points")]


def plot_only(config: Path, out_dir: Path | None = None, runs_dir: Path | None = None) -> list[Path]:
    """Redraw pareto.png, parity.png and fields_recommended.png from a finished study's results.json and run directories,
    without any CFD, plus cutoff_sensitivity.png and local_return.png. The 1D screen (analytic, milliseconds) is re-run because
    its designs are not saved."""
    name, req, space, _ = load(config)
    out_dir = Path(out_dir or Path("out") / name)
    runs_dir = Path(runs_dir or Path("runs") / name)
    saved = json.loads((out_dir / "results.json").read_text())
    sweep = screen_1d(req, fluid_props(req.coolant, req.inlet_temp_c), space)
    feasible = [s for s in sweep if s.feasible]
    front = [feasible[i] for i in pareto_indices([s.dp_pa for s in feasible], [s.r_th_k_per_w for s in feasible])]
    pairs = [(Screen1DResult(**{**d["one_d"], "design": Design(**d["one_d"]["design"])}), CfdResult(**d["cfd"])) for d in saved["designs"]]
    cfd = [c for _, c in pairs]
    rec = next(c for c in cfd if c.design_id == saved["recommended"])
    return plot_selection(selection_extras(pairs, req), pairs, rec, req, out_dir) + [
            plot_pareto(sweep, front, cfd, saved["recommended"], req.dp_budget_pa, out_dir / "pareto.png", req.r_th_threshold_k_per_w,
                        req.r_th_objective_k_per_w),
            plot_parity(pairs, out_dir / "parity.png"),
            plot_fields(runs_dir / saved["recommended"], out_dir / "fields_recommended.png")]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--plot-only", action="store_true", help="redraw the figures from results.json and runs/, no CFD")
    a = ap.parse_args()
    if a.plot_only:
        for p in plot_only(a.config, a.out):
            print(p)
        return 0
    t0 = time.perf_counter()
    study = run_study(a.config, a.out)
    rec = study["recommended"][1]
    out = a.out or Path("out") / study["name"]
    print(f"{study['n_sweep']} designs screened, {len(study['pairs'])} solved with CFD in {time.perf_counter()-t0:.1f} s total")
    print(f"recommended {rec.design_id} ({study['rule']}): R_th {rec.r_th_k_per_w*1e3:.2f} mK/W at {rec.dp_pa/1e3:.2f} kPa")
    for st in study["selection"]["steps"]:
        frm = study["pairs"][st["from"]][1].design_id
        rets = ", ".join(f"{study['pairs'][e['to']][1].design_id} {e['return']:.3f}" for e in st["returns"]) or "none"
        print(f"  from {frm}: returns {rets}")
    if study["selection"]["objective_step"]:
        st = study["selection"]["objective_step"]
        print(f"  objective clause: walk stopped at {study['pairs'][st['from']][1].design_id}; continued to "
              f"{study['pairs'][st['to']][1].design_id}")
    rng = study["cutoff_stable_range"]
    print(f"recommendation unchanged for cut-offs {rng[0]:.2f}-{rng[1]:.2f}" if rng else "cut-off in use not in the sweep")
    print(f"meets threshold: {study['selection']['meets_threshold']}, meets objective: {study['selection']['meets_objective']}")
    print(f"1D vs CFD R_th: mean |diff| {study['parity']['mean_abs_pct']:.1f}%, max {study['parity']['max_abs_pct']:.1f}%")
    print(f"confidence checks: {'all pass' if study['confidence']['all_pass'] else 'FAILURES, see report'}")
    print(f"report: {out / 'report.md'}")
    if study["failed_checks"]:
        print("BLOCKED: confidence checks failed: " + "; ".join(study["failed_checks"]), file=sys.stderr)
        for _, c, r in study["excluded"]:
            print(f"  excluded {c.design_id}: {', '.join(r)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
