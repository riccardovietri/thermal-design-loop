"""run_study: requirements.yaml in, Pareto plot, parity plot and report.md out. No LLM involved.

    uv run python run_study.py examples/gpu_coldplate.yaml
"""
import argparse
import json
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import yaml

from tools.check_confidence import check_confidence, mesh_independence
from tools.fluids import fluid_props
from tools.mesh import mesh
from tools.pareto import capacity_ratios, pareto_indices, recommend_marginal_return, spread_along_front
from tools.postprocess import plot_fields
from tools.report import parity_stats, plot_pareto, plot_parity, write_report
from tools.run_cfd import run_cfd
from tools.screen_1d import screen_1d
from tools.types import DesignSpace, Requirements, to_dict


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
                       _mk(r.get("r_th_threshold_mk_w")), _mk(r.get("r_th_objective_mk_w")), float(r.get("marginal_return_cutoff", 0.25)))
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
    shortlist = [front[i] for i in spread_along_front([s.dp_pa for s in front], opts["n_cfd"])]

    t0 = time.perf_counter()
    pairs, qualities = [], []
    for s in shortlist:
        run = runs_dir / s.design.design_id
        info = mesh(s.design, req.base_thickness_m, run, opts["mesh_cells_across_half_channel"])
        qualities.append(info.min_quality_sicn)
        pairs.append((s, run_cfd(s.design, req, fluid, info.path, run)))
    t_cfd = time.perf_counter() - t0

    cfd = [c for _, c in pairs]
    ok = [(s, c) for s, c in pairs if c.dp_pa <= req.dp_budget_pa]
    if not ok:
        raise SystemExit("No CFD-verified design meets the pressure-drop budget.")
    rec = recommend_marginal_return([(c.dp_pa, c.r_th_k_per_w) for _, c in pairs], req.dp_budget_pa, req.r_th_threshold_k_per_w,
                                    req.r_th_objective_k_per_w, req.marginal_return_cutoff)
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
    confidence = check_confidence(cfd, mesh_study, min(qualities))
    threshold = req.r_th_threshold_k_per_w
    ratios = {c.design_id: capacity_ratios(c.dp_pa, c.r_th_k_per_w, req.heat_load_w, s.design.flow_lpm / 6e4, threshold)
              if threshold is not None else None for s, c in pairs}

    study = {
        "name": name, "requirements": req, "fluid": fluid, "pairs": pairs, "recommended": recommended, "rule": rec["rule"],
        "selection": rec, "references": references, "ratios": ratios, "lowest_r_th": lowest, "flow_sensitivity": sensitivity,
        "n_sweep": len(sweep), "n_feasible": len(feasible), "t_screen_s": t_screen, "t_cfd_s": t_cfd, "confidence": confidence,
        "parity": parity_stats([s.r_th_k_per_w for s, _ in pairs], [c.r_th_k_per_w for _, c in pairs]),
        "dp_max_abs_pct": max(100 * abs(s.dp_pa / c.dp_pa - 1) for s, c in pairs),
    }
    plot_pareto(sweep, front, cfd, recommended[1].design_id, req.dp_budget_pa, out_dir / "pareto.png", req.r_th_threshold_k_per_w,
                req.r_th_objective_k_per_w)
    plot_parity(pairs, out_dir / "parity.png")
    plot_fields(runs_dir / recommended[1].design_id, out_dir / "fields_recommended.png")
    write_report(study, out_dir / "report.md")
    (out_dir / "results.json").write_text(json.dumps({
        "recommended": recommended[1].design_id, "rule": rec["rule"],
        "meets_threshold": rec["meets_threshold"], "meets_objective": rec["meets_objective"],
        "dp_budget_pa": req.dp_budget_pa, "r_th_threshold_k_per_w": threshold,
        "r_th_objective_k_per_w": req.r_th_objective_k_per_w, "marginal_return_cutoff": req.marginal_return_cutoff,
        "selection_steps": [{"from": pairs[st["from"]][1].design_id, "to": pairs[st["to"]][1].design_id if st["to"] is not None else None,
                             "returns": [{"to": pairs[e["to"]][1].design_id, "return": e["return"]} for e in st["returns"]]}
                            for st in rec["steps"]],
        "references": {k: (v[1].design_id if v is not None else None) for k, v in references.items()},
        "parity": study["parity"], "confidence": confidence, "flow_sensitivity": sensitivity,
        "designs": [{"one_d": to_dict(s), "cfd": to_dict(c), "ratios": ratios[c.design_id]} for s, c in pairs],
    }, indent=2))
    return study


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
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
    print(f"meets threshold: {study['selection']['meets_threshold']}, meets objective: {study['selection']['meets_objective']}")
    print(f"1D vs CFD R_th: mean |diff| {study['parity']['mean_abs_pct']:.1f}%, max {study['parity']['max_abs_pct']:.1f}%")
    print(f"confidence checks: {'all pass' if study['confidence']['all_pass'] else 'FAILURES, see report'}")
    print(f"report: {out / 'report.md'}")


if __name__ == "__main__":
    main()
