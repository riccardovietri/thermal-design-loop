"""check_confidence: mesh independence on one design plus per-run solver and energy checks."""
from pathlib import Path

from .mesh import mesh
from .run_cfd import run_cfd
from .types import CfdResult, Design, Fluid, Requirements

LIMITS = {"mesh_r_th_change_pct": 1.0, "mesh_dp_change_pct": 1.0, "linear_residual_rel": 1e-8, "base_flux_error_pct": 1.0}


def mesh_independence(design: Design, req: Requirements, fluid: Fluid, out_dir: Path, levels=(4, 8, 16)) -> dict:
    rows = []
    for cells in levels:
        run = Path(out_dir) / f"cells_{cells:02d}"
        info = mesh(design, req.base_thickness_m, run, cells_across_half_channel=cells)
        r = run_cfd(design, req, fluid, info.path, run)
        rows.append({"cells_across_half_channel": cells, "n_dofs": r.n_dofs, "min_quality_sicn": info.min_quality_sicn,
                     "r_th_k_per_w": r.r_th_k_per_w, "dp_pa": r.dp_pa})
    # production mesh is the middle level; compare it with the finest
    prod, fine = rows[len(rows) // 2], rows[-1]
    return {
        "design_id": design.design_id, "levels": rows,
        "r_th_change_pct": 100 * abs(prod["r_th_k_per_w"] / fine["r_th_k_per_w"] - 1),
        "dp_change_pct": 100 * abs(prod["dp_pa"] / fine["dp_pa"] - 1),
    }


def check_confidence(results: list[CfdResult], mesh_study: dict, min_mesh_quality: float) -> dict:
    worst_lin = max(r.linear_residual_rel for r in results)
    worst_flux = 100 * max(abs(r.base_flux_error_rel) for r in results)
    worst_re = max(r.reynolds for r in results)
    checks = [
        {"check": "Mesh independence, R_th (production vs 2× finer)", "value": f"{mesh_study['r_th_change_pct']:.3f}%",
         "limit": f"< {LIMITS['mesh_r_th_change_pct']}%", "pass": mesh_study["r_th_change_pct"] < LIMITS["mesh_r_th_change_pct"]},
        {"check": "Mesh independence, ΔP (production vs 2× finer)", "value": f"{mesh_study['dp_change_pct']:.3f}%",
         "limit": f"< {LIMITS['mesh_dp_change_pct']}%", "pass": mesh_study["dp_change_pct"] < LIMITS["mesh_dp_change_pct"]},
        {"check": "Linear solver residual, worst run", "value": f"{worst_lin:.1e}",
         "limit": f"< {LIMITS['linear_residual_rel']:.0e}", "pass": worst_lin < LIMITS["linear_residual_rel"]},
        {"check": "Energy balance: flux recovered at chip interface vs applied, worst run", "value": f"{worst_flux:.3f}%",
         "limit": f"< {LIMITS['base_flux_error_pct']}%", "pass": worst_flux < LIMITS["base_flux_error_pct"]},
        {"check": "Laminar assumption: highest Reynolds number", "value": f"{worst_re:.0f}", "limit": "< 2000", "pass": worst_re < 2000},
        {"check": "Mesh quality: lowest element quality (SICN)", "value": f"{min_mesh_quality:.2f}", "limit": "> 0.3", "pass": min_mesh_quality > 0.3},
    ]
    return {"checks": checks, "all_pass": all(c["pass"] for c in checks), "mesh_study": mesh_study}
