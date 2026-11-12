"""Level 2 case: array of confined slot jets impinging on a flat heated base, with return slots in between.

Unit cell (lengths / slot width b): x from the jet centreline (0) to the return-slot centreline (P/2), y from the base (0)
to the cover (H). Jet enters through the cover over 0 <= x <= 0.5, leaves through the cover over P/2 - 0.5 <= x <= P/2.
"""
import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import gmsh
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.tri import Triangulation
from skfem import MeshTri

from tools.fluids import fluid_props
from tools.types import Fluid, Requirements

from .planar import boundary_profile, solve_energy, solve_flow, stream_function


@dataclass(frozen=True)
class SlotJetDesign:
    slot_width_m: float
    gap_m: float
    pitch_m: float
    flow_lpm: float

    @property
    def design_id(self) -> str:
        return f"jet_b{self.slot_width_m*1e6:04.0f}_g{self.gap_m*1e6:04.0f}_p{self.pitch_m*1e6:05.0f}_q{self.flow_lpm*100:03.0f}"


@dataclass(frozen=True)
class SlotJetResult:
    design_id: str
    n_slots: int
    reynolds_slot: float
    jet_velocity_m_s: float
    dp_pa: float
    r_th_k_per_w: float
    r_conv_k_per_w: float
    r_base_k_per_w: float
    h_stagnation_w_m2k: float
    h_mean_w_m2k: float
    base_avg_temp_c: float
    base_max_temp_c: float
    energy_balance_error: float
    newton_iterations: int
    n_unknowns: int


def _mesh(half_pitch: float, gap: float, size: float) -> MeshTri:
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        geo = gmsh.model.geo
        fine, coarse = size, 2.5 * size
        pts = [geo.addPoint(0, 0, 0, fine), geo.addPoint(half_pitch, 0, 0, fine), geo.addPoint(half_pitch, gap, 0, coarse),
               geo.addPoint(half_pitch - 0.5, gap, 0, fine), geo.addPoint(0.5, gap, 0, fine), geo.addPoint(0, gap, 0, coarse)]
        lines = [geo.addLine(pts[i], pts[(i + 1) % 6]) for i in range(6)]
        geo.addPlaneSurface([geo.addCurveLoop(lines)])
        geo.synchronize()
        gmsh.model.mesh.generate(2)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        _, _, conn = gmsh.model.mesh.getElements(2)
    finally:
        gmsh.finalize()
    used = np.unique(conn[0].astype(np.int64))
    remap = np.full(int(tags.max()) + 1, -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    idx = np.zeros(int(tags.max()) + 1, dtype=np.int64)
    idx[tags.astype(np.int64)] = np.arange(len(tags))
    p = coords.reshape(-1, 3)[idx[used], :2].T
    return MeshTri(np.ascontiguousarray(p), np.ascontiguousarray(remap[conn[0].astype(np.int64)].reshape(-1, 3).T))


def run_slot_jet(design: SlotJetDesign, req: Requirements, fluid: Fluid, out_dir: Path | None = None, cells_across: int = 20) -> SlotJetResult:
    b = design.slot_width_m
    half_pitch, gap = design.pitch_m / (2 * b), design.gap_m / b
    n_slots = int((req.footprint_width_m + 1e-12) // design.pitch_m)
    length = req.flow_length_m  # slot length, normal to the solved plane
    v_jet = design.flow_lpm / 60e3 / (n_slots * b * length)
    re = fluid.rho * v_jet * b / fluid.mu
    pr = fluid.cp * fluid.mu / fluid.k

    mesh = _mesh(half_pitch, gap, min(1.0, gap) / cells_across)
    tol = 1e-8
    top = lambda x: x[1] > gap - tol
    inlet = mesh.facets_satisfying(lambda x: top(x) & (x[0] < 0.5), boundaries_only=True)
    outlet = mesh.facets_satisfying(lambda x: top(x) & (x[0] > half_pitch - 0.5), boundaries_only=True)
    cover = mesh.facets_satisfying(lambda x: top(x) & (x[0] > 0.5) & (x[0] < half_pitch - 0.5), boundaries_only=True)
    base = mesh.facets_satisfying(lambda x: x[1] < tol, boundaries_only=True)
    sides = mesh.facets_satisfying(lambda x: (x[0] < tol) | (x[0] > half_pitch - tol), boundaries_only=True)
    jet = lambda x: -1.5 * (1 - (x[0] / 0.5) ** 2)
    steps = [r for r in (5.0, 20.0, 50.0, 100.0, 200.0) if r < re] + [re]
    flow = solve_flow(mesh, steps, {
        0: [(inlet, None), (outlet, None), (sides, None), (cover, None), (base, None)],
        1: [(inlet, jet), (cover, None), (base, None)],
    })
    energy = solve_energy(flow, pr, inlet, base, outlet)

    p_in = flow.p[flow.basis_p.get_dofs(facets=inlet).all()].mean()
    p_out = flow.p[flow.basis_p.get_dofs(facets=outlet).all()].mean()
    dp = (p_in - p_out) * fluid.rho * v_jet**2
    x, theta_w = boundary_profile(energy, base)
    theta_mean = np.trapezoid(theta_w, x) / half_pitch
    area = n_slots * design.pitch_m * length
    flux = req.heat_load_w / area
    scale = flux * b / fluid.k  # K per unit theta
    r_conv = theta_mean * b / (fluid.k * area)
    r_base = req.base_thickness_m / (req.solid_k * area)
    result = SlotJetResult(
        design_id=design.design_id, n_slots=n_slots, reynolds_slot=float(re), jet_velocity_m_s=float(v_jet), dp_pa=float(dp),
        r_th_k_per_w=float(r_conv + r_base), r_conv_k_per_w=float(r_conv), r_base_k_per_w=float(r_base),
        h_stagnation_w_m2k=float(fluid.k / (b * theta_w[0])), h_mean_w_m2k=float(fluid.k / (b * theta_mean)),
        base_avg_temp_c=float(req.inlet_temp_c + req.heat_load_w * (r_conv + r_base)),
        base_max_temp_c=float(req.inlet_temp_c + theta_w.max() * scale + req.heat_load_w * r_base),
        energy_balance_error=energy.energy_balance_error, newton_iterations=flow.newton_iterations,
        n_unknowns=int(flow.basis_u.N + flow.basis_p.N + energy.basis.N),
    )
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "results.json").write_text(json.dumps({"design": asdict(design), "result": asdict(result)}, indent=2))
        psi = stream_function(flow, [(np.concatenate([sides, base]), 0.0), (cover, 0.5)])
        _plot(design, result, fluid, req, flow, energy, psi, x, theta_w, scale, out_dir / "fields.png")
    return result


def _plot(design, result, fluid, req, flow, energy, psi, x, theta_w, scale, out: Path) -> None:
    b_mm = design.slot_width_m * 1e3
    m = flow.mesh
    n = m.p.shape[1]

    def mirrored(ax, values, levels, cmap):
        for sign in (-1, 1):
            cs = ax.tricontourf(Triangulation(sign * m.p[0] * b_mm, m.p[1] * b_mm, m.t.T), values, levels=levels, cmap=cmap)
        return cs

    fig, axes = plt.subplots(3, 1, figsize=(10, 7.6), constrained_layout=True, sharex=True, height_ratios=[1, 1, 1.25])
    speed = flow.speed_at_nodes() * result.jet_velocity_m_s
    cs = mirrored(axes[0], speed, np.linspace(0, speed.max(), 30), "viridis")
    for sign in (-1, 1):
        tri = Triangulation(sign * m.p[0] * b_mm, m.p[1] * b_mm, m.t.T)
        axes[0].tricontour(tri, psi[:n], levels=np.linspace(0.02, 0.48, 9), colors="white", linewidths=0.5)
        axes[0].tricontour(tri, psi[:n], levels=[-0.01, -0.003, -0.0005, 0.5005, 0.503, 0.51], colors="#ff6b6b", linewidths=0.8)
    fig.colorbar(cs, ax=axes[0], label="speed (m/s)")
    axes[0].set_title("Flow: white streamlines, red = recirculation", fontsize=10)
    temp = req.inlet_temp_c + energy.theta[:n] * scale
    cs = mirrored(axes[1], temp, np.linspace(temp.min(), temp.max(), 30), "inferno")
    fig.colorbar(cs, ax=axes[1], label="fluid temperature (°C)")
    axes[1].set_title("Temperature (fluid only; uniform heat flux on the base)", fontsize=10)
    for ax in axes[:2]:
        ax.set_aspect("equal")
        ax.set_ylabel("mm")
    h = fluid.k / (design.slot_width_m * theta_w[1:]) / 1e3
    xs = x[1:] * b_mm
    axes[2].plot(np.concatenate([-xs[::-1], xs]), np.concatenate([h[::-1], h]), color="#c2410c", lw=1.8)
    axes[2].axhline(result.h_mean_w_m2k / 1e3, color="#374151", ls="--", lw=1)
    axes[2].text(0.99, result.h_mean_w_m2k / 1e3, " mean (1/average resistance)", va="bottom", ha="right", fontsize=8,
                 transform=axes[2].get_yaxis_transform())
    axes[2].set_ylabel("local h (kW/m²K)\nreferenced to inlet temperature")
    axes[2].set_xlabel("distance from jet centreline (mm)")
    axes[2].set_ylim(0, None)
    axes[2].grid(alpha=0.25)
    fig.suptitle(
        f"Slot-jet array: {b_mm:.2f} mm slots, {design.gap_m*1e3:.2f} mm gap, {design.pitch_m*1e3:.1f} mm pitch, "
        f"{design.flow_lpm:.1f} LPM {fluid.name}\nRe = {result.reynolds_slot:.0f}, ΔP = {result.dp_pa/1e3:.2f} kPa, "
        f"R_th = {result.r_th_k_per_w*1e3:.1f} mK/W", fontsize=11)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="Solve one slot-jet array design (1 kW, 60x60 mm, PG25 at 40 °C).")
    ap.add_argument("--slot_mm", type=float, required=True)
    ap.add_argument("--gap_mm", type=float, required=True)
    ap.add_argument("--pitch_mm", type=float, required=True)
    ap.add_argument("--flow_lpm", type=float, required=True)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    req = Requirements(1000.0, 0.060, 0.060, 0.002, 390.0, 40.0, 5e3, "PG25", r_th_threshold_k_per_w=20e-3,
                   r_th_objective_k_per_w=15e-3)
    d = SlotJetDesign(a.slot_mm / 1e3, a.gap_mm / 1e3, a.pitch_mm / 1e3, a.flow_lpm)
    r = run_slot_jet(d, req, fluid_props("PG25", 40.0), a.out or Path("runs/level2") / d.design_id)
    print(json.dumps(asdict(r), indent=2))


if __name__ == "__main__":
    main()
