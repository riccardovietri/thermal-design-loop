"""postprocess: figures from a run directory (inputs.json, results.json, fields.npz, mesh.npz)."""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.tri import Triangulation


def plot_fields(run_dir: Path, out_path: Path, n_channels_shown: int = 3) -> Path:
    """Axial velocity and temperature on the cross-section, mirrored to show whole channels (mid-length plane)."""
    run_dir = Path(run_dir)
    inp = json.loads((run_dir / "inputs.json").read_text())
    res = json.loads((run_dir / "results.json").read_text())
    mesh, fields = np.load(run_dir / "mesh.npz"), np.load(run_dir / "fields.npz")
    p, t, fluid = mesh["p"] * 1e3, mesh["t"].T, mesh["fluid"]
    n = p.shape[1]
    w = fields["u"][:n]
    temp = fields["temp"][:n] + inp["requirements"]["inlet_temp_c"] + res["coolant_rise_k"] / 2
    d = inp["design"]
    pitch = (d["w_ch_m"] + d["t_fin_m"]) * 1e3

    fig, axes = plt.subplots(1, 2, figsize=(9, 6.5), constrained_layout=True)
    panels = [(w, "viridis", "Axial velocity (m/s)", ~fluid), (temp, "inferno", "Temperature at mid-length (°C)", None)]
    for ax, (val, cmap, label, mask) in zip(axes, panels):
        levels = np.linspace(val.min(), val.max(), 25)
        for k in range(n_channels_shown):
            for sign in (-1, 1):
                tri = Triangulation(k * pitch + pitch / 2 + sign * (pitch / 2 - p[0]), p[1], t, mask=mask)
                cs = ax.tricontourf(tri, val, levels=levels, cmap=cmap)
        if mask is not None:
            ax.set_facecolor("#b87333")
        ax.set_aspect("equal")
        ax.set_xlim(0, n_channels_shown * pitch)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y from chip interface (mm)")
        fig.colorbar(cs, ax=ax, label=label, shrink=0.8)
    fig.suptitle(
        f"{d['w_ch_m']*1e3:.2f} mm channels, {d['t_fin_m']*1e3:.2f} mm fins, {d['h_fin_m']*1e3:.1f} mm tall, "
        f"{d['flow_lpm']:.1f} LPM {inp['fluid']['name']}\n"
        f"Re = {res['reynolds']:.0f}, ΔP = {res['dp_pa']/1e3:.2f} kPa, R_th = {res['r_th_k_per_w']*1e3:.1f} mK/W",
        fontsize=11,
    )
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    print(plot_fields(a.run_dir, a.out or a.run_dir / "fields.png"))


if __name__ == "__main__":
    main()
