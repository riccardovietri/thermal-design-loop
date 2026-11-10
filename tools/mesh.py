"""make_geometry + mesh: parametric half-fin / half-channel / base unit cell, meshed with gmsh.

Coordinates (metres): x across the fin pitch, y from the chip interface (y = 0) upward.
    solid base : 0 <= y <= t_base
    solid fin  : 0 <= x <= t_fin/2,      y > t_base
    fluid      : t_fin/2 <= x <= pitch/2, y > t_base
"""
import argparse
import json
from pathlib import Path

import gmsh
import numpy as np

from .types import Design, MeshInfo, to_dict


def mesh(design: Design, base_thickness_m: float, out_dir: Path, cells_across_half_channel: int = 8) -> MeshInfo:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mm = 1e3
    t2, w2 = design.t_fin_m / 2 * mm, design.w_ch_m / 2 * mm
    tb, h = base_thickness_m * mm, design.h_fin_m * mm
    size = min(w2, t2 * 2) / cells_across_half_channel

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        occ = gmsh.model.occ
        base = occ.addRectangle(0, 0, 0, t2 + w2, tb)
        fin = occ.addRectangle(0, tb, 0, t2, h)
        chan = occ.addRectangle(t2, tb, 0, w2, h)
        occ.fragment([(2, base)], [(2, fin), (2, chan)])
        occ.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeMax", size)
        gmsh.option.setNumber("Mesh.MeshSizeMin", size)
        gmsh.model.mesh.generate(2)
        gmsh.write(str(out_dir / "mesh.msh"))
        tags, coords, _ = gmsh.model.mesh.getNodes()
        _, etags, conn = gmsh.model.mesh.getElements(2)
        quality = np.asarray(gmsh.model.mesh.getElementQualities(etags[0], "minSICN"))
    finally:
        gmsh.finalize()

    idx = np.zeros(int(tags.max()) + 1, dtype=np.int64)
    idx[tags.astype(np.int64)] = np.arange(len(tags))
    p = coords.reshape(-1, 3)[:, :2].T / mm
    t = idx[conn[0].astype(np.int64)].reshape(-1, 3).T
    centroid = p[:, t].mean(axis=1)
    fluid = (centroid[0] > design.t_fin_m / 2) & (centroid[1] > base_thickness_m)
    path = out_dir / "mesh.npz"
    np.savez(path, p=p, t=t, fluid=fluid)
    info = MeshInfo(str(path), p.shape[1], t.shape[1], int(fluid.sum()), float(quality.min()), float(quality.mean()), size / mm)
    (out_dir / "mesh.json").write_text(json.dumps(to_dict(info), indent=2))
    return info


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    for name in ("w_ch_mm", "t_fin_mm", "h_fin_mm"):
        ap.add_argument(f"--{name}", type=float, required=True)
    ap.add_argument("--base_mm", type=float, default=2.0)
    ap.add_argument("--cells", type=int, default=8)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    d = Design(a.w_ch_mm / 1e3, a.t_fin_mm / 1e3, a.h_fin_mm / 1e3, flow_lpm=0.0)
    print(json.dumps(to_dict(mesh(d, a.base_mm / 1e3, a.out, a.cells)), indent=2))


if __name__ == "__main__":
    main()
