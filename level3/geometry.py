"""Level 3 straight-fin unit cell: bookkeeping on the Level 1 definitions and blockMesh dictionaries for OpenFOAM.

Coordinates (metres): x streamwise (periodic), y across the pitch (y = 0 the fin mid-plane, y = pitch / 2 the channel
mid-plane, both symmetry planes), z up from the chip interface (z = 0, uniform flux) through the 2 mm base to the cover.
    solid: base 0 <= z <= t_base over the whole half pitch, fin 0 <= y <= t/2 above it
    fluid: t/2 <= y <= pitch/2, t_base <= z <= t_base + H (no tip clearance, adiabatic cover)
Two regions (chtMultiRegionSimpleFoam), conformal at the fin side and the channel floor.
"""
import math
from dataclasses import dataclass

from tools.types import Design, Fluid, Requirements

# Cells at level 0; level k multiplies every direction by 2**k. Level 1 is production, level 2 the refinement check.
BASE_CELLS = {"n_channel": 6, "n_base": 6, "n_fluid_z": 24, "n_x": 3}
MESH_LEVELS = {0: 1, 1: 2, 2: 4}
GRADING = {"channel_y": 3.0, "fluid_z": 4.0, "base_z": 0.25}  # last / first cell size within each block (or half block)


@dataclass(frozen=True)
class UnitCell:
    design_id: str
    half_fin_m: float
    half_channel_m: float
    base_m: float
    height_m: float
    flow_length_m: float
    periodic_length_m: float
    n_channels: int
    heat_load_w: float
    flow_m3_s: float
    rho: float
    cp: float
    mu: float
    k_fluid: float
    k_solid: float
    inlet_temp_c: float

    @property
    def half_pitch_m(self) -> float:
        return self.half_fin_m + self.half_channel_m

    @property
    def heated_flux_w_m2(self) -> float:
        """Level 1 spreads the load over n pitches (ASSUMPTIONS M5), not the whole footprint."""
        return self.heat_load_w / (self.n_channels * 2 * self.half_pitch_m * self.flow_length_m)

    @property
    def u_bar(self) -> float:
        return self.flow_m3_s / (self.n_channels * 2 * self.half_channel_m * self.height_m)

    @property
    def coolant_rise_k(self) -> float:
        return self.heat_load_w / (self.rho * self.cp * self.flow_m3_s)

    @property
    def hydraulic_diameter_m(self) -> float:
        w = 2 * self.half_channel_m
        return 2 * w * self.height_m / (w + self.height_m)

    @property
    def reynolds(self) -> float:
        return self.rho * self.u_bar * self.hydraulic_diameter_m / self.mu

    def cell_heat_w(self, length_m: float) -> float:
        """Heat into the modelled half cell over a streamwise length."""
        return self.heated_flux_w_m2 * self.half_pitch_m * length_m

    def cell_flow_m3_s(self) -> float:
        return self.u_bar * self.half_channel_m * self.height_m

    def dp_core_pa(self, dpdx_pa_per_m: float) -> float:
        return dpdx_pa_per_m * self.flow_length_m


def unit_cell(design: Design, req: Requirements, fluid: Fluid, periodic_length_m: float | None = None) -> UnitCell:
    pitch = design.w_ch_m + design.t_fin_m
    return UnitCell(design.design_id, design.t_fin_m / 2, design.w_ch_m / 2, req.base_thickness_m, design.h_fin_m,
                    req.flow_length_m, periodic_length_m or pitch, design.n_channels(req.footprint_width_m), req.heat_load_w,
                    design.flow_lpm / 60e3, fluid.rho, fluid.cp, fluid.mu, fluid.k, req.solid_k, req.inlet_temp_c)


def r_th_from_temperatures(cell: UnitCell, t_base_avg: float, t_bulk: float) -> float:
    """Level 1 definition: (base underside average - local bulk) / Q plus half the coolant rise, i.e. the 60 mm average of
    (base - inlet) under a uniform flux, with fully developed flow."""
    return (t_base_avg - t_bulk + cell.coolant_rise_k / 2) / cell.heat_load_w


def duct_f_re_exact(alpha: float, terms: int = 200) -> float:
    """Fully developed laminar f·Re (Fanning, on Dh) of a rectangular duct, aspect ratio alpha = short / long side; the exact
    series solution (Shah & London 1978, Eq. 341). Used to check the 3D flow solve independently of Level 1."""
    a = min(alpha, 1 / alpha)
    s = sum(math.tanh(n * math.pi / (2 * a)) / n**5 for n in range(1, 2 * terms, 2))
    return 24 / ((1 - 192 * a / math.pi**5 * s) * (1 + a) ** 2)


def grid(cell: UnitCell, factor: int) -> dict:
    n_channel = BASE_CELLS["n_channel"] * factor
    n_fin = max(3, round(BASE_CELLS["n_channel"] * cell.half_fin_m / cell.half_channel_m)) * factor
    n_base, n_fluid_z, n_x = (BASE_CELLS[k] * factor for k in ("n_base", "n_fluid_z", "n_x"))
    lx = cell.periodic_length_m
    return {"factor": factor, "n_fin": n_fin, "n_channel": n_channel, "n_base": n_base, "n_fluid_z": n_fluid_z, "n_x": n_x,
            "lx_m": lx, "dx_m": lx / n_x, "n_cells_fluid": n_x * n_channel * n_fluid_z,
            "n_cells_solid": n_x * (n_fin * n_base + n_channel * n_base + n_fin * n_fluid_z)}


def _header(obj: str, cls: str = "dictionary") -> str:
    return (f"FoamFile\n{{\n    version 2.0;\n    format ascii;\n    class {cls};\n    object {obj};\n}}\n")


def blockmesh_dict(cell: UnitCell, region: str, factor: int) -> str:
    """blockMeshDict for one region. Vertices on the grid x in (0, Lx), y in (0, t/2, pitch/2), z in (0, t_b, t_b + H)."""
    g = grid(cell, factor)
    ys = (0.0, cell.half_fin_m, cell.half_pitch_m)
    zs = (0.0, cell.base_m, cell.base_m + cell.height_m)
    xs = (0.0, g["lx_m"])
    verts, vid = [], {}
    for iz, z in enumerate(zs):
        for iy, y in enumerate(ys):
            for ix, x in enumerate(xs):
                vid[(ix, iy, iz)] = len(verts)
                verts.append(f"    ({x:.9e} {y:.9e} {z:.9e})")

    def hexv(iy, iz):
        c = [(0, iy, iz), (1, iy, iz), (1, iy + 1, iz), (0, iy + 1, iz),
             (0, iy, iz + 1), (1, iy, iz + 1), (1, iy + 1, iz + 1), (0, iy + 1, iz + 1)]
        return " ".join(str(vid[k]) for k in c)

    def face_x(ix, iy, iz):  # streamwise end face of block (iy, iz)
        return f"({vid[(ix, iy, iz)]} {vid[(ix, iy + 1, iz)]} {vid[(ix, iy + 1, iz + 1)]} {vid[(ix, iy, iz + 1)]})"

    def face_y(iy, iz):  # face at y = ys[iy], spanning z from zs[iz] to zs[iz+1]
        return f"({vid[(0, iy, iz)]} {vid[(0, iy, iz + 1)]} {vid[(1, iy, iz + 1)]} {vid[(1, iy, iz)]})"

    def face_z(iy, iz):  # face at z = zs[iz], spanning y from ys[iy] to ys[iy+1]
        return f"({vid[(0, iy, iz)]} {vid[(1, iy, iz)]} {vid[(1, iy + 1, iz)]} {vid[(0, iy + 1, iz)]})"

    gz, gy, gb = GRADING["fluid_z"], GRADING["channel_y"], GRADING["base_z"]
    zgrade = f"((0.5 0.5 {gz:g}) (0.5 0.5 {1 / gz:g}))"
    nx = g["n_x"]
    if region == "fluid":
        blocks = [f"    hex ({hexv(1, 1)}) ({nx} {g['n_channel']} {g['n_fluid_z']}) simpleGrading (1 {gy:g} {zgrade})"]
        ends = [(1, 1)]
        patches = {
            "fluid_to_solid": ("mappedWall", [face_y(1, 1), face_z(1, 1)],
                               "sampleMode nearestPatchFace; sampleRegion solid; samplePatch solid_to_fluid;"),
            "cover": ("wall", [face_z(1, 2)], ""),
            "symFluid": ("symmetryPlane", [face_y(2, 1)], ""),
        }
    elif region == "solid":
        blocks = [f"    hex ({hexv(0, 0)}) ({nx} {g['n_fin']} {g['n_base']}) simpleGrading (1 1 {gb:g})",
                  f"    hex ({hexv(1, 0)}) ({nx} {g['n_channel']} {g['n_base']}) simpleGrading (1 {gy:g} {gb:g})",
                  f"    hex ({hexv(0, 1)}) ({nx} {g['n_fin']} {g['n_fluid_z']}) simpleGrading (1 1 {zgrade})"]
        ends = [(0, 0), (1, 0), (0, 1)]
        patches = {
            "solid_to_fluid": ("mappedWall", [face_y(1, 1), face_z(1, 1)],
                               "sampleMode nearestPatchFace; sampleRegion fluid; samplePatch fluid_to_solid;"),
            "bottom": ("wall", [face_z(0, 0), face_z(1, 0)], ""),
            "finTop": ("wall", [face_z(0, 2)], ""),
            "symFin": ("symmetryPlane", [face_y(0, 0), face_y(0, 1)], ""),  # y = 0
            "symBase": ("symmetryPlane", [face_y(2, 0)], ""),  # y = pitch / 2 under the channel
        }
    else:
        raise ValueError(region)
    patches["front"] = ("cyclic", [face_x(0, *b) for b in ends], "neighbourPatch back;")
    patches["back"] = ("cyclic", [face_x(1, *b) for b in ends], "neighbourPatch front;")
    bnd = []
    for name, (ptype, faces, extra) in patches.items():
        bnd.append(f"    {name}\n    {{\n        type {ptype};\n" + (f"        {extra}\n" if extra else "")
                   + "        faces\n        (\n" + "".join(f"            {f}\n" for f in faces) + "        );\n    }")
    return (_header("blockMeshDict") + "convertToMeters 1;\n\nvertices\n(\n" + "\n".join(verts) + "\n);\n\nblocks\n(\n"
            + "\n".join(blocks) + "\n);\n\nedges ();\n\nboundary\n(\n" + "\n".join(bnd) + "\n);\n\nmergePatchPairs ();\n")
