"""Level 3 stage two: staggered circular pins in 3D (OpenFOAM), bookkeeping, gmsh mesh and case files.

Plan geometry exactly as `level2.planforms.pin_circle` (staggered pins of diameter d, transverse pitch S_T, row pitch S_L,
row B shifted S_T / 2 laterally and S_L streamwise; period 2 S_L). Coordinates in mm inside gmsh, metres in OpenFOAM:
    x streamwise, periodic over [-S_L / 2, 3 S_L / 2] (the cell boundaries pass between pins, so no pin is cut by them);
    y lateral: symmetric cell [0, S_T / 2] (half pin A on y = 0 at x = 0, half pin B on y = S_T / 2 at x = S_L, both
      lateral planes symmetry), or full pitch [-S_T / 2, S_T / 2] (whole pin A, half pins B on lateral periodic planes);
    z up from the chip interface (z = 0, uniform flux) through the base (0 .. t_b) to the cover (t_b + H).
The 2D plan is meshed with quads (gmsh, Frontal-Delaunay for quads, recombined) and extruded in z into hexahedra (and a
few prisms), so the fluid-solid interface is conformal: pin walls and the floor between the pins.

Flow and heat bookkeeping follows the plan view (`level2.planview`): superficial velocity U_s = Q / (W H) over the
60 mm width; heat 1 kW uniform over the whole 60 x 60 mm footprint. The periodic temperature formulation and the sink
are those of stage one (`level3.foam`).
"""
import math
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path

from tools.types import Fluid, Requirements

from .foam import CU_CP, CU_RHO, IMAGE, MOLAR_MASS, T_REF_K, _field, _header, _w

P0_SLAB_MM = 0.05  # P0 cross-check: a thin slab, one cell thick, slip floor and cover (the plan-view problem)
# Plan mesh sizes as fractions of the minimum fluid gap; level k divides all sizes (and multiplies z cells) by 2**k.
PLAN = {"wall": 1 / 32, "core": 1 / 12}
Z_CELLS = {"n_fluid_z": 24, "n_base": 4}  # level 0
GRADING = {"fluid_z": 4.0, "base_z": 0.25}  # last / first in each half of the fin layer; base: top / bottom cell
MESH_LEVELS = {0: 1, 1: 2, 2: 4}


def design_id(params: dict) -> str:
    """Same id as `level2.planview.design_id` for a circular-pin planform."""
    return "pin_circle_" + "_".join(f"{k}{v * 1000:04.0f}" for k, v in params.items())


@dataclass(frozen=True)
class PinCell:
    design_id: str
    d_m: float
    s_t_m: float
    s_l_m: float
    height_m: float
    base_m: float
    footprint_w_m: float
    flow_length_m: float
    heat_load_w: float
    flow_m3_s: float
    rho: float
    cp: float
    mu: float
    k_fluid: float
    k_solid: float
    full_pitch: bool = False

    @property
    def params_mm(self) -> dict:
        return {"d": self.d_m * 1e3, "s_t": self.s_t_m * 1e3, "s_l": self.s_l_m * 1e3, "h_fin": self.height_m * 1e3}

    @property
    def lateral_m(self) -> float:
        return self.s_t_m if self.full_pitch else self.s_t_m / 2

    @property
    def period_m(self) -> float:
        return 2 * self.s_l_m

    @property
    def plan_area_m2(self) -> float:
        return self.lateral_m * self.period_m

    @property
    def pin_area_m2(self) -> float:
        n_pins = 2 if self.full_pitch else 1  # per period: half pin A + half pin B, or a whole A and two half B
        return n_pins * math.pi * self.d_m**2 / 4

    @property
    def porosity(self) -> float:
        return 1 - self.pin_area_m2 / self.plan_area_m2

    @property
    def u_superficial(self) -> float:
        return self.flow_m3_s / (self.footprint_w_m * self.height_m)

    @property
    def cell_flow_m3_s(self) -> float:
        return self.u_superficial * self.height_m * self.lateral_m

    @property
    def heated_flux_w_m2(self) -> float:
        return self.heat_load_w / (self.footprint_w_m * self.flow_length_m)

    @property
    def cell_heat_w(self) -> float:
        return self.heated_flux_w_m2 * self.plan_area_m2

    @property
    def coolant_rise_k(self) -> float:
        return self.heat_load_w / (self.rho * self.cp * self.flow_m3_s)

    @property
    def re_d_max(self) -> float:
        u_max = self.u_superficial * self.s_t_m / (self.s_t_m - self.d_m)
        return self.rho * u_max * self.d_m / self.mu

    def ubar(self, fluid_plan_area_m2: float) -> float:
        """Volume-average streamwise velocity over the fluid (what `meanVelocityForce` holds) for the cell flow:
        Q_cell L / V_fluid = U_s A_plan / A_fluid, with the fluid area of the actual mesh."""
        return self.u_superficial * self.plan_area_m2 / fluid_plan_area_m2

    def dp_core_pa(self, dpdx_pa_per_m: float) -> float:
        return dpdx_pa_per_m * self.flow_length_m


def pin_cell(params: dict, req: Requirements, fluid: Fluid, flow_lpm: float = 1.5, full_pitch: bool = False) -> PinCell:
    p = {k: params[k] for k in ("d", "s_t", "s_l", "h_fin")}
    return PinCell(design_id(p), p["d"] * 1e-3, p["s_t"] * 1e-3, p["s_l"] * 1e-3, p["h_fin"] * 1e-3, req.base_thickness_m,
                   req.footprint_width_m, req.flow_length_m, req.heat_load_w, flow_lpm / 60e3, fluid.rho, fluid.cp,
                   fluid.mu, fluid.k, req.solid_k, full_pitch)


def r_th_pins(cell: PinCell, t_base_avg: float, t_bulk: float) -> float:
    """Stage-one definition: (base underside area average - period-mean bulk) / Q plus half the coolant rise."""
    return (t_base_avg - t_bulk + cell.coolant_rise_k / 2) / cell.heat_load_w


def planview_h(flux_w_m2: float, t_wall: float, t_bulk: float) -> float:
    """Plan-view h: uniform wall flux over (period-mean wall T - period-mean velocity-weighted bulk T)."""
    return flux_w_m2 / (t_wall - t_bulk)


def layer_heights(n: int, ratio: float, symmetric: bool) -> list[float]:
    """Cumulative normalised layer heights for gmsh extrusion. symmetric: geometric growth by `ratio` (last / first) from
    each end to the middle; otherwise one geometric sequence with last / first = ratio."""
    if n == 1:
        return [1.0]
    if symmetric:
        half = n // 2
        g = ratio ** (1 / (half - 1)) if half > 1 else 1.0
        dz = [g**i for i in range(half)]
        dz = dz + ([dz[-1]] if n % 2 else []) + dz[::-1]
    else:
        g = ratio ** (1 / (n - 1))
        dz = [g**i for i in range(n)]
    total = sum(dz)
    out, s = [], 0.0
    for v in dz:
        s += v
        out.append(s / total)
    out[-1] = 1.0
    return out


def pin_centres(params: dict, full_pitch: bool) -> list[tuple]:
    """(x, y, kind) in mm of the pins in the cell: kind 'up' / 'down' = half pin above / below a lateral boundary."""
    s_t, s_l = params["s_t"], params["s_l"]
    if full_pitch:
        return [(0.0, 0.0, "full"), (s_l, -s_t / 2, "up"), (s_l, s_t / 2, "down")]
    return [(0.0, 0.0, "up"), (s_l, s_t / 2, "down")]


def min_gap_mm(params: dict) -> float:
    d, s_t, s_l = params["d"], params["s_t"], params["s_l"]
    return min(s_t - d, math.hypot(s_l, s_t / 2) - d, 2 * s_l - d)


# --------------------------------------------------------------------------------------------------------------------
# gmsh: plan geometry, plan mesh, extrusion


def _plan_geometry(params: dict, full_pitch: bool, with_pins: bool):
    """Build the plan in the geo kernel. Returns dict of tags: fluid surface, pin surfaces, boundary curves."""
    import gmsh
    geo = gmsh.model.geo
    d, s_t, s_l = params["d"], params["s_t"], params["s_l"]
    r = d / 2
    x0, x1 = -s_l / 2, 1.5 * s_l
    y0, y1 = (-s_t / 2, s_t / 2) if full_pitch else (0.0, s_t / 2)
    pins = pin_centres(params, full_pitch)
    pt = {}

    def P(x, y):
        key = (round(x, 9), round(y, 9))
        if key not in pt:
            pt[key] = geo.addPoint(x, y, 0)
        return pt[key]

    arcs, pin_loops = [], []

    def boundary(y, kind):  # segments along y (left to right), split by half pins of `kind` sitting on it
        cs = sorted(p for p in pins if p[2] == kind)
        xs = [x0] + [v for c in cs for v in (c[0] - r, c[0] + r)] + [x1]
        segs, solid = [], []
        for i in range(len(xs) - 1):
            ln = geo.addLine(P(xs[i], y), P(xs[i + 1], y))
            (solid if i % 2 else segs).append(ln)
        return segs, solid, cs

    bot_f, bot_s, bot_pins = boundary(y0, "up")
    top_f, top_s, top_pins = boundary(y1, "down")
    left = geo.addLine(P(x0, y0), P(x0, y1))
    right = geo.addLine(P(x1, y0), P(x1, y1))
    loop = []
    # bottom, left to right: fluid segment, (arc over pin), fluid segment ...
    for i, seg in enumerate(bot_f):
        loop.append(seg)
        if i < len(bot_pins):
            cx = bot_pins[i][0]
            c = P(cx, y0)
            a1 = geo.addCircleArc(P(cx - r, y0), c, P(cx, y0 + r))
            a2 = geo.addCircleArc(P(cx, y0 + r), c, P(cx + r, y0))
            arcs += [a1, a2]
            loop += [a1, a2]
            pin_loops.append([bot_s[i], -a2, -a1])
    loop.append(right)
    # top, right to left
    tarcs = []
    for i, cpin in enumerate(top_pins):
        cx = cpin[0]
        c = P(cx, y1)
        a1 = geo.addCircleArc(P(cx - r, y1), c, P(cx, y1 - r))
        a2 = geo.addCircleArc(P(cx, y1 - r), c, P(cx + r, y1))
        arcs += [a1, a2]
        tarcs.append((a1, a2))
        pin_loops.append([top_s[i], -a2, -a1])
    for i in reversed(range(len(top_f))):
        loop.append(-top_f[i])
        if i > 0:
            a1, a2 = tarcs[i - 1]
            loop += [-a2, -a1]
    loop.append(-left)
    holes = []
    for cx, cy, kind in pins:
        if kind == "full":
            c = P(cx, cy)
            q = [P(cx + r, cy), P(cx, cy + r), P(cx - r, cy), P(cx, cy - r)]
            ring = [geo.addCircleArc(q[k], c, q[(k + 1) % 4]) for k in range(4)]
            arcs += ring
            holes.append(ring)
            pin_loops.append(ring)
    fl = geo.addCurveLoop(loop)
    hl = [geo.addCurveLoop(h) for h in holes]
    fluid = geo.addPlaneSurface([fl] + hl)
    pin_surfs = [geo.addPlaneSurface([geo.addCurveLoop(pl)]) for pl in pin_loops] if with_pins else []
    return {"fluid": fluid, "pins": pin_surfs, "arcs": arcs, "left": left, "right": right, "bot": bot_f + bot_s,
            "top": top_f + top_s, "box": (x0, x1, y0, y1)}


def _size_field(arcs, wall_mm, core_mm, ramp_mm):
    import gmsh
    fd = gmsh.model.mesh.field
    dist = fd.add("Distance")
    fd.setNumbers(dist, "CurvesList", arcs)
    fd.setNumber(dist, "Sampling", 400)
    thr = fd.add("Threshold")
    fd.setNumber(thr, "InField", dist)
    fd.setNumber(thr, "SizeMin", wall_mm)
    fd.setNumber(thr, "SizeMax", core_mm)
    fd.setNumber(thr, "DistMin", 0.0)
    fd.setNumber(thr, "DistMax", ramp_mm)
    fd.setAsBackgroundMesh(thr)
    for opt in ("Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature", "Mesh.MeshSizeExtendFromBoundary"):
        gmsh.option.setNumber(opt, 0)
    gmsh.option.setNumber("Mesh.Algorithm", 8)  # Frontal-Delaunay for quads
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)  # Blossom
    gmsh.option.setNumber("Mesh.RecombineAll", 1)


def _plan_stats(g: dict) -> dict:
    import gmsh
    import numpy as np

    def area(surface):
        tot, n = 0.0, 0
        types, _, conn = gmsh.model.mesh.getElements(2, surface)
        for t, c in zip(types, conn):
            nn = gmsh.model.mesh.getElementProperties(t)[3]
            c = c.reshape(-1, nn)
            for el in c:
                xy = np.array([gmsh.model.mesh.getNode(int(v))[0][:2] for v in el])
                x, y = xy[:, 0], xy[:, 1]
                tot += 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
                n += 1
        return tot, n

    def length(curves):
        tot = 0.0
        for cv in curves:
            types, _, conn = gmsh.model.mesh.getElements(1, abs(cv))
            for t, c in zip(types, conn):
                for a, b in c.reshape(-1, 2):
                    pa, pb = gmsh.model.mesh.getNode(int(a))[0], gmsh.model.mesh.getNode(int(b))[0]
                    tot += math.dist(pa[:2], pb[:2])
        return tot

    def ys(curve):
        tags, coords, _ = gmsh.model.mesh.getNodes(1, curve, includeBoundary=True)
        return sorted(round(v, 9) for v in coords.reshape(-1, 3)[:, 1])

    fa, nf = area(g["fluid"])
    x0, x1, y0, y1 = g["box"]
    fa = float(fa)
    stats = {"fluid_area_mm2": fa, "pin_area_mm2": (x1 - x0) * (y1 - y0) - fa, "pin_perimeter_mm": float(length(g["arcs"])),
             "n_plan_cells": nf, "periodic_nodes_match": ys(g["left"]) == ys(g["right"])}
    types, _, _ = gmsh.model.mesh.getElements(2, g["fluid"])
    stats["plan_element_types"] = [gmsh.model.mesh.getElementProperties(t)[0] for t in types]
    return stats


def _set_periodic(g: dict, full_pitch: bool):
    import gmsh
    x0, x1, y0, y1 = g["box"]
    gmsh.model.mesh.setPeriodic(1, [g["right"]], [g["left"]], [1, 0, 0, x1 - x0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    if full_pitch:
        bot = sorted(g["bot"], key=lambda c: gmsh.model.getBoundingBox(1, c)[0])
        top = sorted(g["top"], key=lambda c: gmsh.model.getBoundingBox(1, c)[0])
        gmsh.model.mesh.setPeriodic(1, top, bot, [1, 0, 0, 0, 0, 1, 0, y1 - y0, 0, 0, 1, 0, 0, 0, 0, 1])


def plan_sizes(params: dict, factor: int) -> tuple[float, float, float]:
    gap = min_gap_mm(params)
    return gap * PLAN["wall"] / factor, gap * PLAN["core"] / factor, gap / 2


def plan_mesh_stats(params: dict, wall_mm: float, core_mm: float, full_pitch: bool = False) -> dict:
    """Mesh the plan only and measure it (used by tests and to size cases)."""
    import gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("plan")
        g = _plan_geometry(params, full_pitch, with_pins=True)
        gmsh.model.geo.synchronize()
        _set_periodic(g, full_pitch)
        _size_field(g["arcs"], wall_mm, core_mm, min_gap_mm(params) / 2)
        gmsh.model.mesh.generate(2)
        return _plan_stats(g)
    finally:
        gmsh.finalize()


def write_mesh(path: Path, params: dict, factor: int, mode: str, full_pitch: bool = False) -> dict:
    """Plan mesh + extrusion, written as gmsh 2.2 (mm). mode 'conjugate': base + pins + fluid, cellZones fluid / solid;
    mode 'p0': the fluid layer only, one cell over a P0_SLAB_MM slab. Physical surface names end in 'Raw'; createPatch
    turns them into the final patches."""
    import gmsh
    wall, core, ramp = plan_sizes(params, factor)
    h_mm, tb_mm = params["h_fin"], 2.0
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("cell")
        geo = gmsh.model.geo
        g = _plan_geometry(params, full_pitch, with_pins=(mode == "conjugate"))
        plan = [g["fluid"]] + g["pins"]
        if mode == "conjugate":
            nb, nz = Z_CELLS["n_base"] * factor, Z_CELLS["n_fluid_z"] * factor
            hb = layer_heights(nb, GRADING["base_z"], symmetric=False)
            hz = layer_heights(nz, GRADING["fluid_z"], symmetric=True)
            ext_b = geo.extrude([(2, s) for s in plan], 0, 0, tb_mm, [1] * nb, hb, recombine=True)
            # extrude returns, per input surface: top surface, volume, lateral surfaces
            idx = [i for i, e in enumerate(ext_b) if e[0] == 3]
            vols_b = [ext_b[i][1] for i in idx]
            tops = [ext_b[i - 1][1] for i in idx]
            ext_f = geo.extrude([(2, s) for s in tops], 0, 0, h_mm, [1] * nz, hz, recombine=True)
            vols_f = [e[1] for e in ext_f if e[0] == 3]
            fluid_vols, solid_vols = [vols_f[0]], vols_b + vols_f[1:]
            z_top = tb_mm + h_mm
        else:
            ext = geo.extrude([(2, g["fluid"])], 0, 0, P0_SLAB_MM, [1], [1.0], recombine=True)
            fluid_vols, solid_vols = [e[1] for e in ext if e[0] == 3], []
            z_top = P0_SLAB_MM
        geo.synchronize()
        _set_periodic(g, full_pitch)
        _size_field(g["arcs"], wall, core, ramp)
        gmsh.model.mesh.generate(2)
        stats = _plan_stats(g)
        gmsh.model.mesh.generate(3)
        x0, x1, y0, y1 = g["box"]
        eps = 1e-7
        groups = {}
        for _, s in gmsh.model.getEntities(2):
            bb = gmsh.model.getBoundingBox(2, s)
            up = gmsh.model.getAdjacencies(2, s)[0]
            if abs(bb[5] - bb[2]) < eps:
                z = bb[2]
                if abs(z) < eps:
                    name = "bottomRaw" if mode == "conjugate" else "floorRaw"
                elif abs(z - z_top) < eps:
                    name = "coverRaw" if any(v in fluid_vols for v in up) else "pinTopRaw"
                else:
                    continue  # internal (floor between base and fin layer)
            elif abs(bb[3] - bb[0]) < eps and abs(bb[0] - x0) < eps:
                name = "frontRaw"
            elif abs(bb[3] - bb[0]) < eps and abs(bb[0] - x1) < eps:
                name = "backRaw"
            elif abs(bb[4] - bb[1]) < eps and abs(bb[1] - y0) < eps:
                name = "loRaw"
            elif abs(bb[4] - bb[1]) < eps and abs(bb[1] - y1) < eps:
                name = "hiRaw"
            elif mode == "p0":
                name = "pinsRaw"
            else:
                continue  # pin wall: interface, created by splitMeshRegions
            groups.setdefault(name, []).append(s)
        for name, ss in groups.items():
            gmsh.model.addPhysicalGroup(2, ss, name=name)
        gmsh.model.addPhysicalGroup(3, fluid_vols, name="fluid")
        if solid_vols:
            gmsh.model.addPhysicalGroup(3, solid_vols, name="solid")
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.SaveAll", 0)
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(path))
        types = {}
        for dim, v in gmsh.model.getEntities(3):
            for t, tags, _ in zip(*gmsh.model.mesh.getElements(3, v)):
                nm = gmsh.model.mesh.getElementProperties(t)[0]
                types[nm] = types.get(nm, 0) + len(tags)
        stats.update({"factor": factor, "wall_mm": wall, "core_mm": core, "mode": mode, "full_pitch": full_pitch,
                      "patches": sorted(groups), "element_counts": types,
                      "n_cells": sum(types.values())})
        if mode == "conjugate":
            stats.update({"n_base_layers": nb, "n_fluid_layers": nz})
        return stats
    finally:
        gmsh.finalize()


# --------------------------------------------------------------------------------------------------------------------
# OpenFOAM case


def _create_patch_dict(entries: list[tuple]) -> str:
    out = []
    for name, info, raw in entries:
        out.append(f"    {{\n        name {name};\n        patchInfo {{ {info} }}\n        constructFrom patches;\n"
                   f"        patches ({raw});\n    }}")
    return _header("createPatchDict") + "pointSync false;\n\npatches\n(\n" + "\n".join(out) + "\n);\n"


def _cyclic_pair(a, b, raw_a, raw_b, sep):
    v = lambda s: "(" + " ".join(f"{float(c):.9g}" for c in s) + ")"
    return [(a, f"type cyclic; neighbourPatch {b}; transform translational; separationVector {v(sep)}; matchTolerance 0.001;", raw_a),
            (b, f"type cyclic; neighbourPatch {a}; transform translational; separationVector {v([-x for x in sep])}; matchTolerance 0.001;", raw_b)]


def patch_layout(cell: PinCell, mode: str) -> dict:
    """Final patches per region: name -> (createPatch info, raw patch, kind)."""
    lx, ly = cell.period_m, cell.lateral_m
    sym = "type symmetryPlane;"
    lat = (_cyclic_pair("sideLo", "sideHi", "loRaw", "hiRaw", [0.0, ly, 0.0]) if cell.full_pitch
           else [("symLo", sym, "loRaw"), ("symHi", sym, "hiRaw")])
    stream = _cyclic_pair("front", "back", "frontRaw", "backRaw", [lx, 0.0, 0.0])
    if mode == "p0":
        return {"fluid": stream + lat + [("floor", sym, "floorRaw"), ("cover", sym, "coverRaw"), ("pins", "type wall;", "pinsRaw")]}
    return {"fluid": stream + lat + [("cover", "type wall;", "coverRaw")],
            "solid": stream + lat + [("bottom", "type wall;", "bottomRaw"), ("pinTop", "type wall;", "pinTopRaw")]}


def _sink(q_cell: float, lx: float, extra: str = "") -> str:
    return f"""
periodicSink
{{
    type            scalarCodedSource;
    selectionMode   all;
    fields          (h);
    name            periodicSink;
    codeInclude
    #{{
    #}};
    codeCorrect
    #{{
    #}};
    codeConstrain
    #{{
    #}};
    codeAddSup
    #{{
        FatalErrorInFunction << "compressible solver expected" << exit(FatalError);
    #}};
    codeAddSupRho
    #{{
        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
        const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
        const scalarField& V = mesh_.V();
        const scalarField ux(U.primitiveField().component(0));
        const scalar sumUV = gSum(ux*V);
        const scalar Qcell = {q_cell!r};
        eqn.source() += Qcell*ux*V/sumUV;
        Info<< "L3MON bulkT " << gSum(ux*T.primitiveField()*V)/sumUV << " flow " << sumUV/{lx!r}
            << " sink " << Qcell << " maxUy " << gMax(mag(U.primitiveField().component(1))) << endl;
    #}};
}}
"""


def write_pin_case(case: Path, cell: PinCell, mesh: dict, mode: str, n_iter: int, transient: dict | None = None,
                   n_procs: int = 1) -> dict:
    """Case files beside an existing mesh.msh. mode 'conjugate' (two regions) or 'p0' (fluid only, flux on pin walls).
    transient: {'end_time': s, 'max_co': c, 'write_interval': s} runs chtMultiRegionFoam from the steady result."""
    case = Path(case)
    lx = cell.period_m
    fluid_area = float(mesh["fluid_area_mm2"]) * 1e-6
    ubar = float(cell.ubar(fluid_area))
    if mode == "p0":
        slab = P0_SLAB_MM * 1e-3
        wall_area = float(mesh["pin_perimeter_mm"]) * 1e-3 * slab
        q_wall = 5.0e4  # any value: the P0 problem is linear in the flux
        q_cell = q_wall * wall_area
    else:
        q_wall, q_cell = None, cell.cell_heat_w
    layout = patch_layout(cell, mode)
    regions = list(layout)
    for r, entries in layout.items():
        _w(case / "system" / r / "createPatchDict", _create_patch_dict([(n, i, raw) for n, i, raw in entries]))
    solid_list = "solid (solid)" if "solid" in regions else "solid ()"
    _w(case / "constant/regionProperties", _header("regionProperties") + f"regions (fluid (fluid) {solid_list});\n")
    pr = cell.cp * cell.mu / cell.k_fluid
    _w(case / "constant/fluid/thermophysicalProperties", _header("thermophysicalProperties") + f"""
thermoType {{ type heRhoThermo; mixture pureMixture; transport const; thermo hConst; equationOfState rhoConst;
             specie specie; energy sensibleEnthalpy; }}
mixture
{{
    specie {{ molWeight {MOLAR_MASS}; }}
    equationOfState {{ rho {cell.rho!r}; }}
    thermodynamics {{ Cp {cell.cp!r}; Hf 0; }}
    transport {{ mu {cell.mu!r}; Pr {pr!r}; }}
}}
""")
    _w(case / "constant/fluid/turbulenceProperties", _header("turbulenceProperties") + "simulationType laminar;\n")
    _w(case / "constant/fluid/radiationProperties", _header("radiationProperties") + "radiation off;\nradiationModel none;\n")
    _w(case / "constant/g", _header("g", "uniformDimensionedVectorField") + "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);\n")
    _w(case / "constant/fluid/fvOptions", _header("fvOptions") + f"""
meanFlow
{{
    type            meanVelocityForce;
    selectionMode   all;
    fields          (U);
    Ubar            ({ubar!r} 0 0);
}}
""" + _sink(q_cell, lx))
    noslip = "        type noSlip;\n"
    cyc, sym = "        type cyclic;\n", "        type symmetryPlane;\n"
    coupled = ("        type compressible::turbulentTemperatureCoupledBaffleMixed;\n"
               "        Tnbr T;\n        kappaMethod {km};\n        value uniform " + f"{T_REF_K};\n")
    pfix = {"p_rgh": "        type fixedFluxPressure;\n        value uniform 1e5;\n",
            "p": "        type calculated;\n        value uniform 1e5;\n"}

    def kind(info):
        return cyc if "cyclic" in info else (sym if "symmetryPlane" in info else None)

    def fluid_bc(field):
        out = {}
        for name, info, _ in layout["fluid"]:
            k = kind(info)
            if k:
                out[name] = k
            elif name == "pins":  # P0 heated pin wall
                out[name] = {"U": noslip, "T": ("        type externalWallHeatFluxTemperature;\n        mode flux;\n"
                                                f"        q uniform {q_wall!r};\n        kappaMethod fluidThermo;\n"
                                                f"        value uniform {T_REF_K};\n")}.get(field, pfix.get(field))
            else:  # cover (wall)
                out[name] = {"U": noslip, "T": "        type zeroGradient;\n"}.get(field, pfix.get(field))
        if "solid" in layout:
            out["fluid_to_solid"] = {"U": noslip, "T": coupled.format(km="fluidThermo")}.get(field, pfix.get(field))
        return out

    z0 = "0.orig"
    u0 = f"({ubar!r} 0 0)"
    _w(case / z0 / "fluid/U", _field("U", "volVectorField", "[0 1 -1 0 0 0 0]", u0, fluid_bc("U")))
    _w(case / z0 / "fluid/T", _field("T", "volScalarField", "[0 0 0 1 0 0 0]", f"{T_REF_K}", fluid_bc("T")))
    _w(case / z0 / "fluid/p_rgh", _field("p_rgh", "volScalarField", "[1 -1 -2 0 0 0 0]", "1e5", fluid_bc("p_rgh")))
    _w(case / z0 / "fluid/p", _field("p", "volScalarField", "[1 -1 -2 0 0 0 0]", "1e5", fluid_bc("p")))
    if "solid" in layout:
        _w(case / "constant/solid/thermophysicalProperties", _header("thermophysicalProperties") + f"""
thermoType {{ type heSolidThermo; mixture pureMixture; transport constIso; thermo hConst; equationOfState rhoConst;
             specie specie; energy sensibleEnthalpy; }}
mixture
{{
    specie {{ molWeight 63.5; }}
    transport {{ kappa {cell.k_solid!r}; }}
    thermodynamics {{ Hf 0; Cp {CU_CP}; }}
    equationOfState {{ rho {CU_RHO}; }}
}}
""")
        _w(case / "constant/solid/radiationProperties", _header("radiationProperties") + "radiation off;\nradiationModel none;\n")
        flux = ("        type externalWallHeatFluxTemperature;\n        mode flux;\n"
                f"        q uniform {cell.heated_flux_w_m2!r};\n        kappaMethod solidThermo;\n        value uniform {T_REF_K};\n")
        spat = {}
        for name, info, _ in layout["solid"]:
            spat[name] = kind(info) or (flux if name == "bottom" else "        type zeroGradient;\n")
        spat["solid_to_fluid"] = coupled.format(km="solidThermo")
        _w(case / z0 / "solid/T", _field("T", "volScalarField", "[0 0 0 1 0 0 0]", f"{T_REF_K}", spat))
        _w(case / z0 / "solid/p", _field("p", "volScalarField", "[1 -1 -2 0 0 0 0]", "1e5",
                                         {k: (v if v in (cyc, sym) else "        type calculated;\n        value uniform 1e5;\n")
                                          for k, v in spat.items()}))
    # --- numerics ---
    steady = transient is None
    ddt = "steadyState" if steady else "backward"
    _w(case / "system/fluid/fvSchemes", _header("fvSchemes") + f"""
ddtSchemes {{ default {ddt}; }}
gradSchemes {{ default Gauss linear; }}
divSchemes
{{
    default none;
    div(phi,U) bounded Gauss linearUpwind grad(U);
    div(phi,h) bounded Gauss linearUpwind grad(h);
    div(phi,K) bounded Gauss linear;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}}
laplacianSchemes {{ default Gauss linear limited corrected 0.5; }}
interpolationSchemes {{ default linear; }}
snGradSchemes {{ default limited corrected 0.5; }}
""")
    if steady:
        algo = """SIMPLE
{
    momentumPredictor yes;
    nNonOrthogonalCorrectors 0;
    pRefCell 0;
    pRefValue 1e5;
}
relaxationFactors
{
    fields { rho 1; p_rgh 0.3; }
    equations { U 0.7; h 0.9; }
}
"""
    else:
        algo = """PIMPLE
{
    momentumPredictor yes;
    nOuterCorrectors 1;
    nCorrectors 2;
    nNonOrthogonalCorrectors 0;
    pRefCell 0;
    pRefValue 1e5;
}
"""
    final = "" if steady else '    "(U|h)Final" { solver PBiCGStab; preconditioner DILU; tolerance 1e-12; relTol 0; }\n' \
        '    p_rghFinal { solver GAMG; smoother GaussSeidel; tolerance 1e-10; relTol 0; }\n' \
        '    rhoFinal { solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0; }\n'
    _w(case / "system/fluid/fvSolution", _header("fvSolution") + f"""
solvers
{{
    rho {{ solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0; }}
    p_rgh {{ solver GAMG; smoother GaussSeidel; tolerance 1e-12; relTol 0.01; maxIter 50; }}
    "(U|h)" {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-14; relTol 0.01; }}
{final}}}
{algo}""")
    if "solid" in layout:
        _w(case / "system/solid/fvSchemes", _header("fvSchemes") + f"""
ddtSchemes {{ default {ddt}; }}
gradSchemes {{ default Gauss linear; }}
divSchemes {{ default none; }}
laplacianSchemes {{ default Gauss linear limited corrected 0.5; }}
interpolationSchemes {{ default linear; }}
snGradSchemes {{ default limited corrected 0.5; }}
""")
        sfinal = "" if steady else ' hFinal { solver GAMG; smoother GaussSeidel; tolerance 1e-12; relTol 0; }'
        salgo = "SIMPLE { nNonOrthogonalCorrectors 0; }\nrelaxationFactors { equations { h 1; } }\n" if steady else \
            "PIMPLE { nNonOrthogonalCorrectors 0; }\n"
        _w(case / "system/solid/fvSolution", _header("fvSolution") + f"""
solvers {{ h {{ solver GAMG; smoother DICGaussSeidel; tolerance 1e-14; relTol 0.1; maxIter 20; }}{sfinal} }}
{salgo}""")
    # top level: read by the mesh utilities (splitMeshRegions) only
    _w(case / "system/fvSchemes", _header("fvSchemes") + "ddtSchemes { default steadyState; }\ngradSchemes { default Gauss linear; }\n"
       "divSchemes { default none; }\nlaplacianSchemes { default Gauss linear corrected; }\n"
       "interpolationSchemes { default linear; }\nsnGradSchemes { default corrected; }\n")
    _w(case / "system/fvSolution", _header("fvSolution") + ("" if steady else "PIMPLE { nOuterCorrectors 1; }\n"))
    if n_procs > 1:
        dec = _header("decomposeParDict") + f"numberOfSubdomains {n_procs};\nmethod scotch;\n"
        _w(case / "system/decomposeParDict", dec)
        for r in regions:
            _w(case / "system" / r / "decomposeParDict", dec)
    zmid = (cell.base_m + cell.height_m / 2) if mode == "conjugate" else P0_SLAB_MM * 1e-3 / 2
    mon = []
    if mode == "conjugate":
        mon.append("""    baseT
    {
        type surfaceFieldValue; libs (fieldFunctionObjects); region solid; regionType patch; name bottom;
        operation areaAverage; fields (T); writeFields false; writeControl timeStep; writeInterval 1; log false;
    }
    residualsSolid
    {
        type solverInfo; libs (utilityFunctionObjects); region solid; fields (h); writeResidualFields no;
    }
    heatFluid
    {
        type wallHeatFlux; libs (fieldFunctionObjects); region fluid; patches (fluid_to_solid);
        writeControl writeTime; writeFields false;
    }
    heatSolid
    {
        type wallHeatFlux; libs (fieldFunctionObjects); region solid; patches (bottom solid_to_fluid);
        writeControl writeTime; writeFields false;
    }
    floorT
    {
        type surfaceFieldValue; libs (fieldFunctionObjects); region fluid; regionType patch; name fluid_to_solid;
        operation areaAverage; fields (T); writeFields false; writeControl writeTime; log false;
    }
    planes
    {
        type surfaces; libs (sampling); region fluid; surfaceFormat raw; fields (U T); interpolationScheme cell;
        writeControl writeTime;
        surfaces { midZ { type cuttingPlane; point (0 0 %r); normal (0 0 1); interpolate false; } }
    }
    solidPlanes
    {
        type surfaces; libs (sampling); region solid; surfaceFormat raw; fields (T); interpolationScheme cell;
        writeControl writeTime;
        surfaces { lo { type patch; patches (%s); interpolate false; } wall { type patch; patches (solid_to_fluid); interpolate false; } }
    }
""" % (zmid, "symLo" if not cell.full_pitch else "sideLo"))
    else:
        mon.append("""    wallT
    {
        type surfaceFieldValue; libs (fieldFunctionObjects); region fluid; regionType patch; name pins;
        operation areaAverage; fields (T); writeFields false; writeControl timeStep; writeInterval 1; log false;
    }
    heatFluid
    {
        type wallHeatFlux; libs (fieldFunctionObjects); region fluid; patches (pins);
        writeControl writeTime; writeFields false;
    }
""")
    if steady:
        time_ctl = f"""application     chtMultiRegionSimpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {n_iter};
deltaT          1;
writeControl    timeStep;
writeInterval   {n_iter};
purgeWrite      1;
"""
    else:
        time_ctl = f"""application     chtMultiRegionFoam;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {transient['end_time']!r};
deltaT          {transient['dt0']!r};
adjustTimeStep  yes;
maxCo           {transient['max_co']!r};
maxDi           1e6;  // solid is implicit: the time step is set by the fluid Courant number
writeControl    adjustable;
writeInterval   {transient['write_interval']!r};
purgeWrite      3;
"""
    _w(case / "system/controlDict", _header("controlDict") + time_ctl + f"""writeFormat     ascii;
writePrecision  12;
timePrecision   10;
runTimeModifiable false;

functions
{{
    residualsFluid
    {{
        type solverInfo; libs (utilityFunctionObjects); region fluid; fields (U p_rgh h); writeResidualFields no;
    }}
{''.join(mon)}}}
""")
    return {"ubar_m_s": ubar, "q_cell_w": q_cell, "q_wall_w_m2": q_wall, "fluid_area_m2": fluid_area, "regions": regions}


def allrun(regions: list[str], n_procs: int = 1, application: str = "chtMultiRegionSimpleFoam", mesh: bool = True) -> str:
    """Shell script run inside the container."""
    lines = ["set -e", "cd /case"]
    if mesh:
        lines += ["gmshToFoam mesh.msh > log.gmshToFoam 2>&1",
                  "transformPoints -scale '(0.001 0.001 0.001)' > log.transformPoints 2>&1"]
        if "solid" in regions:
            lines += ["splitMeshRegions -cellZonesOnly -overwrite > log.splitMeshRegions 2>&1", "rm -rf constant/polyMesh"]
        else:
            lines += ["mkdir -p constant/fluid && rm -rf constant/fluid/polyMesh && mv constant/polyMesh constant/fluid/polyMesh"]
        for r in regions:
            lines.append(f"createPatch -region {r} -overwrite > log.createPatch.{r} 2>&1")
            lines.append(f"checkMesh -region {r} > log.checkMesh.{r} 2>&1")
        lines += ["rm -rf 0 && cp -r 0.orig 0"]
    if n_procs > 1:
        lines += ["rm -rf processor*", "decomposePar -allRegions > log.decomposePar 2>&1",
                  f"mpirun --oversubscribe -np {n_procs} {application} -parallel > log.solver 2>&1 || true",
                  "reconstructPar -allRegions -latestTime > log.reconstructPar 2>&1 || true"]
    else:
        lines += [f"{application} > log.solver 2>&1 || true"]
    lines += ["cat /sys/fs/cgroup/memory.peak > memory.peak 2>/dev/null || true"]
    return "\n".join(lines) + "\n"


def docker_cmd(case: Path, cpus: int = 1, detach: bool = False, name: str | None = None,
               source: Path | None = None) -> list[str]:
    """`source`: another case mounted read-only at /src (to map a converged solution onto a finer mesh)."""
    case = Path(case).resolve()
    cmd = ["docker", "run"] + (["-d"] if detach else ["--rm"]) + (["--name", name] if name else [])
    src = ["-v", f"{Path(source).resolve()}:/src:ro"] if source else []
    return cmd + ["--cpus", str(cpus), "--shm-size", "1g", "--user", f"{os.getuid()}:{os.getgid()}", "-e", "HOME=/tmp",
                  "-v", f"{case}:/case"] + src + ["-w", "/case", IMAGE, "-c", "bash /case/Allrun"]


def run_docker(case: Path, cpus: int = 1) -> None:
    subprocess.run(docker_cmd(case, cpus), check=True)
