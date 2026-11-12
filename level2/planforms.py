"""Plan-view fin geometries for the Level 2 plan-view model (all lengths in mm).

The plan view is the plane parallel to the base; fins are extruded shapes of infinite height in the 2D solve.
Each builder returns a `Planform`: the geometric quantities the resistance network needs, the manufacturing checks,
and a gmsh description of the fluid region so `mesh_planform` can mesh either
  * a streamwise-periodic cell of one or more periods (production: developed-periodic values), or
  * a developing strip: uniform-inflow lead-in, N periods, outlet buffer (cross-check of periodic development).

Topologies and solved domains (x streamwise, y lateral):
  straight     half channel, fin wall at y = 0, symmetry at y = w/2; lateral unit p/2, p = w + t.
  sinusoidal   one full channel between two fin walls, centreline y = A cos(2 pi x / lam) (a crest at each period
               boundary), walls offset +-w/2 along the local normal; lateral pitch p chosen so the thinnest point of
               the fin (normal thickness) is exactly t.
  zigzag       as sinusoidal with a triangular-wave centreline of amplitude A (sharp-corner value); corners rounded with
               an inner fluid-wall radius r_in = 0.05 mm (centreline radius r_in + w/2); p = (w + t) / cos(phi).
  offset_strip strip y in [0, p_t/2], symmetry on both sides; rectangular strips l x t, row A centred on y = 0,
               row B on y = p_t/2 shifted by l + g; period 2 (l + g).
  pin_circle   strip y in [0, S_T/2], symmetry on both sides; staggered circular pins of diameter d, row A on y = 0,
               row B on y = S_T/2 shifted by S_L; period 2 S_L.
  pin_diamond  as pin_circle with squares of side s rotated 45 deg (vertices on the x and y axes).
"""
import math
from dataclasses import dataclass, field

import gmsh
import numpy as np
from skfem import MeshTri

MIN_SOLID_MM = 0.15
MIN_GAP_MM = 0.15
MAX_FIN_ASPECT = 10.0
MAX_FIN_HEIGHT_MM = 4.0
ZIGZAG_INNER_RADIUS_MM = 0.05


class GeometryRejected(ValueError):
    """The design breaks a manufacturing limit or cannot be built; `reasons` lists why."""

    def __init__(self, reasons: list[str]):
        super().__init__("; ".join(reasons))
        self.reasons = reasons


@dataclass(frozen=True)
class Planform:
    topology: str
    params: dict
    kind: str  # "half_channel" | "channel" | "strip"
    period_mm: float
    lateral_unit_mm: float  # plan width (y) that the solved domain represents
    wetted_perimeter_mm: float  # fin wall length in the solved domain, per period
    fluid_area_mm2: float  # fluid plan area in the solved domain, per period
    fin_section_perimeter_mm: float  # P_c of one fin element (per period for continuous fins)
    fin_section_area_mm2: float  # A_c of one fin element
    min_solid_mm: float
    min_gap_mm: float
    min_fin_thickness_mm: float
    ref_gap_mm: float  # gap used to size the mesh
    violations: tuple = ()
    notes: tuple = ()
    geo: dict = field(default_factory=dict, repr=False, compare=False)

    @property
    def perimeter_density_per_mm(self) -> float:
        """Wetted fin perimeter per unit plan area (1/mm)."""
        return self.wetted_perimeter_mm / (self.period_mm * self.lateral_unit_mm)

    @property
    def porosity(self) -> float:
        return self.fluid_area_mm2 / (self.period_mm * self.lateral_unit_mm)

    @property
    def dh_mm(self) -> float:
        """Plan-view hydraulic diameter 4 A_fluid / P_wetted (2w for straight fins)."""
        return 4 * self.fluid_area_mm2 / self.wetted_perimeter_mm

    def summary(self) -> dict:
        return {"topology": self.topology, "params_mm": self.params, "kind": self.kind, "period_mm": self.period_mm,
                "lateral_unit_mm": self.lateral_unit_mm, "wetted_perimeter_mm": self.wetted_perimeter_mm,
                "fluid_area_mm2": self.fluid_area_mm2, "fin_section_perimeter_mm": self.fin_section_perimeter_mm,
                "fin_section_area_mm2": self.fin_section_area_mm2, "perimeter_density_per_mm": self.perimeter_density_per_mm,
                "porosity": self.porosity, "dh_mm": self.dh_mm, "min_solid_mm": self.min_solid_mm,
                "min_gap_mm": self.min_gap_mm, "min_fin_thickness_mm": self.min_fin_thickness_mm,
                "violations": list(self.violations), "notes": list(self.notes)}


def _limits(h_fin: float, min_solid: float, min_gap: float, min_thick: float) -> list[str]:
    v = []
    if min_solid < MIN_SOLID_MM - 1e-9:
        v.append(f"min_solid_feature {min_solid:.3f} mm < {MIN_SOLID_MM} mm")
    if min_gap < MIN_GAP_MM - 1e-9:
        v.append(f"min_fluid_gap {min_gap:.3f} mm < {MIN_GAP_MM} mm")
    if h_fin > MAX_FIN_HEIGHT_MM + 1e-9:
        v.append(f"fin_height {h_fin:.2f} mm > {MAX_FIN_HEIGHT_MM} mm")
    if h_fin / min_thick > MAX_FIN_ASPECT + 1e-9:
        v.append(f"fin_aspect H/t_min = {h_fin / min_thick:.2f} > {MAX_FIN_ASPECT:g}")
    return v


def _check(pf: Planform, strict: bool) -> Planform:
    if strict and pf.violations:
        raise GeometryRejected(list(pf.violations))
    return pf


def _channel_aspect_note(w: float, h_fin: float) -> tuple:
    return (f"channel aspect H/w = {h_fin / w:.1f} > 10 (Level 1 skiving limit, not applied here)",) if h_fin / w > 10 + 1e-9 else ()


# --------------------------------------------------------------------------------------------------------------------
# Polygon helpers (exact distances between polygons; circles as fine polygons)


def _poly_dist(a: np.ndarray, b: np.ndarray) -> float:
    """Minimum distance between two closed polygons (n x 2), assuming they do not intersect."""

    def pt_seg(p, s0, s1):
        d = s1 - s0
        t = np.clip(((p[:, None, :] - s0[None]) * d[None]).sum(-1) / (d * d).sum(-1)[None], 0, 1)
        q = s0[None] + t[..., None] * d[None]
        return np.sqrt(((p[:, None, :] - q) ** 2).sum(-1)).min()

    return float(min(pt_seg(a, b, np.roll(b, -1, 0)), pt_seg(b, a, np.roll(a, -1, 0))))


def _circle(c, r, n=720):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.column_stack([c[0] + r * np.cos(t), c[1] + r * np.sin(t)])


def _diamond(c, a):
    return np.array([[c[0] + a, c[1]], [c[0], c[1] + a], [c[0] - a, c[1]], [c[0], c[1] - a]])


def _rect(x0, x1, y0, y1):
    return np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])


# --------------------------------------------------------------------------------------------------------------------
# Builders


def straight(w: float, t: float, h_fin: float, strict: bool = True) -> Planform:
    """Straight continuous fins: channel width w, fin thickness t. Period length is arbitrary (set to w)."""
    v = _limits(h_fin, t, w, t)
    pf = Planform("straight", {"w": w, "t": t, "h_fin": h_fin}, "half_channel", period_mm=w, lateral_unit_mm=(w + t) / 2,
                  wetted_perimeter_mm=w, fluid_area_mm2=w * w / 2, fin_section_perimeter_mm=2 * w, fin_section_area_mm2=t * w,
                  min_solid_mm=t, min_gap_mm=w, min_fin_thickness_mm=t, ref_gap_mm=w, violations=tuple(v),
                  notes=_channel_aspect_note(w, h_fin), geo={"half_width": w / 2})
    return _check(pf, strict)


def _sin_walls(w, amp, lam, n=400):
    x = np.linspace(0, lam, n + 1)
    k = 2 * np.pi / lam
    yc, dy = amp * np.cos(k * x), -amp * k * np.sin(k * x)
    s = np.sqrt(1 + dy**2)
    nx, ny = -dy / s, 1 / s
    lower = np.column_stack([x - 0.5 * w * nx, yc - 0.5 * w * ny])
    upper = np.column_stack([x + 0.5 * w * nx, yc + 0.5 * w * ny])
    lower[[0, -1], 0], upper[[0, -1], 0] = (0, lam), (0, lam)  # exact on the period boundaries (crests)
    return lower, upper


def _polyline_length(p):
    return float(np.sum(np.hypot(*np.diff(p, axis=0).T)))


def sinusoidal(w: float, t: float, amp: float, lam: float, h_fin: float, strict: bool = True) -> Planform:
    """Sinusoidal channel of normal gap w; fin of minimum normal thickness t."""
    k = 2 * np.pi / lam
    r_min = 1 / (amp * k * k)  # centreline radius of curvature at the crests
    v = []
    if r_min - w / 2 < ZIGZAG_INNER_RADIUS_MM - 1e-9:
        v.append(f"inner wall radius {r_min - w / 2:.3f} mm < {ZIGZAG_INNER_RADIUS_MM} mm (offset wall would fold)")
        pf = Planform("sinusoidal", {"w": w, "t": t, "amp": amp, "lam": lam, "h_fin": h_fin}, "channel", lam, w + t, 0, 0, 0, 0,
                      t, w, t, w, violations=tuple(v + _limits(h_fin, t, w, t)))
        return _check(pf, strict)
    lower, upper = _sin_walls(w, amp, lam, 2000)
    lc = _polyline_length(np.column_stack([np.linspace(0, lam, 4001), amp * np.cos(k * np.linspace(0, lam, 4001))]))

    def fin_thickness(p):  # min distance between this channel's upper wall and the next channel's lower wall
        a = np.vstack([upper[:-1] - [lam, 0], upper[:-1], upper[:-1] + [lam, 0]])
        b = lower[::4] + [0, p]
        from scipy.spatial import cKDTree
        return float(cKDTree(a).query(b)[0].min())

    phi = math.atan(amp * k)
    lo, hi = w + 0.5 * t, (w + t) / math.cos(phi) + t
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if fin_thickness(mid) < t else (lo, mid)
    p = 0.5 * (lo + hi)
    v += _limits(h_fin, t, w, t)
    pf = Planform("sinusoidal", {"w": w, "t": t, "amp": amp, "lam": lam, "h_fin": h_fin}, "channel", period_mm=lam,
                  lateral_unit_mm=p, wetted_perimeter_mm=_polyline_length(lower) + _polyline_length(upper),
                  fluid_area_mm2=w * lc, fin_section_perimeter_mm=2 * lc, fin_section_area_mm2=p * lam - w * lc,
                  min_solid_mm=t, min_gap_mm=w, min_fin_thickness_mm=t, ref_gap_mm=w, violations=tuple(v),
                  notes=_channel_aspect_note(w, h_fin) + (f"lateral pitch {p:.4f} mm; max wall angle {math.degrees(phi):.1f} deg",),
                  geo={"w": w, "amp": amp, "lam": lam})
    return _check(pf, strict)


def _zigzag_parts(w, amp, lam, r_in):
    """Centreline geometry of one period starting at a (rounded) crest: slope angle, centreline radius, arc centres."""
    phi = math.atan(4 * amp / lam)
    r = r_in + w / 2
    top_c = (0.0, amp - r / math.cos(phi))
    bot_c = (lam / 2, -amp + r / math.cos(phi))
    return phi, r, top_c, bot_c


def zigzag(w: float, t: float, amp: float, lam: float, h_fin: float, strict: bool = True,
           r_in: float = ZIGZAG_INNER_RADIUS_MM) -> Planform:
    """Triangular-wave channel (amplitude amp of the sharp wave, wavelength lam), corners rounded to inner radius r_in."""
    phi, r, top_c, bot_c = _zigzag_parts(w, amp, lam, r_in)
    v = []
    # the two corner arcs must not overlap: tangent points must lie on the segment between corners
    seg = math.hypot(lam / 2, 2 * amp)
    tangent_len = r * math.tan(math.pi / 2 - (math.pi - 2 * phi) / 2)  # r * tan(half the turning angle) = r tan(phi)
    if 2 * tangent_len > seg - 1e-9:
        v.append("corner fillets overlap (amplitude too small for the rounding radius)")
    p = (w + t) / math.cos(phi)
    lc = 2 * (seg - 2 * tangent_len) + 2 * r * (2 * phi)  # per period: two straight runs, two arcs of angle 2 phi
    v += _limits(h_fin, t, w, t)
    pf = Planform("zigzag", {"w": w, "t": t, "amp": amp, "lam": lam, "h_fin": h_fin, "r_in": r_in}, "channel", period_mm=lam,
                  lateral_unit_mm=p, wetted_perimeter_mm=2 * lc, fluid_area_mm2=w * lc, fin_section_perimeter_mm=2 * lc,
                  fin_section_area_mm2=p * lam - w * lc, min_solid_mm=t, min_gap_mm=w, min_fin_thickness_mm=t, ref_gap_mm=w,
                  violations=tuple(v), notes=_channel_aspect_note(w, h_fin) + (
                      f"lateral pitch {p:.4f} mm; wall angle {math.degrees(phi):.1f} deg; corners rounded, inner radius {r_in} mm",),
                  geo={"w": w, "amp": amp, "lam": lam, "r_in": r_in})
    return _check(pf, strict)


def offset_strip(l: float, t: float, p_t: float, g: float, h_fin: float, strict: bool = True) -> Planform:
    """Offset strip fins: strips l x t, lateral pitch p_t, alternate rows shifted p_t/2 laterally and l + g streamwise."""
    half = p_t / 2
    lat_gap = half - t
    v = []
    if lat_gap <= 0:
        v.append("strips of adjacent rows overlap laterally (p_t/2 <= t); not supported")
        min_gap = g
    else:
        min_gap = min(p_t - t, math.hypot(g, lat_gap), l + 2 * g)
    v += _limits(h_fin, min(t, l), min_gap, t)
    period = 2 * (l + g)
    pf = Planform("offset_strip", {"l": l, "t": t, "p_t": p_t, "g": g, "h_fin": h_fin}, "strip", period_mm=period,
                  lateral_unit_mm=half, wetted_perimeter_mm=2 * (l + t), fluid_area_mm2=period * half - l * t,
                  fin_section_perimeter_mm=2 * (l + t), fin_section_area_mm2=l * t, min_solid_mm=min(t, l), min_gap_mm=min_gap,
                  min_fin_thickness_mm=t, ref_gap_mm=max(lat_gap, MIN_GAP_MM / 2), violations=tuple(v),
                  notes=(f"in-row channel width {p_t - t:.3f} mm, lateral gap between rows {lat_gap:.3f} mm",),
                  geo={"l": l, "t": t, "half": half, "g": g})
    return _check(pf, strict)


def _pins(shape: str, size: float, s_t: float, s_l: float, h_fin: float, strict: bool) -> Planform:
    half = s_t / 2
    if shape == "circle":
        a = size / 2  # half extent along x and y
        mk = lambda c: _circle(c, a)
        perim, area, thick = math.pi * size, math.pi * size**2 / 4, size
    else:
        a = size / math.sqrt(2)  # half diagonal
        mk = lambda c: _diamond(c, a)
        perim, area, thick = 4 * size, size**2, size
    ref = mk((0.0, 0.0))
    neighbours = [(0.0, s_t), (s_l, half), (2 * s_l, 0.0)]
    gaps = [_poly_dist(ref, mk(c)) for c in neighbours] if shape == "diamond" else \
        [s_t - size, math.hypot(s_l, half) - size, 2 * s_l - size]
    v = []
    if s_l <= a + 1e-9:
        v.append("row-B pins reach the periodic cell boundary (S_L <= half pin extent); not supported")
    if a >= half - 1e-9:
        v.append("pins of the same row touch (S_T <= pin extent)")
    min_gap = float(min(gaps))
    v += _limits(h_fin, size, min_gap, thick)
    topo = "pin_circle" if shape == "circle" else "pin_diamond"
    key = "d" if shape == "circle" else "s"
    pf = Planform(topo, {key: size, "s_t": s_t, "s_l": s_l, "h_fin": h_fin}, "strip", period_mm=2 * s_l, lateral_unit_mm=half,
                  wetted_perimeter_mm=perim, fluid_area_mm2=2 * s_l * half - area, fin_section_perimeter_mm=perim,
                  fin_section_area_mm2=area, min_solid_mm=size, min_gap_mm=min_gap, min_fin_thickness_mm=thick,
                  ref_gap_mm=max(min_gap, MIN_GAP_MM / 2), violations=tuple(v),
                  notes=(f"gaps: transverse {gaps[0]:.3f}, diagonal {gaps[1]:.3f}, streamwise {gaps[2]:.3f} mm",),
                  geo={"shape": shape, "a": a, "half": half, "s_l": s_l})
    return _check(pf, strict)


def pin_circle(d: float, s_t: float, s_l: float, h_fin: float, strict: bool = True) -> Planform:
    return _pins("circle", d, s_t, s_l, h_fin, strict)


def pin_diamond(s: float, s_t: float, s_l: float, h_fin: float, strict: bool = True) -> Planform:
    return _pins("diamond", s, s_t, s_l, h_fin, strict)


BUILDERS = {"straight": straight, "sinusoidal": sinusoidal, "zigzag": zigzag, "offset_strip": offset_strip,
            "pin_circle": pin_circle, "pin_diamond": pin_diamond}


def build(topology: str, params: dict, strict: bool = True) -> Planform:
    return BUILDERS[topology](**params, strict=strict)


# --------------------------------------------------------------------------------------------------------------------
# Meshing


@dataclass
class PlanMesh:
    mesh: MeshTri
    inlet: np.ndarray  # periodic: master boundary
    outlet: np.ndarray  # periodic: slave boundary
    walls: np.ndarray  # all no-slip fin walls
    heated: np.ndarray  # heated fin walls (excludes the lead-in of a developing channel)
    symmetry: np.ndarray
    cuts: np.ndarray  # x of period boundaries (mesh conforms to them)
    periodic: bool
    flow_rate: float  # nondimensional, per unit depth, with velocity scale = superficial velocity, length scale = 1 mm
    inlet_width: float
    x_in: float
    x_out: float


def _add_poly(occ, pts):
    tags = [occ.addPoint(x, y, 0) for x, y in pts]
    lines = [occ.addLine(tags[i], tags[(i + 1) % len(tags)]) for i in range(len(tags))]
    return occ.addPlaneSurface([occ.addCurveLoop(lines)])


def _channel_surfaces(pf: Planform, occ, x_start: float, n_periods: int):
    """Fluid surfaces of a wavy channel over n periods starting at x_start (a crest)."""
    g = pf.geo
    lam, w = g["lam"], g["w"]
    surfaces = []
    for k in range(n_periods):
        x0 = x_start + k * lam
        if pf.topology == "sinusoidal":
            lower, upper = _sin_walls(w, g["amp"], lam, 200)
            lo = [occ.addPoint(x + x0, y, 0) for x, y in lower]
            up = [occ.addPoint(x + x0, y, 0) for x, y in upper]
            c_lo, c_up = occ.addSpline(lo), occ.addSpline(up)
            loop = [c_lo, occ.addLine(lo[-1], up[-1]), -c_up, occ.addLine(up[0], lo[0])]
            surfaces.append(occ.addPlaneSurface([occ.addCurveLoop(loop)]))
        else:
            surfaces.append(_zigzag_surface(occ, pf, x0))
    return surfaces


def _zigzag_surface(occ, pf, x0):
    g = pf.geo
    w, amp, lam, r_in = g["w"], g["amp"], g["lam"], g["r_in"]
    phi, r, top_c, bot_c = _zigzag_parts(w, amp, lam, r_in)
    sp, cp = math.sin(phi), math.cos(phi)

    def wall(side):  # side = +1 upper wall, -1 lower wall; returns curve tags from x0 to x0 + lam and end point tags
        rt = r + side * w / 2  # top corner: upper wall is outside (centre below)
        rb = r - side * w / 2  # bottom corner: upper wall is inside (centre above)
        pts = []
        ct = occ.addPoint(x0 + top_c[0], top_c[1], 0)
        cb = occ.addPoint(x0 + bot_c[0], bot_c[1], 0)
        ct2 = occ.addPoint(x0 + lam + top_c[0], top_c[1], 0)
        p0 = occ.addPoint(x0, top_c[1] + rt, 0)
        p1 = occ.addPoint(x0 + top_c[0] + rt * sp, top_c[1] + rt * cp, 0)
        p2 = occ.addPoint(x0 + bot_c[0] - rb * sp, bot_c[1] - rb * cp, 0)
        p3 = occ.addPoint(x0 + bot_c[0] + rb * sp, bot_c[1] - rb * cp, 0)
        p4 = occ.addPoint(x0 + lam - rt * sp, top_c[1] + rt * cp, 0)
        p5 = occ.addPoint(x0 + lam, top_c[1] + rt, 0)
        curves = [occ.addCircleArc(p0, ct, p1), occ.addLine(p1, p2), occ.addCircleArc(p2, cb, p3), occ.addLine(p3, p4),
                  occ.addCircleArc(p4, ct2, p5)]
        pts = (p0, p5)
        return curves, pts

    lo, (lo0, lo1) = wall(-1)
    up, (up0, up1) = wall(+1)
    loop = lo + [occ.addLine(lo1, up1)] + [-c for c in reversed(up)] + [occ.addLine(up0, lo0)]
    return occ.addPlaneSurface([occ.addCurveLoop(loop)])


def _strip_solids(pf: Planform, x_from: float, x_to: float) -> list[np.ndarray]:
    """Solid polygons (or circle specs) of a strip topology with centres in [x_from, x_to]."""
    g = pf.geo
    out = []
    if pf.topology == "offset_strip":
        l, t, half, gap = g["l"], g["t"], g["half"], g["g"]
        period = 2 * (l + gap)
        k0, k1 = int(math.floor(x_from / period)) - 1, int(math.ceil(x_to / period)) + 1
        for k in range(k0, k1 + 1):
            for xc, yc in ((k * period, 0.0), (k * period + l + gap, half)):
                if x_from - 1e-9 <= xc <= x_to + 1e-9:
                    out.append(("poly", _rect(xc - l / 2, xc + l / 2, yc - t / 2, yc + t / 2)))
    else:
        a, half, s_l = g["a"], g["half"], g["s_l"]
        k0, k1 = int(math.floor(x_from / (2 * s_l))) - 1, int(math.ceil(x_to / (2 * s_l))) + 1
        for k in range(k0, k1 + 1):
            for xc, yc in ((2 * k * s_l, 0.0), (2 * k * s_l + s_l, half)):
                if x_from - 1e-9 <= xc <= x_to + 1e-9:
                    out.append(("circle", (xc, yc, a)) if g["shape"] == "circle" else ("poly", _diamond((xc, yc), a)))
    return out


def _to_meshtri() -> MeshTri:
    tags, coords, _ = gmsh.model.mesh.getNodes()
    _, _, conn = gmsh.model.mesh.getElements(2)
    conn = np.concatenate([c.astype(np.int64) for c in conn])
    used = np.unique(conn)
    remap = np.full(int(tags.max()) + 1, -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    idx = np.zeros(int(tags.max()) + 1, dtype=np.int64)
    idx[tags.astype(np.int64)] = np.arange(len(tags))
    p = coords.reshape(-1, 3)[idx[used], :2].T
    t = remap[conn].reshape(-1, 3).T
    a, b, c = p[:, t[0]], p[:, t[1]], p[:, t[2]]
    det = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    t[1, det < 0], t[2, det < 0] = t[2, det < 0].copy(), t[1, det < 0].copy()  # counter-clockwise
    return MeshTri(np.ascontiguousarray(p), np.ascontiguousarray(t))


def mesh_planform(pf: Planform, n_periods: int = 1, periodic: bool = True, cells: float = 8.0,
                  lead_in_mm: float = 1.0, outlet_mm: float | None = None) -> PlanMesh:
    """Mesh the fluid region. Wall size = ref_gap / (2 cells), core size = ref_gap / cells, graded over ref_gap / 2.

    Periodic: x in [0, n_periods * period], master at x = 0, slave at the end, internal period lines embedded.
    Developing: lead-in (adiabatic, no fins for strips; an unheated straight channel for channels), n_periods heated
    periods with embedded period lines at x = 0, period, ..., and an outlet buffer (strips: fins end at the last row,
    then empty; channels: an unheated straight channel).
    """
    period = pf.period_mm
    x_end = n_periods * period
    outlet_mm = period / 2 if outlet_mm is None else outlet_mm
    h_wall, h_core = pf.ref_gap_mm / (2 * cells), pf.ref_gap_mm / cells
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("planform")
        occ = gmsh.model.occ
        cuts = [k * period for k in range(n_periods + 1)]
        if pf.kind == "half_channel":
            y0, y1 = 0.0, pf.geo["half_width"]
            x_in, x_out = (0.0, x_end) if periodic else (-lead_in_mm, x_end + outlet_mm)
            surfaces = [occ.addRectangle(x_in, y0, 0, x_out - x_in, y1 - y0)]
        elif pf.kind == "channel":
            g = pf.geo
            surfaces = _channel_surfaces(pf, occ, 0.0, n_periods)
            if pf.topology == "sinusoidal":
                y_lo, y_hi = g["amp"] - g["w"] / 2, g["amp"] + g["w"] / 2
            else:
                _, r, top_c, _ = _zigzag_parts(g["w"], g["amp"], g["lam"], g["r_in"])
                y_lo, y_hi = top_c[1] + r - g["w"] / 2, top_c[1] + r + g["w"] / 2
            y0, y1 = y_lo, y_hi
            x_in, x_out = 0.0, x_end
            if not periodic:
                x_in, x_out = -lead_in_mm, x_end + outlet_mm
                surfaces += [occ.addRectangle(x_in, y_lo, 0, lead_in_mm, y_hi - y_lo),
                             occ.addRectangle(x_end, y_lo, 0, outlet_mm, y_hi - y_lo)]
        else:
            y0, y1 = 0.0, pf.lateral_unit_mm
            if periodic:
                x_in, x_out = 0.0, x_end
                solids = _strip_solids(pf, -period, x_end + period)
            else:
                ext = pf.geo["l"] / 2 if pf.topology == "offset_strip" else pf.geo["a"]
                x_in, x_out = -ext - lead_in_mm, x_end + ext + outlet_mm
                solids = _strip_solids(pf, 0.0, x_end)
            rect = occ.addRectangle(x_in, y0, 0, x_out - x_in, y1 - y0)
            tools = [(2, occ.addDisk(s[0], s[1], 0, s[2], s[2])) if kind == "circle" else (2, _add_poly(occ, s))
                     for kind, s in solids]
            out, _ = occ.cut([(2, rect)], tools)
            surfaces = [tag for dim, tag in out if dim == 2]
        occ.synchronize()
        internal = [c for c in cuts if (c > x_in + 1e-9 and c < x_out - 1e-9)]
        lines = []
        for c in internal:
            a, b = occ.addPoint(c, y0 - 1.0, 0), occ.addPoint(c, y1 + 1.0, 0)
            lines.append((1, occ.addLine(a, b)))
        occ.fragment([(2, s) for s in surfaces], lines)
        occ.synchronize()
        # drop curve pieces outside the fluid
        on_surface = set()
        for _, s in gmsh.model.getEntities(2):
            on_surface.update(abs(c) for c in gmsh.model.getAdjacencies(2, s)[1])
        stray = [(1, c) for _, c in gmsh.model.getEntities(1) if c not in on_surface]
        if stray:
            gmsh.model.removeEntities(stray, recursive=True)
        eps = 1e-6
        if periodic:
            def vertical(x):
                ents = gmsh.model.getEntitiesInBoundingBox(x - eps, y0 - 10, -eps, x + eps, y1 + 10, eps, dim=1)
                return sorted((t for _, t in ents), key=lambda t: gmsh.model.getBoundingBox(1, t)[1])
            m_tags, s_tags = vertical(x_in), vertical(x_out)
            if len(m_tags) != len(s_tags):
                raise GeometryRejected(["periodic boundaries do not match"])
            affine = [1, 0, 0, x_out - x_in, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
            gmsh.model.mesh.setPeriodic(1, s_tags, m_tags, affine)
        wall_curves = []
        for _, c in gmsh.model.getEntities(1):
            bb = gmsh.model.getBoundingBox(1, c)
            is_vertical_special = abs(bb[3] - bb[0]) < eps and any(abs(bb[0] - x) < eps for x in [x_in, x_out] + internal)
            is_sym = abs(bb[4] - bb[1]) < eps and ((pf.kind == "strip" and (abs(bb[1] - y0) < eps or abs(bb[1] - y1) < eps))
                                                    or (pf.kind == "half_channel" and abs(bb[1] - y1) < eps))
            if not is_vertical_special and not is_sym:
                wall_curves.append(c)
        fd = gmsh.model.mesh.field
        dist = fd.add("Distance")
        fd.setNumbers(dist, "CurvesList", wall_curves)
        fd.setNumber(dist, "Sampling", 400)
        thr = fd.add("Threshold")
        fd.setNumber(thr, "InField", dist)
        fd.setNumber(thr, "SizeMin", h_wall)
        fd.setNumber(thr, "SizeMax", h_core)
        fd.setNumber(thr, "DistMin", 0.0)
        fd.setNumber(thr, "DistMax", pf.ref_gap_mm / 2)
        fd.setAsBackgroundMesh(thr)
        for opt in ("Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature", "Mesh.MeshSizeExtendFromBoundary"):
            gmsh.option.setNumber(opt, 0)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.model.mesh.generate(2)
        mesh = _to_meshtri()
    finally:
        gmsh.finalize()

    tol = 1e-6
    bnd = mesh.boundary_facets()
    mid = mesh.p[:, mesh.facets[:, bnd]].mean(axis=1)
    inlet = bnd[np.abs(mid[0] - x_in) < tol]
    outlet = bnd[np.abs(mid[0] - x_out) < tol]
    if pf.kind == "strip":
        sym_mask = (np.abs(mid[1] - y0) < tol) | (np.abs(mid[1] - y1) < tol)
    elif pf.kind == "half_channel":
        sym_mask = np.abs(mid[1] - y1) < tol
    else:
        sym_mask = np.zeros(len(bnd), bool)
    rest = ~sym_mask & (np.abs(mid[0] - x_in) >= tol) & (np.abs(mid[0] - x_out) >= tol)
    walls = bnd[rest]
    if periodic or pf.kind == "strip":
        heated = walls
    else:
        wm = mid[0][rest]
        heated = walls[(wm > 0) & (wm < x_end)]
    inlet_width = float(np.sum(np.abs(np.diff(mesh.p[1, mesh.facets[:, inlet]], axis=0))))
    return PlanMesh(mesh, inlet, outlet, walls, heated, bnd[sym_mask], np.array(cuts), periodic, pf.lateral_unit_mm,
                    inlet_width, x_in, x_out)
