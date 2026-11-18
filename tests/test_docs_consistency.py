"""Published text (README.md, docs/index.html) against the committed results.

Every expected string is derived from the committed results files (or, for the screen counts, from re-running the analytic
1D screen on the committed requirements file, as `run_study.py --plot-only` does). Nothing the study produces is hard-coded
here except the short list of superseded values that must not reappear.
"""
import json
import re
from functools import lru_cache
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
INDEX = ROOT / "docs" / "index.html"
LEVEL1 = ROOT / "examples" / "output" / "gpu_coldplate" / "results.json"
PLANVIEW = ROOT / "examples" / "output" / "level2_planview" / "results.json"
ENTRANCE = ROOT / "examples" / "output" / "entrance_effects" / "results.json"
G3 = ROOT / "validation" / "level2_planview_straight.json"
CONFIG = ROOT / "examples" / "gpu_coldplate.yaml"

SUPERSEDED = ("14.93", "672", "258")  # values of earlier studies that must not reappear in the published text
INDEX_LABELS = {"pin_circle": "Circular pins", "straight": "Straight", "zigzag": "Zigzag", "sinusoidal": "Sinusoidal",
                "pin_diamond": "Diamond pins", "offset_strip": "Offset strip"}


def _flat(path: Path) -> str:
    return " ".join(path.read_text().split())


@lru_cache(maxsize=None)
def ctx() -> dict:
    from run_study import load
    from tools.fluids import fluid_props
    from tools.screen_1d import screen_1d

    l1 = json.loads(LEVEL1.read_text())
    pv = json.loads(PLANVIEW.read_text())
    en = json.loads(ENTRANCE.read_text())
    g3 = json.loads(G3.read_text())
    _, req, space, _ = load(CONFIG)
    sweep = screen_1d(req, fluid_props(req.coolant, req.inlet_temp_c), space)
    rec = next(d for d in l1["designs"] if d["cfd"]["design_id"] == l1["recommended"])
    return {"l1": l1, "pv": pv, "en": en, "g3": g3, "rec": rec, "n_sweep": len(sweep),
            "n_feasible": sum(s.feasible for s in sweep), "objective": l1["r_th_objective_k_per_w"]}


def _rec_geometry(c):
    g = c["rec"]["one_d"]["design"]
    return f"{g['w_ch_m']*1e3:.2f} mm channels, {g['t_fin_m']*1e3:.2f} mm fins, {g['h_fin_m']*1e3:.2f} mm tall"


def _flow_row_html(c, q):
    r = next(x for x in c["l1"]["flow_sensitivity"] if abs(x["flow_lpm"] - q) < 1e-9)
    return (f'<td class="n">{q:.1f} LPM</td><td class="n">{r["dp_pa"]/1e3:.2f} kPa</td>'
            f'<td class="n">{r["r_th_k_per_w"]*1e3:.1f} mK/W</td><td class="n">{r["base_avg_temp_c"]:.1f} °C</td>'
            f'<td class="n">{r["base_max_temp_c"]:.1f} °C</td>')


def _flow_readme(c, q):
    r = next(x for x in c["l1"]["flow_sensitivity"] if abs(x["flow_lpm"] - q) < 1e-9)
    return f"{r['r_th_k_per_w']*1e3:.1f} mK/W ({r['base_avg_temp_c']:.1f} °C)"


def _pv_pick(c, topology=None):
    s = c["pv"]["summary"]
    p = s["recommended"] if topology is None else s["per_topology"][topology]["recommended"]
    return p["r_th_k_per_w"] * 1e3, p["dp_pa"] / 1e3


def _wetted_ratios(c):
    """Heated wetted area (channel floor plus both fin sides) over the modelled heated footprint (n x pitch x length)."""
    out = []
    for d in c["l1"]["designs"]:
        g = d["one_d"]["design"]
        w, t, h = g["w_ch_m"], g["t_fin_m"], g["h_fin_m"]
        out.append((d["cfd"]["design_id"], (w + 2 * h) / (w + t)))
    return out


def _heated_width_mm(c):
    g, n = c["rec"]["one_d"]["design"], c["rec"]["one_d"]["n_channels"]
    return n * (g["w_ch_m"] + g["t_fin_m"]) * 1e3


def _entrance_rec(c):
    en = c["en"]
    row = next(r for r in en["one_d"] if r["design_id"] == en["level1_recommended"])
    return row["cfd_anchored"]


def _entrance_widest_r_th(c):
    res = c["en"]["resolved"]
    return 100 * res[0]["comparison"]["r_th_contribution_disagreement"]


BOTH, R, I = ("README.md", "docs/index.html"), ("README.md",), ("docs/index.html",)

# (label, files, expected-string builder)
CHECKS = [
    ("recommended geometry", BOTH, _rec_geometry),
    ("recommended R_th", BOTH, lambda c: f"{c['rec']['cfd']['r_th_k_per_w']*1e3:.1f} mK/W at {c['rec']['cfd']['dp_pa']/1e3:.2f} kPa"),
    ("recommended ΔP", BOTH, lambda c: f"{c['rec']['cfd']['dp_pa']/1e3:.2f} kPa"),
    ("objective margin", I, lambda c: f"{100*(1 - c['rec']['cfd']['r_th_k_per_w']/c['objective']):.1f}% margin"),
    ("designs screened", BOTH, lambda c: f"{c['n_sweep']:,}"),
    ("designs feasible", I, lambda c: f"{c['n_feasible']:,} of them"),
    ("designs solved with CFD", BOTH, lambda c: f"{len(c['l1']['designs'])} designs"),
    ("1D vs CFD mean difference", BOTH, lambda c: f"{c['l1']['parity']['mean_abs_pct']:.1f}% on average"),
    ("1D vs CFD worst difference", BOTH, lambda c: f"{c['l1']['parity']['max_abs_pct']:.1f}% at worst"),
    ("1D resistance split total", I, lambda c: f"{c['rec']['one_d']['r_th_k_per_w']*1e3:.1f} mK/W"),
    ("1D split: fins and convection", I, lambda c: f"{c['rec']['one_d']['r_conv']*1e3:.1f} for fins and convection"),
    ("1D split: coolant", I, lambda c: f"{c['rec']['one_d']['r_caloric']*1e3:.1f} for the coolant"),
    ("1D split: base", I, lambda c: f"{c['rec']['one_d']['r_base']*1e3:.1f} for conduction through the base"),
    ("flow sensitivity 1.0 LPM", R, lambda c: _flow_readme(c, 1.0)),
    ("flow sensitivity 2.0 LPM", R, lambda c: _flow_readme(c, 2.0)),
    ("flow sensitivity row 1.0 LPM", I, lambda c: _flow_row_html(c, 1.0)),
    ("flow sensitivity row 1.5 LPM", I, lambda c: _flow_row_html(c, 1.5)),
    ("flow sensitivity row 2.0 LPM", I, lambda c: _flow_row_html(c, 2.0)),
    ("cross-topology pick", R, lambda c: "{:.1f} mK/W at {:.2f} kPa".format(*_pv_pick(c))),
    ("cross-topology pick", I, lambda c: '{:.1f} mK/W</td><td class="n">{:.2f} kPa'.format(*_pv_pick(c))),
    ("plan-view tall-fin bias on R_th", I, lambda c: f"{-100*c['g3']['g3_r_th_k_per_w']['difference']:.1f}%"),
    ("plan-view tall-fin bias on ΔP", I, lambda c: f"{-100*c['g3']['g3_dp_pa']['difference']:.1f}% low"),
    ("wetted area ratio range", I, lambda c: "{:.1f} to {:.1f} times".format(min(r for _, r in _wetted_ratios(c)),
                                                                             max(r for _, r in _wetted_ratios(c)))),
    ("wetted area ratio, recommended", I, lambda c: "{:.1f} for the recommended design".format(
        dict(_wetted_ratios(c))[c["l1"]["recommended"]])),
    ("entrance-corrected R_th", BOTH, lambda c: f"{_entrance_rec(c)['r_th_combined']*1e3:.1f} mK/W"),
    ("entrance-corrected objective margin", BOTH,
     lambda c: f"{100*(1 - _entrance_rec(c)['r_th_combined']/c['objective']):.1f}%"),
    ("entrance check, widest channel miss on R_th", I, lambda c: f"{abs(_entrance_widest_r_th(c)):.1f}%"),
]

PER_TOPOLOGY = [(f"plan-view pick, {t}", I[0], (lambda t: lambda c: '<td>{}</td><td>'.format(INDEX_LABELS[t]))(t),
                 (lambda t: lambda c: '{:.1f} mK/W</td><td class="n">{:.2f} kPa</td>'.format(*_pv_pick(c, t)))(t))
                for t in INDEX_LABELS]

TEXT = {"README.md": README, "docs/index.html": INDEX}


@pytest.mark.parametrize("label, name, build", [(lab, f, b) for lab, files, b in CHECKS for f in files],
                         ids=[f"{lab} | {f}" for lab, files, _ in CHECKS for f in files])
def test_text_quotes_results(label, name, build):
    expected = build(ctx())
    assert expected in _flat(TEXT[name]), f"{name}: {label} should read '{expected}' (from the committed results)"


@pytest.mark.parametrize("label, name, row_start, cells", PER_TOPOLOGY, ids=[p[0] for p in PER_TOPOLOGY])
def test_index_plan_view_table(label, name, row_start, cells):
    c = ctx()
    html = _flat(TEXT[name])
    start = html.find(row_start(c))
    assert start >= 0, f"{name}: no table row for {label}"
    row = html[start:html.find("</tr>", start)]
    assert cells(c) in row, f"{name}: {label} should read '{cells(c)}' in its table row; found '{row}'"


@pytest.mark.parametrize("name", TEXT)
@pytest.mark.parametrize("value", SUPERSEDED)
def test_no_superseded_values(name, value):
    hits = re.findall(rf"(?<![\d.]){re.escape(value)}(?![\d])", _flat(TEXT[name]))
    assert not hits, f"{name}: superseded value '{value}' from an earlier study still appears"


def test_index_has_no_em_or_en_dash():
    text = INDEX.read_text()
    bad = [(i + 1, line.strip()[:80]) for i, line in enumerate(text.splitlines()) if "—" in line or "–" in line]
    assert not bad, f"docs/index.html: em or en dash on lines {bad}; expected none (write ranges as 'a to b')"


def test_index_images_exist():
    refs = re.findall(r'src="(img/[^"]+)"', INDEX.read_text())
    assert refs, "docs/index.html: no img/ references found"
    missing = [r for r in refs if not (INDEX.parent / r).exists()]
    assert not missing, f"docs/index.html: referenced images missing from docs/: {missing}"


def test_readme_images_exist():
    refs = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", README.read_text())
    missing = [r for r in refs if not (ROOT / r).exists()]
    assert not missing, f"README.md: referenced images missing: {missing}"
