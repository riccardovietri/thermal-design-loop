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

TEXT = {"README.md": README, "docs/index.html": INDEX, "ASSUMPTIONS.md": ROOT / "ASSUMPTIONS.md"}


def load_texts() -> dict:
    return {name: path.read_text() for name, path in TEXT.items()}


def failures(texts: dict) -> list[str]:
    """Every disagreement between the published text (name -> raw text) and the committed results."""
    c, out = ctx(), []
    flat = {name: " ".join(t.split()) for name, t in texts.items()}
    for label, files, build in CHECKS:
        expected = build(c)
        for name in files:
            if expected not in flat[name]:
                out.append(f"{name}: {label} should read '{expected}' (from the committed results)")
    for label, name, row_start, cells in PER_TOPOLOGY:
        html = flat[name]
        start = html.find(row_start(c))
        if start < 0:
            out.append(f"{name}: no table row for {label}")
            continue
        row = html[start:html.find("</tr>", start)]
        if cells(c) not in row:
            out.append(f"{name}: {label} should read '{cells(c)}' in its table row; found '{row}'")
    for name in ("README.md", "docs/index.html"):
        for value in SUPERSEDED:
            if re.search(rf"(?<![\d.]){re.escape(value)}(?![\d])", flat[name]):
                out.append(f"{name}: superseded value '{value}' from an earlier study still appears")
    bad = [i + 1 for i, line in enumerate(texts["docs/index.html"].splitlines()) if "—" in line or "–" in line]
    if bad:
        out.append(f"docs/index.html: em or en dash on lines {bad}; expected none (write ranges as 'a to b')")
    return out


def test_published_text_matches_results():
    found = failures(load_texts())
    assert not found, "\n".join(found)


# Corruptions of the published text that the checks must catch: (id, file, exact old text, new text).
MUTATIONS = [
    ("README hero R_th gains a leading digit", "README.md", "3.75 mm tall. 14.0 mK/W at 1.60 kPa", "3.75 mm tall. 114.0 mK/W at 1.60 kPa"),
    ("README entrance margin 88.9%", "README.md", "objective margin 8.9%", "objective margin 88.9%"),
    ("README second parity claim disagrees", "README.md", "to 1.0% on average and 2.9% at worst", "to 1.0% on average and 3.9% at worst"),
    ("README entrance dP 1.627 kPa", "README.md", "from 1.600 to 1.617 kPa", "from 1.600 to 1.627 kPa"),
    ("README entrance R_th 13.56 at two decimals", "README.md", "from 14.05 to 13.66 mK/W", "from 14.05 to 13.56 mK/W"),
    ("README pin pitch 0.8 x 1.0 mm", "README.md", "at 0.8 × 0.9 mm pitch", "at 0.8 × 1.0 mm pitch"),
    ("README screened count", "README.md", "screens 4,485 fin geometries", "screens 4,486 fin geometries"),
    ("README CFD count in hero table", "README.md", "2.9% at worst (18 designs)", "2.9% at worst (19 designs)"),
    ("README duct bound below the saved error", "README.md", "Nu to within 0.008%", "Nu to within 0.007%"),
    ("README reference row R_th", "README.md", "| 0.42 / 0.39 / 3.75 mm | 14.9 mK/W | 1.29 kPa |",
     "| 0.42 / 0.39 / 3.75 mm | 13.9 mK/W | 1.29 kPa |"),
    ("README 2.0 LPM flow sensitivity", "README.md", "falls to 12.8 mK/W (52.8 °C)", "falls to 12.3 mK/W (52.8 °C)"),
    ("README plan-view straight pick", "README.md", "(13.4 mK/W at 1.31 kPa)", "(13.4 mK/W at 1.21 kPa)"),
    ("index entrance dP 1.72 kPa", "docs/index.html", "from 1.60 to 1.62 kPa", "from 1.60 to 1.72 kPa"),
    ("index duct bound below the saved error", "docs/index.html", "within 0.008%", "within 0.007%"),
    ("index energy-balance bound below the saved error", "docs/index.html", "within 0.76%", "within 0.75%"),
    ("index mesh change", "docs/index.html", "0.002% change", "0.02% change"),
    ("index wetted-area range", "docs/index.html", "6.2 to 10.3 times", "6.2 to 8.3 times"),
    ("index widest-channel miss", "docs/index.html", "misses by 10.8%", "misses by 1.8%"),
    ("index walk table row", "docs/index.html", '0.46 / 0.44 / 4.0 mm</td><td class="n">15.9 mK/W',
     '0.46 / 0.44 / 4.0 mm</td><td class="n">16.9 mK/W'),
    ("index flow table 1.0 LPM resistance", "docs/index.html", '1.07 kPa</td><td class="n">16.5 mK/W',
     '1.07 kPa</td><td class="n">15.5 mK/W'),
    ("index pin pitch", "docs/index.html", "0.4 mm pins, 0.8 × 0.9 mm pitch", "0.4 mm pins, 0.9 × 0.9 mm pitch"),
    ("index plan-view bias in prose", "docs/index.html", "resistance 5.6% lower than", "resistance 4.6% lower than"),
    ("index 1D split part", "docs/index.html", "about 7.5 for fins", "about 7.8 for fins"),
    ("index pin advantage range", "docs/index.html", "2.6 to 5.8% lower resistance", "2.6 to 6.8% lower resistance"),
    ("index margin in the introduction", "docs/index.html", "with a 6.3% margin, but", "with a 6.8% margin, but"),
    ("index allowable rise for capacity", "docs/index.html", "at the 20 K allowable", "at the 25 K allowable"),
    ("index entrance gate agreement", "docs/index.html", "to within 3.1% and 2.9%", "to within 3.1% and 1.9%"),
    ("ASSUMPTIONS heated width", "ASSUMPTIONS.md", "59.64 mm (77 × 0.775 mm)", "59.84 mm (77 × 0.775 mm)"),
]


@pytest.mark.parametrize("name, old, new", [m[1:] for m in MUTATIONS], ids=[m[0] for m in MUTATIONS])
def test_checks_catch_mutation(name, old, new):
    texts = load_texts()
    assert texts[name].count(old) == 1, f"mutation target not found exactly once in {name}: {old!r}"
    texts[name] = texts[name].replace(old, new)
    assert failures(texts), f"{name}: corrupting {old!r} to {new!r} was not caught"


def test_index_images_exist():
    refs = re.findall(r'src="(img/[^"]+)"', INDEX.read_text())
    assert refs, "docs/index.html: no img/ references found"
    missing = [r for r in refs if not (INDEX.parent / r).exists()]
    assert not missing, f"docs/index.html: referenced images missing from docs/: {missing}"


def test_readme_images_exist():
    refs = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", README.read_text())
    missing = [r for r in refs if not (ROOT / r).exists()]
    assert not missing, f"README.md: referenced images missing: {missing}"
