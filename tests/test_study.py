import json

import pytest
import yaml

from run_study import plot_only, run_study
from tools.plotstyle import round_ticks
from tools.report import parity_stats


def test_parity_stats():
    s = parity_stats([1.1, 0.9, 1.0], [1.0, 1.0, 1.0])
    assert s["mean_pct"] == pytest.approx(0.0, abs=1e-9)
    assert s["mean_abs_pct"] == pytest.approx(20 / 3)
    assert s["max_abs_pct"] == pytest.approx(10.0)


def test_run_study_end_to_end(tmp_path):
    cfg = yaml.safe_load(open("examples/gpu_coldplate.yaml"))
    cfg["name"] = "smoke"
    cfg["design_space"].update({"w_ch_mm": [0.3, 0.6], "t_fin_mm": [0.3], "h_fin_mm": [2.0, 3.0], "flow_lpm": [1.5]})
    cfg["study"] = {"n_cfd": 3, "mesh_cells_across_half_channel": 4, "flow_sensitivity_lpm": [1.0, 2.0]}
    path = tmp_path / "smoke.yaml"
    path.write_text(yaml.safe_dump(cfg))
    study = run_study(path, tmp_path / "out", tmp_path / "runs")
    for name in ("pareto.png", "parity.png", "fields_recommended.png", "report.md", "results.json"):
        assert (tmp_path / "out" / name).stat().st_size > 0
    assert study["confidence"]["all_pass"]
    assert study["parity"]["max_abs_pct"] < 10
    results = json.loads((tmp_path / "out" / "results.json").read_text())
    rec = results["recommended"]
    assert (tmp_path / "runs" / rec / "results.json").exists()
    assert all(d["cfd"]["dp_pa"] <= cfg["requirements"]["dp_budget_kpa"] * 1e3 for d in results["designs"])
    assert results["rule"] in ("marginal_return", "fallback_lowest_r_th")
    assert results["r_th_threshold_k_per_w"] == pytest.approx(cfg["requirements"]["r_th_threshold_mk_w"] / 1e3)
    assert results["marginal_return_cutoff"] == cfg["requirements"]["marginal_return_cutoff"]
    if results["rule"] == "marginal_return":
        assert results["selection_steps"][-1]["from"] == rec and results["selection_steps"][-1]["to"] is None
    for d in results["designs"]:
        r = d["ratios"]
        assert r["capacity_kw"] == pytest.approx(results["r_th_threshold_k_per_w"] / d["cfd"]["r_th_k_per_w"])
        assert r["kw_per_pumping_w"] == pytest.approx(1.0 / (d["cfd"]["dp_pa"] * 1.5e-3 / 60))
    assert [r["flow_lpm"] for r in results["flow_sensitivity"]] == [1.0, 1.5, 2.0]
    assert f"`{results['rule']}`" in (tmp_path / "out" / "report.md").read_text()
    before = (tmp_path / "out" / "results.json").read_bytes()
    for name in ("pareto.png", "parity.png", "fields_recommended.png"):
        (tmp_path / "out" / name).unlink()
    plot_only(path, tmp_path / "out", tmp_path / "runs")
    for name in ("pareto.png", "parity.png", "fields_recommended.png"):
        assert (tmp_path / "out" / name).stat().st_size > 0
    assert (tmp_path / "out" / "results.json").read_bytes() == before


def test_round_ticks():
    assert round_ticks(41.92, 54.93) == [42.5, 45.0, 47.5, 50.0, 52.5]
    assert round_ticks(0.0, 0.194, min_ticks=3) == pytest.approx([0.0, 0.05, 0.1, 0.15])
    assert all(40.0 <= t <= 188.0 for t in round_ticks(40.0, 188.0))
