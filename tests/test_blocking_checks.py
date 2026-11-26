"""Issue #21: a failed confidence check must block the recommendation and give a non-zero exit status.

Failures are injected through monkeypatching (a solver result with a bad residual or energy balance, a poor mesh, a mesh
study that moves), never by loosening a limit.
"""
import json
import shutil
import sys
from dataclasses import replace
from pathlib import Path

import pytest
import yaml

import run_study as rs
from tools.check_confidence import LIMITS, run_failures
from tools.types import CfdResult

GOOD = CfdResult(design_id="d", n_channels=77, reynolds=150.0, dp_pa=1600.0, pump_power_w=0.04, f_re_fanning=20.0, nu_eff=5.0,
                 r_th_k_per_w=0.014, r_th_max_k_per_w=0.019, base_avg_temp_c=54.0, base_max_temp_c=59.0, coolant_rise_k=10.0,
                 energy_residual_rel=1e-14, n_dofs=1000, linear_residual_rel=1e-12, base_flux_error_rel=1e-12)


def test_run_failures_passing_run_has_none():
    assert run_failures(GOOD, 0.9) == []


@pytest.mark.parametrize("change, quality, word", [
    ({"linear_residual_rel": 1e-3}, 0.9, "residual"),
    ({"base_flux_error_rel": 0.05}, 0.9, "energy"),
    ({"base_flux_error_rel": -0.05}, 0.9, "energy"),
    ({"reynolds": 2500.0}, 0.9, "Reynolds"),
    ({}, 0.1, "mesh quality"),
])
def test_run_failures_names_each_failed_check(change, quality, word):
    reasons = run_failures(replace(GOOD, **change), quality)
    assert len(reasons) == 1 and word in reasons[0]


def test_run_failures_uses_the_existing_limits():
    assert run_failures(replace(GOOD, linear_residual_rel=0.5 * LIMITS["linear_residual_rel"]), 0.9) == []
    assert run_failures(replace(GOOD, linear_residual_rel=2 * LIMITS["linear_residual_rel"]), 0.9) != []


def _config(tmp_path: Path) -> Path:
    cfg = yaml.safe_load(open("examples/gpu_coldplate.yaml"))
    cfg["name"] = "smoke"
    cfg["design_space"].update({"w_ch_mm": [0.3, 0.6], "t_fin_mm": [0.3], "h_fin_mm": [2.0, 3.0], "flow_lpm": [1.5]})
    cfg["study"] = {"n_cfd": 3, "mesh_cells_across_half_channel": 4, "flow_sensitivity_lpm": [1.0, 2.0]}
    path = tmp_path / "smoke.yaml"
    path.write_text(yaml.safe_dump(cfg))
    return path


@pytest.fixture(scope="module")
def passing(tmp_path_factory):
    """The smoke study with nothing injected: its recommendation, report and results."""
    tmp = tmp_path_factory.mktemp("passing")
    cfg = _config(tmp)
    study = rs.run_study(cfg, tmp / "out", tmp / "runs")
    return {"cfg": cfg, "study": study, "rec": study["recommended"][1].design_id,
            "results": json.loads((tmp / "out" / "results.json").read_text()), "report": (tmp / "out" / "report.md").read_text()}


def test_passing_study_records_no_block(passing):
    assert passing["study"]["confidence"]["all_pass"]
    assert "status" not in passing["results"] and "excluded_designs" not in passing["results"]
    assert passing["report"].startswith("# Cold-plate design study: smoke\n\n**Recommended design:**")
    assert "checks failed" not in passing["report"].lower()


def _inject_run(monkeypatch, bad_id: str, **change):
    real = rs.run_cfd

    def fake(design, *args, **kw):
        r = real(design, *args, **kw)
        return replace(r, **change) if r.design_id == bad_id else r
    monkeypatch.setattr(rs, "run_cfd", fake)


def _blocked(tmp_path, cfg):
    study = rs.run_study(cfg, tmp_path / "out", tmp_path / "runs")
    return study, json.loads((tmp_path / "out" / "results.json").read_text()), (tmp_path / "out" / "report.md").read_text()


@pytest.mark.parametrize("change, word", [({"linear_residual_rel": 1e-3}, "residual"), ({"base_flux_error_rel": 0.05}, "energy")])
def test_failed_run_check_excludes_design_and_blocks(passing, tmp_path, monkeypatch, change, word):
    bad = passing["rec"]
    _inject_run(monkeypatch, bad, **change)
    study, res, report = _blocked(tmp_path, passing["cfg"])
    assert res["status"] == "blocked"
    assert any(word in c.lower() for c in res["failed_checks"])
    assert [e["design_id"] for e in res["excluded_designs"]] == [bad]
    assert word in res["excluded_designs"][0]["reasons"][0].lower()
    assert res["recommended"] != bad and bad not in [d["cfd"]["design_id"] for d in res["designs"]]
    head = report.split("\n\n")[1]
    assert head.startswith("**Confidence checks failed") and word in head.lower()
    assert "**Recommended design:**" not in report
    assert f"`{bad}`" in report.split("## Excluded designs")[1]


def test_poor_mesh_quality_excludes_design(passing, tmp_path, monkeypatch):
    bad = passing["rec"]
    real = rs.mesh

    def fake(design, *args, **kw):
        info = real(design, *args, **kw)
        return replace(info, min_quality_sicn=0.1) if design.design_id == bad else info
    monkeypatch.setattr(rs, "mesh", fake)
    _, res, _ = _blocked(tmp_path, passing["cfg"])
    assert res["status"] == "blocked" and any("quality" in c.lower() for c in res["failed_checks"])
    assert [e["design_id"] for e in res["excluded_designs"]] == [bad]


def test_failed_mesh_study_blocks_without_excluding(passing, tmp_path, monkeypatch):
    real = rs.mesh_independence

    def fake(*args, **kw):
        return {**real(*args, **kw), "r_th_change_pct": 5.0}
    monkeypatch.setattr(rs, "mesh_independence", fake)
    _, res, report = _blocked(tmp_path, passing["cfg"])
    assert res["status"] == "blocked" and res["excluded_designs"] == []
    assert any("Mesh independence, R_th" in c for c in res["failed_checks"])
    assert res["recommended"] == passing["rec"]  # the pick is shown, but only as provisional
    assert report.split("\n\n")[1].startswith("**Confidence checks failed")


def test_every_design_failing_stops_the_study(passing, tmp_path, monkeypatch):
    real = rs.run_cfd
    monkeypatch.setattr(rs, "run_cfd", lambda *a, **k: replace(real(*a, **k), linear_residual_rel=1.0))
    with pytest.raises(SystemExit) as e:
        rs.run_study(passing["cfg"], tmp_path / "out", tmp_path / "runs")
    assert "check" in str(e.value.code).lower()


def _main(monkeypatch, tmp_path, cfg) -> int:
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["run_study.py", str(cfg), "--out", str(tmp_path / "out")])
    return rs.main()


def test_main_exit_status(passing, tmp_path, monkeypatch):
    for sub in ("a", "b"):
        (tmp_path / sub).mkdir()
    assert _main(monkeypatch, tmp_path / "a", passing["cfg"]) == 0
    _inject_run(monkeypatch, passing["rec"], base_flux_error_rel=0.05)
    assert _main(monkeypatch, tmp_path / "b", passing["cfg"]) != 0


# --- level2/entrance.py: a failed gate gives a non-zero exit status ---

def _entrance_copy(tmp_path: Path) -> Path:
    out = tmp_path / "entrance"
    out.mkdir()
    shutil.copy("examples/output/entrance_effects/results.json", out / "results.json")
    return out


def test_entrance_committed_results_exit_zero(tmp_path):
    from level2.entrance import run
    out = _entrance_copy(tmp_path)
    assert run(plot_only=True, out_dir=out) == 0
    data = json.loads((out / "results.json").read_text())
    outside = [r for r in data["resolved"] if abs(r["comparison"]["r_th_contribution_disagreement"]) > 0.10]
    assert outside and all(r["design_id"] != data["level1_recommended"] for r in outside)  # reported, not gated


@pytest.mark.parametrize("field", ["energy", "agreement"])
def test_entrance_failed_gate_exits_nonzero(tmp_path, field):
    from level2.entrance import run
    out = _entrance_copy(tmp_path)
    data = json.loads((out / "results.json").read_text())
    rec = next(r for r in data["resolved"] if r["design_id"] == data["level1_recommended"])
    if field == "energy":
        rec["uniform"]["energy_balance_error"] = 0.05
    else:
        rec["comparison"]["r_th_contribution_disagreement"] = 0.25
    (out / "results.json").write_text(json.dumps(data))
    assert run(plot_only=True, out_dir=out) != 0


# --- level2/planview_sweep.py: failed energy balance or no convergence excludes a design ---

def _planview():
    return json.loads(Path("examples/output/level2_planview/results.json").read_text())


def test_planview_committed_picks_unchanged():
    from level2.planview_sweep import summarise
    data = _planview()
    summ = json.loads(json.dumps(summarise(data), default=float))  # as written to results.json
    assert summ == data["summary"]
    assert "excluded_checks" not in summ


def test_planview_energy_failure_and_non_convergence_excluded():
    from level2.planview_sweep import ENERGY_LIMIT, summarise
    data = _planview()
    pick = data["summary"]["recommended"]["design_id"]
    rows = data["designs"]
    k = next(i for i, r in enumerate(rows) if r["design_id"] == pick)
    rows[k]["result"]["energy_balance_error"] = 2 * ENERGY_LIMIT
    j = next(i for i, r in enumerate(rows) if r["topology"] == "straight" and i != k)
    nc = rows[j]["design_id"]
    rows[j] = {key: rows[j][key] for key in ("topology", "params_mm", "flow_lpm", "design_id")} | {
        "status": "newton_failed", "failed_at_reynolds": 100.0}
    summ = summarise(data)
    reasons = {e["design_id"]: e["reason"] for e in summ["excluded_checks"]}
    assert "energy balance" in reasons[pick] and "not converged" in reasons[nc]
    assert summ["recommended"]["design_id"] != pick
    assert pick not in summ["per_topology"]["pin_circle"]["recommended"]["front"]
    assert all(v["design_id"] != pick for v in summ["vs_straight_front"])
