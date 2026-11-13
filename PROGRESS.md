# PROGRESS

## Resume here
SPEC Amendment 3 (threshold 20 / objective 15 mK/W, 5 kPa budget, marginal-return rule with cut-off 0.25, capacity ratios,
H / gap ≤ 10 in the plan view) is on `feat/threshold-objective-criteria` as a pull request awaiting Riccardo's review.
Headline: Level 1 still recommends 0.42 / 0.44 / 4.0 mm (14.93 mK/W at 1.29 kPa, meets the objective); across plan-view topologies the
rule picks circular pins d 0.4, S_T 0.8, S_L 0.9 mm (12.65 mK/W at 1.42 kPa), 2.6–5.8% below straight fins at equal ΔP, inside the
plan view's 6% tall-fin bias. No slot-jet case meets the 20 mK/W threshold.
Next: (1) Riccardo reviews the pull request and the open questions below; (2) blog post on the fidelity ladder;
(3) final repo name, make public; (4) LinkedIn post.

## Open questions for Riccardo
- Plan view: the pin advantage (2.6–5.8% at equal ΔP) is smaller than the plan view's tall-fin bias measured on straight fins (−5.7% on
  R_th, no floor or cover). Treat it as unconfirmed until a 3D or end-wall check?
- Pins at Re_d,max 50–136 may shed vortices; steady symmetric solutions were computed. Accept, or add an unsteady check?
- Marginal-return rule on the offset strips stops at the first eligible point (all returns 0.19–0.23 < 0.25), so the per-topology
  pick (18.4 mK/W) misses the objective while other offset designs meet it. Acceptable as stated, or should a pick that misses the
  objective walk on to the lowest-ΔP point meeting it?
- Non-straight planforms have no known volume process at 0.15 mm in copper (ASSUMPTIONS V11). Keep them as a modelling comparison only?
- Next topology: manifold microchannels via a depth-averaged 2D model? Needs Riccardo's OK as new scope.

## Done
- 2026-11-10: M0 spec frozen, toolchain validated (0.007% vs Shah & London).
- 2026-11-10: M1 one design end to end (`tools/mesh.py`, `tools/run_cfd.py`).
- 2026-11-10: M2 1D straight-fin screen (`tools/screen_1d.py`), Pareto tools (`tools/pareto.py`).
- 2026-11-10: M3/M4 batch study, confidence checks, plots, report (`run_study.py`, `tools/check_confidence.py`,
  `tools/report.py`). `uv run pytest` → 16 passed.
- 2026-11-10: Level 2 solver (`level2/planar.py`), validation (`validation/level2_thermal_entrance.py`), slot-jet case and sweep
  (`level2/slot_jet.py`, `level2/sweep.py`). `uv run pytest` → 18 passed.
- 2026-11-10: Level 2 prototype (`prototypes/level2_wavy_channel.py`), findings in `docs/build-log.md`.
- 2026-11-10: SPEC Amendment 1: rated flow, 6 kPa budget, 15 mK/W target, `tools.pareto.recommend`, flow sensitivity,
  region-of-interest plots, 18-case slot-jet sweep at rated flow with a mesh study. `uv run pytest` → 21 passed.

- 2026-11-11: SPEC Amendment 2, plan-view Level 2 for fin planforms (`level2/planforms.py`, `level2/planview.py`,
  `level2/planview_sweep.py`, `validation/level2_planview_straight.py`); gates G1–G3 pass; 63 designs, all converged.
  `uv run pytest` → 30 passed.
- 2026-11-11: SPEC Amendment 3: two-tier R_th (20 / 15 mK/W), 5 kPa budget, `tools.pareto.recommend_marginal_return` (cut-off 0.25)
  and `capacity_ratios`, new report columns and step table, plots with threshold and objective windows; plan view: H / gap ≤ 10
  applied (36 designs excluded), 15 designs added (78 converged), per-topology and cross-topology picks in
  `examples/output/level2_planview/summary.md`. `uv run pytest` → 38 passed.

- 2026-11-12: Readable figures (`feat/readable-figures`): shared style `tools/plotstyle.py`, one-line titles with captions in
  the markdown, plan-view comparison split into an overview (`pareto_topologies.png`) and per-planform panels
  (`pareto_topologies_panels.png`), `--plot-only` for `run_study.py`, `--fields` for the plan-view sweep, `--cells` for the slot-jet CLI.
  No numeric result changed.

## Result of the example study
Rated flow 1.5 LPM, 672 designs screened (258 feasible within 5 kPa), 8 solved. 1D vs CFD on R_th: mean |diff| 1.1%, max 2.9%.
All confidence checks pass.
Recommended (`marginal_return`, cut-off 0.25): 0.42 mm channels, 0.44 mm fins, 4.0 mm tall → 14.93 mK/W at 1.29 kPa, meets the
objective, interface 54.9 °C average. Steps: from 0.50 mm (17.12 mK/W, 0.86 kPa) return 0.331 → 0.42 mm; further returns 0.141 and 0.159 → stop.
Lowest ΔP meeting the threshold: 0.50 / 0.44 / 4.0 mm → 17.12 mK/W at 0.86 kPa. Lowest R_th in budget: 0.30 / 0.26 / 2.5 mm → 12.58 mK/W at 3.79 kPa.
Slot jets at 1.5 LPM with ≥ 0.1 mm features: best 23.8 mK/W at 0.49 kPa; none meets the 20 mK/W threshold.

## Open items
- Merge policy (Riccardo agreed 11/10): branch per milestone, pull request into `main`, merge when CI passes.
- README, ASSUMPTIONS.md, AGENTS.md and CI (`.github/workflows/ci.yml`) added 11/10.
- CLIs exist for screen_1d, mesh, postprocess and run_study; run_cfd, check_confidence, pareto and report are functions only.
- The energy-balance check passes at 1e-12 because the base temperature field is linear near the chip interface;
  it is a real check but a weak one.
- Dow lists cp about 4% above CoolProp for PG25; unresolved, low impact.
- Blog post format: Riccardo wants depth balanced with demo value, modelled on Anthropic/OpenAI engineering posts,
  linking to the detailed report. Draft after Level 2.

## Notes for a fresh session
- Run everything with `uv run`. In Claude Code's sandbox, `uv`, `git` and `gh` writes may need the sandbox disabled.
- Fable usage: 7% used after M1 on 11/10.
- GitHub: private repo riccardovietri/coldplate-design-loop; `main` pushed once with Riccardo's OK, work goes on branches.
