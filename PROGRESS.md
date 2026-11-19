# PROGRESS

## Resume here
SPEC Amendment 4 (#17) is merged; Riccardo accepted the new Level 1 pick. README and the case study (`docs/index.html`, `docs/img/`)
match the committed results as of `docs/case-study-update`. Headline: Level 1 recommends 0.39 / 0.39 / 3.75 mm straight fins
(14.05 mK/W at 1.60 kPa, 6.3% margin on the 15 mK/W objective; 4,485 designs screened, 18 solved, 1D vs CFD 1.0% mean, 2.9% max);
across plan-view topologies the rule picks circular pins d 0.4, S_T 0.8, S_L 0.9 mm (12.65 mK/W at 1.42 kPa), 2.6–5.8% below
straight fins at equal ΔP, inside the plan view's tall-fin bias. No slot-jet case meets the 20 mK/W threshold.
Entrance effects (issue #12, pull request #18, open): with the re-anchored correction the pick reads 13.66 mK/W at 1.617 kPa and the
tall-fin gate passes (+3.1% R_th, −2.9% ΔP); the case study already quotes these numbers, so merge #18 with or before this branch.
Next: (1) Riccardo's call on the open questions below; (2) issue #10 (3D straight fins against pins), gate values from results.json;
(3) issue #11 (measured plate); (4) blog post and LinkedIn post from the case study.

## Open questions for Riccardo
- Entrance effects (issue #12): with the re-anchored correction applied, the selection rule would stop at 0.42 / 0.39 / 3.75 mm
  (14.39 mK/W at 1.31 kPa) instead of 0.39 / 0.39 / 3.75 mm (return 0.248 against 0.25). Keep fully developed values for selection?
- Plan view: the pin advantage (2.6–5.8% at equal ΔP) is smaller than the plan view's tall-fin bias measured on straight fins (−5.7% on
  R_th, no floor or cover). Treat it as unconfirmed until a 3D or end-wall check?
- Pins at Re_d,max 50–136 may shed vortices; steady symmetric solutions were computed. Accept, or add an unsteady check?
- Cross-topology pick (pins 0.8 × 0.9 mm) holds only at a cut-off of exactly 0.25 among the swept values; circular pins win at every
  cut-off. Is the 0.25 cut-off still the right call, given that it now decides which pin array is named?
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
- 2026-11-15: Entrance effects (issue #12, `feat/entrance-effects`): 1D option `entrance_effects=True` (off by default),
  `level2/entrance.py`, `examples/output/entrance_effects/summary.md`. Recommended design 14.93 → 14.44 mK/W, core ΔP
  1.294 → 1.320 kPa; objective still met; recommendation unchanged. Resolved-check gate on the 1D increments fails
  (−36% R_th, +41% ΔP, tall-fin limit), reported untuned. `uv run pytest` → 64 passed.
- 2026-11-15: Entrance correction re-anchored (approved) and rerun on the 18-design Level 1 set: thermal blend to Nu_H1,
  K(∞) = 0.680 (Ray, Durst & Ray 2017), inlet/outlet separate. Recommended 0.39 / 0.39 / 3.75 mm: 14.05 → 13.66 mK/W, core ΔP
  1.600 → 1.617 kPa, objective margin 8.9%. Gate now passes (+3.1% R_th, −2.9% ΔP). On corrected values the rule would stop at
  0.42 / 0.39 / 3.75 mm (return 0.248 vs 0.25); not applied. `uv run pytest` → 81 passed.

- 2026-11-15: SPEC Amendment 4 (issue #13): `tools.pareto.recommend_marginal_return_objective` (objective clause, ΔP-ratio guard
  1.10), denser Level 1 study (4,485 designs, all 18 front designs solved, about 5 s), cut-off sensitivity and local-return figures;
  plan view re-selected from saved results (offset strips now 14.58 mK/W at 1.57 kPa via the clause). Details in `docs/build-log.md`.
- 2026-11-15: README and case study brought in line with Amendment 4 (`docs/case-study-update`): fin-shape comparison as the first
  README figure, every case-study number and figure from the committed results, cut-off and local-return paragraph, entrance-effect
  paragraph (numbers from #18), indicative measured-plate comparison (issue #11), work-plan gate for issue #10 reads results.json.
- 2026-11-18: Review corrections (`docs/review-corrections`), after an external adversarial review (Codex, GPT-6) of main at
  d85f40c: `tests/test_docs_consistency.py` (README and case study against the committed results, written failing first);
  case study and README corrected (cross-topology pick 12.6 not 12.7 mK/W, 1D split labelled as 1D with its 13.9 total,
  wetted area 6.2 to 10.3 times the heated footprint, ratio wording, entrance gate scope, energy-balance and duct-validation
  scope, heated width 59.64 mm); report, plan-view and entrance generators reworded and regenerated without solves. No
  results JSON changed.

## Result of the example study
Rated flow 1.5 LPM, 4,485 designs screened (1,783 feasible within 5 kPa), all 18 on the 1D front solved. 1D vs CFD on R_th: mean
|diff| 1.0%, max 2.9%. All confidence checks pass.
Recommended (`marginal_return`, cut-off 0.25, ΔP-ratio guard 1.10): 0.39 mm channels, 0.39 mm fins, 3.75 mm tall → 14.05 mK/W at
1.60 kPa, meets the objective (margin 6.3%), interface 54.0 °C average, 59.0 °C peak. Walk: 0.55 / 0.44 / 4.0 (18.47, 0.70 kPa) →
0.50 (0.378) → 0.46 (0.346) → 0.42 / 0.39 / 3.75 (0.327) → 0.39 / 0.39 / 3.75 (0.277) → stop (best further 0.187).
Pick holds for cut-offs 0.20–0.25; 0.30–0.60 → 0.42 / 0.39 / 3.75 mm (14.90 mK/W at 1.29 kPa); 0.15 → 2.64 kPa; ≤ 0.10 → 4.66 kPa.
Flow sensitivity: 1.0 LPM 16.55 mK/W (misses the objective), 2.0 LPM 12.80 mK/W.
Lowest ΔP meeting the threshold: 0.55 / 0.44 / 4.0 mm → 18.47 mK/W at 0.70 kPa. Lowest R_th in budget: 0.27 / 0.26 / 2.5 mm → 11.92 mK/W at 4.66 kPa.
Slot jets at 1.5 LPM with ≥ 0.1 mm features: best 23.8 mK/W at 0.49 kPa; none meets the 20 mK/W threshold.

## Open items
- Code findings of the 2026-11-18 review, each with a failing test to write first: #20 (entrance option: peak temperature
  uses length-mean h), #21 (failed confidence checks must block recommendations), #22 (selection guard depends on front
  sampling), #23 (plan-view reuse ignores operating conditions and solver settings), #24 (validate the production
  conjugate solver through `run_cfd()`), #25 (input validation and laminar limit at the solver boundary).
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
- GitHub: public repo riccardovietri/thermal-design-loop (renamed from coldplate-design-loop on 2026-11-12); work goes on branches and merges by pull request.
