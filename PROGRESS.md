# PROGRESS

## Resume here
Levels 1 and 2 are merged to `main` (pull requests #1 and #2). The requirements amendment (SPEC Amendment 1: 1.5 LPM rated flow,
6 kPa core budget, 15 mK/W target, target-based recommendation rule, 0.1 mm jet features, region-of-interest plots) is on
`feat/rated-flow-requirements` as a pull request awaiting Riccardo's review.
Next: (1) Riccardo reviews and merges the amendment; (2) blog post on the fidelity ladder (levels 1 and 2), format modelled on
Anthropic/OpenAI engineering posts, linking to the reports; (3) final repo name, make public; (4) LinkedIn post.

## Open questions for Riccardo
- Next topology: manifold microchannels (jets feeding short fin channels) via a depth-averaged 2D model, validated against Level 1
  in the long-channel limit? Supported by docs/literature.md; needs Riccardo's OK as new scope.
- The recommended design meets 15 mK/W with a 0.5% margin (14.93 mK/W) and misses it at 1.0 LPM (17.4 mK/W). Is a margin
  needed on the target (e.g. for entrance effects, manifolds or tolerance), which would move the pick to 0.36 mm channels (13.8 mK/W at 2.3 kPa)?
- Slot jets on a flat base cannot reach 15 mK/W at ≥ 0.1 mm gaps (best 23.8 mK/W). Drop them, or try jets onto a finned base?

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

## Result of the example study
Rated flow 1.5 LPM, 672 designs screened (273 feasible within 6 kPa), 9 solved. 1D vs CFD on R_th: mean |diff| 1.1%, max 2.9%.
All confidence checks pass.
Recommended (`meets_target_lowest_dp`): 0.42 mm channels, 0.44 mm fins, 4.0 mm tall → 14.9 mK/W at 1.29 kPa, interface 54.9 °C average.
Knee (secondary): 0.50 / 0.44 / 4.0 mm → 17.1 mK/W at 0.86 kPa. Lowest R_th in budget: 0.25 / 0.26 / 2.5 mm → 11.3 mK/W at 5.75 kPa.
Slot jets at 1.5 LPM with ≥ 0.1 mm features: best 23.8 mK/W at 0.49 kPa; none meets the target.

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
