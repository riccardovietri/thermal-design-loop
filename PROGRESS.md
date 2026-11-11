# PROGRESS

## Resume here
V1 (Level 1) is functionally complete on branch `feat/v1-study`:
`uv run python run_study.py examples/gpu_coldplate.yaml` → `out/gpu_coldplate/` (pareto.png, parity.png,
fields_recommended.png, report.md, results.json) in about 3 s. A committed copy is in `examples/output/gpu_coldplate/`.
Next: (1) Riccardo reviews the report and the recommendation rule; (2) README + hero figure (M5);
(3) Level 2: thermal validation case, then slot-jet impingement; (4) blog post on the fidelity ladder.

## Done
- 2026-11-10: M0 spec frozen, toolchain validated (0.007% vs Shah & London).
- 2026-11-10: M1 one design end to end (`tools/mesh.py`, `tools/run_cfd.py`).
- 2026-11-10: M2 1D straight-fin screen (`tools/screen_1d.py`), Pareto tools (`tools/pareto.py`).
- 2026-11-10: M3/M4 batch study, confidence checks, plots, report (`run_study.py`, `tools/check_confidence.py`,
  `tools/report.py`). `uv run pytest` → 16 passed.
- 2026-11-10: Level 2 prototype (`prototypes/level2_wavy_channel.py`), findings in `docs/build-log.md`.

## Result of the example study
1008 designs screened (530 feasible), 10 solved. 1D vs CFD on R_th: mean |diff| 0.9%, max 2.9%. All confidence checks pass.
Knee design: 0.44 mm channels, 0.47 mm fins, 4.0 mm tall, 1.5 LPM → 15.7 mK/W at 1.2 kPa.
Lowest R_th in budget: 0.15 / 0.15 / 1.5 mm, 1.5 LPM → 9.3 mK/W at 26.4 kPa.

## Open items
- Branches `feat/m1-single-design` and `feat/v1-study` are pushed, not merged to `main` (needs a PR; Riccardo to say).
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
