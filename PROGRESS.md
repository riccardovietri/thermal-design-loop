# PROGRESS

## Resume here
Level 1 is merged to `main` (pull request #1). Level 2 (planar Navier-Stokes, slot-jet arrays) is on `feat/level2-planar`
as pull request #2: solver validated, one example case with a field plot, a 24-case sweep and a comparison plot against fins.
Next: (1) Riccardo answers the open questions below; (2) blog post on the fidelity ladder (levels 1 and 2), format modelled on
Anthropic/OpenAI engineering posts, linking to the reports; (3) final repo name, make public; (4) LinkedIn post.

## Open questions for Riccardo
- ΔP budget: change 30 kPa to 10 kPa core-only? (ASSUMPTIONS.md R5)
- Minimum manufacturable jet gap? Competitive slot-jet designs need 0.05–0.1 mm. (ASSUMPTIONS.md L8)
- Compare jets and fins at equal flow rate, and extend the jet sweep to pitches under 2 mm?

## Done
- 2026-11-10: M0 spec frozen, toolchain validated (0.007% vs Shah & London).
- 2026-11-10: M1 one design end to end (`tools/mesh.py`, `tools/run_cfd.py`).
- 2026-11-10: M2 1D straight-fin screen (`tools/screen_1d.py`), Pareto tools (`tools/pareto.py`).
- 2026-11-10: M3/M4 batch study, confidence checks, plots, report (`run_study.py`, `tools/check_confidence.py`,
  `tools/report.py`). `uv run pytest` → 16 passed.
- 2026-11-10: Level 2 solver (`level2/planar.py`), validation (`validation/level2_thermal_entrance.py`), slot-jet case and sweep
  (`level2/slot_jet.py`, `level2/sweep.py`). `uv run pytest` → 18 passed.
- 2026-11-10: Level 2 prototype (`prototypes/level2_wavy_channel.py`), findings in `docs/build-log.md`.

## Result of the example study
1008 designs screened (530 feasible), 10 solved. 1D vs CFD on R_th: mean |diff| 0.9%, max 2.9%. All confidence checks pass.
Knee design: 0.44 mm channels, 0.47 mm fins, 4.0 mm tall, 1.5 LPM → 15.7 mK/W at 1.2 kPa.
Lowest R_th in budget: 0.15 / 0.15 / 1.5 mm, 1.5 LPM → 9.3 mK/W at 26.4 kPa.

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
