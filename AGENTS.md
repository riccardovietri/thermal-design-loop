# AGENTS.md

Instructions for any coding agent (or person) picking this project up cold.

## Read in this order
0. `docs/work-plan.md` — scoped instructions and order for the open issues.
1. `PROGRESS.md` — "Resume here" says what to do next. Update it at the end of every working session.
2. `SPEC.md` — frozen scope. Do not change it without Riccardo's OK; propose scope changes in one line with the reason.
3. `PLAN.md` — milestones and acceptance checks.
4. `ASSUMPTIONS.md` — every external number and modelling choice, with its basis. Add a row before relying on a new one.
5. `DECISIONS.md` — one line per decision: date, decision, why, alternatives rejected.
6. `docs/literature.md` — literature survey with confidence flags; snippet-level numbers are unverified.
7. `docs/build-log.md` — what was tried and what the numbers were; raw material for the write-up.

## Commands
- `uv sync` then `uv run pytest` (must pass before any commit).
- `uv run python run_study.py examples/gpu_coldplate.yaml --out examples/output/gpu_coldplate` — the whole Level 1 study, about 5 s (regenerates the committed example output).
  `--plot-only` redraws its three figures from `results.json` and `runs/` without CFD.
- `uv run python validation/duct_fully_developed.py` — Level 1 solver validation against Shah & London.
- `uv run python validation/level2_thermal_entrance.py` — Level 2 solver validation (developing flow between plates), about 20 s.
- `uv run python -m level2.slot_jet --slot_mm 0.3 --gap_mm 0.6 --pitch_mm 6 --flow_lpm 1.5` — one slot-jet case with a field plot
  (`--cells 25 --out docs/figures/slot_jet_example` reproduces the docs example exactly).
- `uv run python validation/level2_planview_straight.py` — plan-view gates G1–G3 on straight fins, about 3 s.
- `uv run python -m level2.planview_sweep` — plan-view checks (mesh, periodicity, developing strip, field figures, about 6 min) and the
  78-design topology sweep with the comparison plot and `summary.md` (about 5 min); `--checks`, `--sweep`, `--plot-only` run parts;
  `--extend` solves only designs not yet in `results.json`; `--fields` re-solves only the six field-figure designs (about 20 s).
- `uv run python -m level2.sweep --plot-only` — redraw the slot-jet comparison from saved results.
- `uv run python -m tools.screen_1d --w_ch_mm 0.42 --t_fin_mm 0.44 --h_fin_mm 4 --flow_lpm 1.5 --entrance` — one design with
  developing flow and inlet/outlet losses (the option is off by default; `analyze_1d(..., entrance_effects=True)`).
- `uv run python -m level2.entrance` — entrance effects on every CFD design in `examples/output/gpu_coldplate/results.json`
  and the resolved tall-fin check on five picked from it (widest, recommended and neighbours, narrowest); writes
  `examples/output/entrance_effects/` (about 5 min). Rerun after the Level 1 results change; `--level1`, `--out` select
  other files; `--plot-only` redraws the figure and summary.
- Level 3 (3D conjugate, issue #10 stage one) needs Docker: `docker pull opencfd/openfoam-default:2406` (about 1.9 GB, native
  arm64; the code pins it by digest, `level3.foam.IMAGE`). `uv run python -m level3.straight` solves the Level 1 recommended and
  widest-channel designs at three mesh levels in parallel (about 10 min, six one-core containers, under 0.7 GB each) and writes
  `runs/level3/<design_id>/level<k>/` and `examples/output/level3_straight_validation/`; `--design <id> --levels 1 2` selects cases,
  `--summary-only` rewrites the summary from `runs/level3/`. `uv run python validation/level3_straight.py` checks gates G1–G4 on
  the committed results in seconds (`--run` re-solves first).
- Level 3 stage two (circular pins against straight fins, same Docker image): `uv run python -m level3.pins_study --p0`
  (plan-view cross-check, about 2 min), `--parallel-check` (coarse pick serial against four processes, about 5 min),
  `--steady` (the three pin designs on three processes each plus the plan-view straight pick, about 70 min),
  `--summary-only` (rewrite `examples/output/level3_pins_vs_straight/` from `runs/level3/`, seconds), `--launch-long`
  (detached mesh-gate and transient jobs in containers `l3pins_mesh_gate` and `l3pins_transient`, kept awake with
  `caffeinate -i`; hours). `uv run python -m level3.status` says whether each detached job is running, finished or failed
  and prints its latest R_th and ΔP; `--summary-only` folds finished jobs into the summary.
  `uv run python validation/level3_pins.py` checks the recorded stage-two gates in seconds.
- Figures share one style, `tools/plotstyle.py` (8 in wide, 200 dpi, text ≥ 10 pt so it stays legible at 700 px). Titles are one
  short line; conditions and caveats go in the markdown caption next to the figure.
- `uv run python -m level2.sweep` — slot-jet mesh study, the 18-case sweep at rated flow and the comparison plot with fins, 10–20 minutes.

## Rules
- Every tool in `tools/` is a plain function with typed, JSON-serialisable inputs and outputs. No LLM calls in the pipeline.
- A physics change needs: a row in `ASSUMPTIONS.md` or `DECISIONS.md`, and a test against a hand calculation or a published value.
- Do not weaken a test or a confidence-check limit to make a run pass.
- Every CFD number must trace to `runs/<study>/<design_id>/`.
- Work on a branch; merge to `main` by pull request once tests pass. Conventional commit prefixes (`feat:`, `fix:`, `docs:`, `test:`, `examples:`).
- Do not edit `../thermal-mcp-server` or `../coldplate-opt`. Do not write to Notion.
- `prototypes/` is exploratory code: not covered by tests, not part of the pipeline.
- Test first for pure functions, physics formulas and selection logic: write the failing test from a hand calculation or a
  published value before the implementation. Field solves are validated against known solutions instead.
- External review (for example Codex, read-only) is optional input, not a merge gate. Use it when Riccardo asks for it, or for a
  change to physics, validation or selection logic where a second opinion is worth the time. Triage the findings, fix what is
  real, and stop after one round unless something major turns up.
- Text that quotes results (README.md, docs/index.html) must pass `tests/test_docs_consistency.py`; add a check there when
  a new number is quoted.
