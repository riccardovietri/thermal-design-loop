# AGENTS.md

Instructions for any coding agent (or person) picking this project up cold.

## Read in this order
1. `PROGRESS.md` — "Resume here" says what to do next. Update it at the end of every working session.
2. `SPEC.md` — frozen scope. Do not change it without Riccardo's OK; propose scope changes in one line with the reason.
3. `PLAN.md` — milestones and acceptance checks.
4. `ASSUMPTIONS.md` — every external number and modelling choice, with its basis. Add a row before relying on a new one.
5. `DECISIONS.md` — one line per decision: date, decision, why, alternatives rejected.
6. `docs/build-log.md` — what was tried and what the numbers were; raw material for the write-up.

## Commands
- `uv sync` then `uv run pytest` (must pass before any commit).
- `uv run python run_study.py examples/gpu_coldplate.yaml` — the whole Level 1 study, about 3 s.
- `uv run python validation/duct_fully_developed.py` — solver validation against Shah & London.

## Rules
- Every tool in `tools/` is a plain function with typed, JSON-serialisable inputs and outputs. No LLM calls in the pipeline.
- A physics change needs: a row in `ASSUMPTIONS.md` or `DECISIONS.md`, and a test against a hand calculation or a published value.
- Do not weaken a test or a confidence-check limit to make a run pass.
- Every CFD number must trace to `runs/<study>/<design_id>/`.
- Work on a branch; merge to `main` by pull request once tests pass. Conventional commit prefixes (`feat:`, `fix:`, `docs:`, `test:`, `examples:`).
- Do not edit `../thermal-mcp-server` or `../coldplate-opt`. Do not write to Notion.
- `prototypes/` is exploratory code: not covered by tests, not part of the pipeline.
