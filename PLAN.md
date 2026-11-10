# PLAN

| # | Day | Milestone | Acceptance check | Status |
|---|---|---|---|---|
| M0 | Tue 11/10 | Spec, toolchain decision, validation case | `validation/duct_fully_developed.py` within 1% of Shah & London; SPEC approved | done 11/10 |
| M1 | Wed 11/11 | One design end to end: geometry → mesh → solve → R_th, ΔP | Energy-balance error < 0.5%; limiting case (k_s → ∞, tall thin channel) reproduces duct validation; result within plausible range of 1D | done 11/10 (f·Re exact; parallel-plate Nu limit within 5%) |
| M2 | Wed/Thu | 1D straight-fin screen | Tests pass: polynomial vs Shah & London table, fin-efficiency hand calc, caloric term | |
| M3 | Thu 11/12 | Batch over ~10 Pareto-spread designs, `runs/<id>/` traceability | All runs complete from one command; each run dir has inputs, log, results JSON | |
| M4 | Fri 11/13 | Pareto + parity plots, confidence checks, `report.md` | `run_study.py examples/gpu_coldplate.yaml` produces all three; mesh-independence change < 1% | |
| M5 | Sat 11/14 | README, hero figure, publish | Fresh clone + `uv sync` + one command reproduces the report | |
| M6 | Sun 11/15 | Buffer; LinkedIn post | | |

Gate: if M1 is not working Wednesday night, stop and propose a simpler path.
Stretch (only after M4): MCP agent wrapper; streamwise developing-flow case.
