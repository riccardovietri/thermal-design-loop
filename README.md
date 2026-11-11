# Cold-plate design loop

An engineer's design loop for a liquid cold plate, automated end to end: a 1D physics model screens a thousand designs,
a CFD model checks the ten that matter, and the output is a Pareto front, a recommendation and the evidence for trusting it.

![Thermal resistance vs pressure drop](examples/output/gpu_coldplate/pareto.png)

One command, about three seconds on a laptop, no LLM in the loop:

```bash
uv sync
uv run python run_study.py examples/gpu_coldplate.yaml
```

## The problem

A 1 kW chip on a 60 × 60 mm copper cold plate, cooled by 25% propylene glycol at 40 °C. At a rated flow of 1.5 LPM, with a channel-core
pressure budget of 6 kPa and a target of 15 mK/W, which fin geometry meets the target for the least pressure drop, and how far can a hand-calculation-level model be trusted to answer that?

## How it works

```mermaid
flowchart LR
    A[Requirements<br/>requirements.yaml] --> B[1D screen<br/>672 designs, 2 ms]
    B --> C[CFD on the shortlist<br/>up to 10 designs along the 1D Pareto front]
    C --> D[Confidence checks<br/>mesh, residuals, energy balance]
    D --> E[Pareto front, parity plot,<br/>recommendation, report.md]
```

1. **Screen.** A 1D resistance network (Shah & London laminar duct correlations, fin efficiency, coolant heating) evaluates every
   combination of channel width, fin thickness and fin height at rated flow, and discards designs that break the pressure or manufacturing limits.
2. **Verify.** Up to ten designs spread along the 1D Pareto front are solved with a 2D conjugate finite-element model of the channel cross-section
   (gmsh + scikit-fem): fluid velocity, and temperature in the fluid, fin and base together.
3. **Decide.** Mesh independence, solver residual and energy balance are checked and reported. Among CFD-verified designs that meet both the pressure budget and the
   thermal target, the one with the lowest pressure drop is recommended; if none meets the target, the lowest resistance within budget.
   The knee of the front is reported as secondary information.

## What it found

| | Channel / fin / height | R_th | ΔP |
|---|---|---|---|
| Recommended (lowest ΔP meeting 15 mK/W) | 0.42 / 0.44 / 4.0 mm | 14.9 mK/W | 1.29 kPa |
| Knee of the CFD front (secondary) | 0.50 / 0.44 / 4.0 mm | 17.1 mK/W | 0.86 kPa |
| Lowest resistance in budget | 0.25 / 0.26 / 2.5 mm | 11.3 mK/W | 5.75 kPa |

All at the 1.5 LPM rated flow. The recommended design averages 54.9 °C at the chip interface; at 1.0 LPM it rises to 17.4 mK/W
(57.4 °C), at 2.0 LPM it falls to 13.7 mK/W (53.7 °C).
The 1D screen and the CFD agree on thermal resistance to 1.1% on average and 2.9% at worst.
For straight fins in laminar flow, a correct 1D model is enough; the CFD is the evidence for that statement, not a replacement for it.

![1D vs CFD](examples/output/gpu_coldplate/parity.png)

![Velocity and temperature in the recommended design](examples/output/gpu_coldplate/fields_recommended.png)

Full output: [`examples/output/gpu_coldplate/report.md`](examples/output/gpu_coldplate/report.md).

Level 2 (`level2/`, 2D planar Navier-Stokes) tests confined slot-jet arrays on a flat base at the same flow, with jet gaps and slots
of at least 0.1 mm: the best reaches 23.8 mK/W at 0.49 kPa, so none meets the target
([comparison](examples/output/level2_slot_jet/pareto_compare.png)).

## Why trust the solver

`uv run python validation/duct_fully_developed.py` solves laminar flow and heat transfer in rectangular ducts at four aspect ratios
and reproduces the Shah & London reference values for f·Re and Nu to within 0.007%.

## What it does not do

The CFD is a 2D cross-section with fully developed flow. It does not capture entrance effects, manifolds, flow maldistribution,
or any geometry that varies along the flow. Every assumption is listed with its basis in [`ASSUMPTIONS.md`](ASSUMPTIONS.md).

## Layout

- `run_study.py` — the deterministic orchestrator.
- `tools/` — one function per engineering step: `screen_1d`, `mesh`, `run_cfd`, `postprocess`, `check_confidence`, `pareto`, `report`.
- `examples/` — the requirements file and a committed copy of its output.
- `validation/` — solver validation against analytical results.
- `prototypes/` — exploratory work toward the next fidelity level (2D Navier-Stokes).
- `SPEC.md`, `PLAN.md`, `DECISIONS.md`, `ASSUMPTIONS.md`, `PROGRESS.md` — scope, milestones, decisions, assumptions and current state.

MIT licence.
