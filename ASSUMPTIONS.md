# ASSUMPTIONS

Every number or modelling choice that is not derived inside this repo. One row each, so nothing is implicit.
Status: **sourced** (citation), **chosen** (Riccardo's requirement), **assumed** (judgment, no source), **open** (needs a decision).
Change a row and DECISIONS.md together.

## Requirements
| ID | Assumption | Value | Status | Basis | Effect if wrong |
|---|---|---|---|---|---|
| R1 | Heat load | 1,000 W | chosen | Riccardo, placeholder for a GPU-class device | R_th unchanged; temperatures scale linearly |
| R2 | Heat flux distribution | Uniform over 60 × 60 mm | assumed | Simplest case | A smaller die adds spreading resistance in the base; thicker base would then help |
| R3 | Coolant and inlet | PG25, 40 °C | chosen | Riccardo; data-centre practice | Water would cut ΔP roughly in half and lower R_th |
| R4 | Flow rate | 1–2 LPM | chosen / sourced | Industry guidance of 1.2–1.5 LPM per kW (Electronics Cooling 2025; CoolIT) | Sets the coolant-heating floor on R_th |
| R5 | ΔP budget, channel core | 30 kPa | **open** | Riccardo, "somewhat arbitrary". Published context: 124 kPa rack budget for 80 kW with 34 kPa in hoses and fittings; a 4 kW plate at 55 kPa for the full loop with fittings; a measured plate at 4 kPa. A core-only budget of 5–10 kPa looks more representative | Moves the right-hand end of the Pareto front; does not move the knee |
| R6 | Plate material | Copper, k = 390 W/m·K | assumed | Typical for C110 | ~1% on R_th |
| R7 | Base thickness | 2 mm | assumed | Literature spans 0.1–4.5 mm, no typical value | Fixed 1.4 mK/W offset; ranking unchanged |

## Design space and manufacturing
| ID | Assumption | Value | Status | Basis | Effect if wrong |
|---|---|---|---|---|---|
| M1 | Minimum fin and channel width | 0.15 mm | sourced | Boyd: copper skived fins 0.1–0.8 mm, 0.1 mm minimum gap; ASME study: 100–200 µm | Lower limit moves the low-R_th end of the front |
| M2 | Maximum aspect ratio (height / fin, height / channel) | 10 | sourced, conservative | ASME study reports 10–20 | Higher limit allows taller fins at fine pitch |
| M3 | Fin height range | 1–4 mm | assumed | Fin efficiency falls off above ~3 mm for thin fins | Front hits the 4 mm bound at the low-ΔP end; see report |
| M4 | Topology | Straight rectangular fins, parallel channels, no tip clearance | chosen | V1 scope | Other topologies need Level 2 or 3 |

## Physics and numerics
| ID | Assumption | Value | Status | Basis | Effect if wrong |
|---|---|---|---|---|---|
| P1 | Flow regime | Laminar, steady | checked per run | Re 80–170 on the shortlist; runs with Re > 2000 are rejected | — |
| P2 | Fully developed flow and heat transfer | Both 1D and CFD | assumed | x* = L/(Dh·Re·Pr) is 0.02–0.2 across the shortlist, so the thermal entrance region is not negligible for the widest channels | Real R_th lower and ΔP higher than reported; Level 2 quantifies it |
| P3 | Fluid properties | Constant, at inlet temperature, CoolProp `INCOMP::MPG[0.26]` | sourced | Melinder data via CoolProp; within a few % of DOWFROST LC 25 on μ and k; Dow cp is ~4% higher | Viscosity falls ~15% over a 7 K rise, so ΔP is slightly overstated |
| P4 | 1D Nusselt number | Shah & London Nu_H1, four heated walls, with fin efficiency | sourced | Shah & London 1978 polynomial; cover is actually adiabatic | Part of the ≤ 2.9% parity gap |
| P5 | Flow distribution | Equal flow in every channel | assumed | No manifold model | Maldistribution raises peak temperature |
| P6 | Losses outside the channels | None (inlet, outlet, manifolds, fittings) | assumed | ΔP budget is defined as core-only | System ΔP is higher than reported |
| P7 | R_th definition | (average interface temperature − inlet) / Q; excludes TIM and junction-to-case | chosen | Isolates what the plate controls | — |
| P8 | Numerics | P2 triangles, 8 cells across the half channel, direct solve | checked per study | Mesh change < 0.01% at 2× refinement | — |
