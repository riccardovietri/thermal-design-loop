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
| R4 | Rated flow | 1.5 LPM (1.5 LPM per kW); not a design variable; 1.0 and 2.0 LPM reported as sensitivity on the recommended design | chosen | Riccardo, SPEC Amendment 1; industry guidance of 1.2–1.5 LPM per kW (Electronics Cooling 2025; CoolIT) | Sets the coolant-heating floor on R_th (5.0 mK/W at 1.5 LPM) |
| R5 | ΔP budget, channel core | 6 kPa at rated flow (was 30 kPa) | chosen | Riccardo's allocation, SPEC Amendment 1; an allocation, not a sourced specification. Consistent with published context: 124 kPa rack budget for 80 kW with 34 kPa in hoses and fittings; a 4 kW plate at 55 kPa for the full loop with fittings; a measured plate at 4 kPa | Sets the right-hand end of the Pareto front and which designs are eligible for recommendation |
| R6 | Plate material | Copper, k = 390 W/m·K | assumed | Typical for C110 | ~1% on R_th |
| R7 | Base thickness | 2 mm | assumed | Literature spans 0.1–4.5 mm, no typical value | Fixed 1.4 mK/W offset; ranking unchanged |
| R8 | Thermal-resistance target | R_th ≤ 15 mK/W | chosen | Riccardo, SPEC Amendment 1: 55 °C average interface at 1 kW and 40 °C inlet | Decides which rule fires; the recommended design is the lowest-ΔP one meeting it |

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
| P1 | Flow regime | Laminar, steady | checked per run | Re 110–180 on the shortlist at rated flow (80–160 in the flow-sensitivity runs); runs with Re > 2000 are rejected | — |
| P2 | Fully developed flow and heat transfer | Both 1D and CFD | assumed | x* = L/(Dh·Re·Pr) is 0.02–0.1 across the shortlist, so the thermal entrance region is not negligible for the widest channels | Real R_th lower and ΔP higher than reported; Level 2 quantifies it |
| P3 | Fluid properties | Constant, at inlet temperature, CoolProp `INCOMP::MPG[0.26]` | sourced | Melinder data via CoolProp; within a few % of DOWFROST LC 25 on μ and k; Dow cp is ~4% higher | Viscosity falls ~15% over a 7 K rise, so ΔP is slightly overstated |
| P4 | 1D Nusselt number | Shah & London Nu_H1, four heated walls, with fin efficiency | sourced | Shah & London 1978 polynomial; cover is actually adiabatic | Part of the ≤ 2.9% parity gap |
| P5 | Flow distribution | Equal flow in every channel | assumed | No manifold model | Maldistribution raises peak temperature |
| P6 | Losses outside the channels | None (inlet, outlet, manifolds, fittings) | assumed | ΔP budget is defined as core-only | System ΔP is higher than reported |
| P7 | R_th definition | (average interface temperature − inlet) / Q; excludes TIM and junction-to-case | chosen | Isolates what the plate controls | — |
| P8 | Numerics | P2 triangles, 8 cells across the half channel, direct solve | checked per study | Mesh change < 0.01% at 2× refinement | — |

## Level 2 (2D planar Navier-Stokes, slot-jet arrays)
| ID | Assumption | Value | Status | Basis | Effect if wrong |
|---|---|---|---|---|---|
| L1 | Flow is two-dimensional, steady and laminar | Slot length (60 mm) much greater than slot width; Re on slot width 5–20 in the rated-flow sweep | assumed | Re is far below reported transition for confined slot jets; no stability analysis done | End effects and any unsteadiness are missed |
| L2 | Slot-width Reynolds number depends only on slot count | Re = ρQ / (μ N L) | derived | Continuity | — |
| L3 | Heated surface | Flat, uniform flux applied directly to the fluid boundary | assumed | Simplest case; no fins under the jets | A conducting base would flatten the stagnation peak and lower the peak temperature; average changes little |
| L4 | Base conduction | Added as t/(kA) in series | assumed | Same 2 mm copper base as Level 1 | See L3 |
| L5 | Jet inlet profile | Parabolic (developed slot flow), temperature uniform at inlet value | assumed | Nozzle plate thicker than a few slot widths | Flat profile would lower stagnation h slightly |
| L6 | Pressure drop | Inlet-slot face to return-slot face; nozzle-plate and manifold losses excluded | assumed | Same core-only convention as Level 1 | System ΔP is higher |
| L7 | Return path | One return slot midway between jets, same width as the jet | chosen | Makes the unit cell periodic | Other return layouts change the wall-jet length |
| L8 | Minimum jet gap and slot width | 0.1 mm | chosen | Riccardo: 50 µm laser-cut features exist but only at low volume. A gap is set by a spacer or etch depth, not by fin skiving; tolerance and fouling risk not assessed | Gaps of 0.05 mm (excluded) gave about 16 mK/W in the earlier sweep; at ≥ 0.1 mm no slot-jet case reaches the 15 mK/W target |
| L9 | Numerics | Taylor-Hood P2–P1, P2 temperature, no stabilisation, 20 cells across min(slot, gap) | checked | Validated on thermally developing channel flow: ≤ 0.64% on local Nu, 0.03% energy balance; slot-jet mesh study (0.2 mm slot, 0.1 mm gap, 2 mm pitch) changes R_th by 0.002% and ΔP by 0.25% from 20 to 40 cells | — |
