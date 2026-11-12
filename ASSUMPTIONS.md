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
| R4 | Rated flow | 1.5 LPM (1.5 LPM per kW); not a design variable; 1.0 and 2.0 LPM reported as sensitivity on the recommended design | sourced | OCP OAI System Liquid Cooling Guidelines (March 2023): typical design target 10 K coolant rise, "1.5 LPM/kW", range 1.25–2.0 LPM/kW for PG25 (docs/literature.md, second pass); SPEC Amendments 1 and 3 | Sets the coolant-heating floor on R_th (5.0 mK/W at 1.5 LPM) |
| R5 | ΔP budget, channel core | 5 kPa at rated flow (was 30, then 6 kPa) | chosen | One third of a 15 kPa whole-plate allocation (SPEC Amendment 3); an allocation, not a sourced specification. Context: a measured skived copper plate (0.2 mm fins and channels, 4 mm tall, 43 mm long; Ortega et al., ASME InterPACK 2022) showed 43 kPa at 2.6 LPM for the whole plate including manifolds; an estimate of that plate's channel core alone at 1.5 LPM is about 4.6 kPa assuming a 50 mm plate width (the width is not given in the paper, so this is indicative); a vendor listing gives under 20 kPa at 1.8 L/min for a GB200-class plate (unverified); OCP cold-plate documents define what to report, not limits | Sets the right-hand end of the Pareto front and which designs are eligible for recommendation |
| R6 | Plate material | Copper, k = 390 W/m·K | assumed | Typical for C110 | ~1% on R_th |
| R7 | Base thickness | 2 mm | assumed | Literature spans 0.1–4.5 mm, no typical value | Fixed 1.4 mK/W offset; ranking unchanged |
| R8 | Thermal-resistance requirement, two tiers | Threshold R_th ≤ 20 mK/W (must meet); objective R_th ≤ 15 mK/W (want) | chosen | SPEC Amendment 3. Threshold: the OCP OAI guidelines report under 0.02 °C/W case-to-inlet including the interface material for a 1 kW module on a 30 × 30 mm heated test vehicle, PG25, 2–5 LPM, "ideal conditions"; a plate alone must beat that. Objective: Ortega et al. measured 12 mK/W base-to-inlet at 2.6 LPM on a skived plate; 15 leaves roughly 5 mK/W for the interface material against the OCP figure | Moves the acceptable window |
| R9 | Marginal-return cut-off for the recommendation rule | 0.25 (a doubling of ΔP must buy at least about 16% lower R_th) | chosen | Judgment; Riccardo approved (SPEC Amendment 3); `marginal_return_cutoff` in the requirements file | A higher cut-off recommends a lower-ΔP design |

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
| L8 | Minimum jet gap and slot width | 0.1 mm | chosen | Riccardo: 50 µm laser-cut features exist but only at low volume. A gap is set by a spacer or etch depth, not by fin skiving; tolerance and fouling risk not assessed | Gaps of 0.05 mm (excluded) gave about 16 mK/W in the earlier sweep; at ≥ 0.1 mm no slot-jet case reaches the 20 mK/W threshold (best 23.8 mK/W) |
| L9 | Numerics | Taylor-Hood P2–P1, P2 temperature, no stabilisation, 20 cells across min(slot, gap) | checked | Validated on thermally developing channel flow: ≤ 0.64% on local Nu, 0.03% energy balance; slot-jet mesh study (0.2 mm slot, 0.1 mm gap, 2 mm pitch) changes R_th by 0.002% and ΔP by 0.25% from 20 to 40 cells | — |

## Level 2 (plan view: fin planforms in the tall-fin limit)
| ID | Assumption | Value | Status | Basis | Effect if wrong |
|---|---|---|---|---|---|
| V1 | Tall-fin limit | 2D solve in the plane of the fin layer; floor and cover drag and heat transfer are not in the 2D solve; the floor gets the side-wall h | assumed, measured | Gate G3: for the Level 1 recommended straight design (0.422 / 0.444 / 4.0 mm) the plan view gives ΔP 7.1% and R_th 5.7% below the Level 1 CFD, which includes floor and cover | Plan-view ΔP and R_th read low by roughly that amount for every topology; larger for short fins (H/gap small), where the floor and cover matter more |
| V2 | Lateral symmetry, steady laminar flow | Strip topologies (offset strips, pins) solved on half a lateral pitch with symmetry planes; channels solved as one channel; steady Newton solution | assumed | Re_Dh 71–216 (interstitial velocity, plan Dh) at rated flow; pins Re_d,max 50–136 (gap velocity, pin diameter or diamond side); a steady symmetric solution exists and Newton converged for all 63 designs, but no stability analysis is done | Staggered pins and strips may shed vortices at these Re (an isolated cylinder sheds above Re_d ≈ 47; arrays delay it); shedding would raise h and ΔP; 3D Dean vortices in wavy channels are absent in 2D, so wavy-fin benefit is understated |
| V3 | Wall heat flux | Uniform flux on all fin walls in the 2D solve; one mean h per design = q'' / (period-mean wall T − period-mean velocity-weighted bulk T) | assumed | Standard developed-periodic definition (Patankar, Liu & Sparrow 1977); reduces to Nu = 8.235 for straight fins | Real fins have non-uniform flux (more near the base); local h variation around pins is lumped into one value |
| V4 | Floor heat transfer | Floor (plan fluid area) gets the same h as the fin sides, η = 1 | assumed | Same structure as Level 1 | The floor is a different surface (boundary layer of the end wall); effect is part of V1 |
| V5 | Fin efficiency | Adiabatic tip, η = tanh(mH)/(mH), m = sqrt(h P_c / (k_s A_c)) with the plan-section perimeter and area of one fin | assumed | Standard for pins and strips; for a continuous wall it is the Level 1 formula | Ignores spreading in the base under sparse pins and 2D conduction at strip ends; pins at H = 4 mm have η 0.4–0.6, so R_th is sensitive to it |
| V6 | Developed-periodic values over the whole length | No entrance effect; checked with a developing strip from a uniform inlet for one design per topology | assumed, measured | Pins and offset strips reach the periodic h within 2% in 3–5 periods (4–10 mm); wavy channels need 10 periods (20 mm); straight channels from the 20–25 mm window | Real h is higher near the inlet, so R_th is conservative; the length-mean h of a straight channel over 60 mm is about 8% above the developed value (`validation/level2_planview_straight.py`); comparable for wavy channels, smaller for pins and strips |
| V7 | Manufacturing limits | Solid feature ≥ 0.15 mm, fluid gap ≥ 0.15 mm (true minimum distance between solids), H ≤ 4 mm, H / minimum fin thickness ≤ 10; Level 1's H / gap ≤ 10 is now applied too (SPEC Amendment 3), with gap = the true minimum fluid gap: designs that break it are solved but excluded from fronts and the recommendation and listed with the reason | chosen | Riccardo, SPEC Amendments 2 and 3 | Excludes 36 of 78 converged designs, including every design with gaps under 0.4 mm at H = 4 mm |
| V8 | Planform construction | Wavy walls offset ± w/2 along the normal (constant normal gap); lateral pitch set so the thinnest fin section is t; zigzag corners rounded to a 0.05 mm inner wall radius; offset strips with g = 0 and rows not overlapping laterally; pin cells cut through row-A pin centres | chosen | Simplest constructions that meet the brief; sharp-corner zigzags would put a singularity in the offset wall | Fins are thicker than t at wavy crests (more solid, fewer channels per width) |
| V9 | Lateral layout | Continuum: flow per lateral unit = Q × (unit width) / 60 mm; no integer channel count, no edge channels | assumed | Same footprint and total flow for every topology | For straight fins the Level 1 integer count (69 channels in 59.75 mm) differs by 0.4% |
| V10 | Numerics | Streamwise-periodic cell (one period), Taylor-Hood P2–P1, P2 temperature, no stabilisation; wall cell ref_gap/20, core ref_gap/10, graded over ref_gap/2 | checked | Straight channel: f·Re 24.0000, Nu 8.2353 (exact 8.2353). Mesh 10 → 20 cells on one design per topology: dp/dx ≤ 0.43%, h ≤ 0.55%. Energy balance ≤ 0.75% (diamond pins, sharp corners; falls to 0.24% at 2× refinement) and ≤ 0.24% elsewhere; thermal undershoot ≤ 0.3% of the wall-to-bulk difference. Diamond pins: the two windows of a two-period cell differ by 1.7% in h (corner sensitivity; the cell means agree to 0.09%) | — |
| V11 | Manufacturability of non-straight planforms | Wavy, zigzag, offset-strip and pin planforms are assumed makeable at the same 0.15 mm limits as skived straight fins | open | They cannot be made by skiving. OCP lists pin and offset fins with 0.1–1.0 mm internal features as options, but no volume process or feature limit at 0.15 mm in copper was found (docs/literature.md, second pass) | A non-straight recommendation carries unassessed manufacturing risk |
