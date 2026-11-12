# SPEC — Cold-plate design loop (V1)

Status: **FROZEN 2026-11-10 (approved by Riccardo).** Changes need his sign-off.

## Purpose
Automate an engineer's design loop for a direct-to-chip cold plate:
need → design space → 1D screen → CFD on the shortlist → confidence checks → Pareto front → recommendation.
Every step is a plain function/CLI that runs without an LLM. One command reproduces the study.

## Requirements case (`examples/gpu_coldplate.yaml`)
| Item | Value |
|---|---|
| Heat load | 1,000 W, uniform flux over the finned footprint (27.8 W/cm²) |
| Footprint | 60 × 60 mm, flow length 60 mm |
| Coolant | PG25 (25 vol% ≈ 26 mass% propylene glycol), 40 °C inlet, constant properties at 40 °C from CoolProp `INCOMP::MPG[0.26]`: ρ = 1010.6 kg/m³, cp = 3958 J/kg·K, k = 0.481 W/m·K, μ = 1.394 mPa·s (Pr = 11.5) |
| Flow rate | 1–2 LPM (design variable) |
| ΔP budget | ≤ 30 kPa, channel core only (no manifolds/QDs), at 40 °C inlet |
| Plate | Copper, k = 390 W/m·K, 2 mm base (fixed), cover adiabatic, no tip clearance |

## Design variables
| Variable | Range |
|---|---|
| Channel width `w_ch` | 0.15–1.0 mm |
| Fin thickness `t_fin` | 0.15–1.0 mm |
| Fin height `h_fin` | 1.0–4.0 mm |
| Flow rate | 1.0–2.0 LPM |

Channel count = floor(60 mm / (w_ch + t_fin)).

## Constraints
- Core ΔP ≤ 30 kPa.
- Manufacturing (skived copper): w_ch ≥ 0.15 mm, t_fin ≥ 0.15 mm, h_fin / t_fin ≤ 10, h_fin / w_ch ≤ 10.

## Objectives
- Minimise R_th and core ΔP. **Hero plot: R_th vs ΔP.** Hydraulic pumping power is reported, not optimised.
- R_th = (area-average base temperature at the chip interface − T_inlet) / Q.
  Includes base conduction, fin/convection and caloric rise; excludes TIM and R_jc.
  The max-temperature version (outlet end) is also reported.

## Physics
### "CFD" (2D conjugate cross-section, fully developed)
Unit cell: half fin + half channel + base, symmetry on both sides.
- Momentum: μ∇²w = dp/dx in the fluid, no-slip on fin, base and cover.
- Energy: fluid ρ·cp·w·∂T/∂x = k_f∇²T; solid k_s∇²T = 0; ∂T/∂x constant from the overall energy balance.
- BCs: uniform flux on the base underside, adiabatic cover, symmetry sides, continuity at fluid–solid interfaces.
- Laminar only (expected Re ≈ 100–700); a design with Re > 2000 is rejected, not solved.
- Not captured: hydrodynamic/thermal entrance effects, inlet/outlet losses, flow maldistribution,
  property variation with temperature, axial variation of heat flux.

### 1D screen (in this repo; upstream to thermal-mcp-server after the demo)
- f·Re(α): Shah & London rectangular-duct polynomial.
- Nu(α): Shah & London Nu_H1 polynomial (four heated walls), applied with adiabatic-tip fin efficiency
  on the fin sides and η = 1 on the channel floor; cover adiabatic.
- R_th = R_base + R_fin/conv + caloric (Q / 2ṁcp). Fully developed, to match the CFD.
- thermal-mcp-server supplies the resistance-network structure; its circular-tube constants are not used.

## Pipeline
1. `screen_1d`: full-factorial/LHS sweep of ≥ 300 designs, filter by constraints.
2. `pareto`: 1D Pareto front; pick ~10 designs spread along it.
3. `make_geometry` + `mesh`: parametric unit cell in gmsh; mesh-quality metrics recorded.
4. `run_cfd`: scikit-fem solve; writes `runs/<design_id>/` (inputs, logs, results).
5. `check_confidence`: mesh independence (3 levels, one design), linear-solver residual, energy-balance error.
6. `postprocess` + `report`: CFD Pareto plot, 1D-vs-CFD parity plot, recommendation, `report.md`.

## Deliverables
- `uv run python run_study.py examples/gpu_coldplate.yaml` → Pareto plot, parity plot, `report.md`.
- `tools/` with typed, JSON-serialisable functions and a CLI each.
- Tests for the 1D screen and post-processing; confidence checks in the report.
- Public GitHub repo (MIT) with README: problem, 3-step diagram, hero figure, how to run.

## Out of scope (V1)
3D, developing flow in the main study, FEA, web UI, cloud runs, renaming thermal-mcp-server,
autonomous agent. Stretch: MCP wrapper over the same tools; one streamwise case to bound entrance error.

## Amendment 1 (2026-11-10, approved by Riccardo)
Supersedes the flow-rate, ΔP-budget and recommendation items above where they conflict; everything else stands.
- **Rated flow** is 1.5 LPM (1.5 LPM per kW). The Pareto study is run at rated flow only; flow is no longer a design variable.
  Flow sensitivity (1.0 and 2.0 LPM) is reported for the recommended design.
- **Pressure budget:** channel core ΔP ≤ 6 kPa at rated flow (replaces 30 kPa).
- **Thermal target:** R_th ≤ 15 mK/W.
- **Recommendation rule:** among CFD-verified designs that meet both the ΔP budget and the R_th target, recommend the one with
  the lowest ΔP. If none meets the target, recommend the lowest R_th within budget and say so. The knee (existing definition:
  the point on the CFD-verified front farthest from the chord between its end points in log-log space) is still computed and
  reported as secondary information.
- **Plots** show the region of interest only.

## Amendment 2 (2026-11-11, requested by Riccardo)
Adds to scope; everything above stands.
- **Scope now includes** a plan-view Level 2 evaluation of sinusoidal, zigzag, offset-strip and pin-fin (circular and diamond)
  planforms, with straight fins computed the same way as the reference, against the Amendment 1 requirements: 1 kW uniform over
  60 × 60 mm, PG25 at 40 °C, rated flow 1.5 LPM, channel-core ΔP ≤ 6 kPa, R_th ≤ 15 mK/W, copper k = 390 W/m·K, 2 mm base.
- **Manufacturing limits for these planforms:** minimum solid feature 0.15 mm, minimum fluid gap 0.15 mm, fin height ≤ 4 mm,
  fin height / minimum fin thickness ≤ 10.
- **Method:** steady laminar 2D Navier-Stokes + energy in the plane of the fin layer (tall-fin limit), developed-periodic
  pressure gradient and heat-transfer coefficient, put into the Level 1 resistance network (base + convection with fin efficiency
  + coolant heating). Entrance effects excluded, as in Level 1.
- **Level 1 pipeline unchanged** (`run_study.py`, `tools/`).

## Amendment 3 (2026-11-11, approved by Riccardo)
Supersedes the ΔP-budget, thermal-target and recommendation items of Amendment 1 (and their use in Amendment 2); everything else stands.
1. **Rated flow** stays 1.5 LPM.
2. **Two-tier thermal requirement:** threshold R_th ≤ 20 mK/W (must meet) and objective R_th ≤ 15 mK/W (want). This replaces the
   single 15 mK/W target.
3. **Pressure budget:** channel-core ΔP ≤ 5 kPa at rated flow (replaces 6 kPa), one third of a 15 kPa whole-plate allocation.
4. **Recommendation by marginal return** (replaces "lowest ΔP that meets the target" and the knee):
   - Candidates are the CFD-verified designs with ΔP ≤ budget and R_th ≤ threshold. Take their non-dominated front, sorted by ΔP
     ascending.
   - Start at the lowest-ΔP front point. From the current point i, look at every front point j further along and compute the
     return e(i→j) = −ln(R_j / R_i) / ln(ΔP_j / ΔP_i). If any j has e ≥ the cut-off, move to the nearest such j and repeat.
     Stop when no further point gives e ≥ the cut-off relative to the current point, and recommend the point where the walk stopped.
   - The cut-off is 0.25 (a doubling of ΔP must buy at least about 16% lower R_th); it is set in the requirements file as
     `marginal_return_cutoff`.
   - The report states whether the recommended design meets the objective.
   - If no design meets the threshold within budget, the lowest R_th within budget is recommended and the report says clearly that
     the threshold is not met. If no design is within budget, the study stops with an error.
   - The rule is a pure function in `tools/pareto.py` that returns the index, the rule name and the steps taken with their returns.
     The old `recommend` and `knee_index` functions stay in the code but are no longer used in the reports.
5. **Extra columns for every CFD-verified design** in the reports and results JSON: capacity_kw = threshold × heat load / R_th / 1000
   (the heat load the design could carry at the threshold base-to-inlet rise, 20 K at 1 kW); capacity_kw_per_kpa = capacity_kw / ΔP
   in kPa; kw_per_pumping_w = heat load / (ΔP × volumetric flow) / 1000, using hydraulic power. The definitions are stated in the
   report, with a note that both ratios always favour the lowest ΔP inside the window, which is why they are reported and not used
   to choose.
6. **Level 1 aspect-ratio limit H / gap ≤ 10 applies to the plan-view results** as well. Designs that break it are excluded from the
   fronts and the recommendation and listed as excluded with the reason. Solves are not repeated for this; a topology left with
   fewer than 5 eligible designs gets added designs to restore about 8.
