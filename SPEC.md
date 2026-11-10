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
