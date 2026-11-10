# DECISIONS

Format: date — decision — why — alternatives rejected.

- 2026-11-10 — **2D formulation: conjugate cross-section unit cell, fully developed flow** — keeps all four design variables and tests exactly what the 1D model approximates (f·Re(α), Nu(α), fin efficiency, base conduction); exact analytical validation exists — rejected: streamwise parallel-plate plane (fin height drops out as a variable).
- 2026-11-10 — **Toolchain: gmsh (geometry + mesh) + scikit-fem (solve), pip-native via uv, Python 3.12** — installs natively on the M3 Pro with no compilers or Docker, fully scriptable, small enough to read, validated to 0.007% on f·Re and Nu_H1 (four aspect ratios, 1.6 s total) — rejected: OpenFOAM (needs Docker or a case-sensitive volume on macOS; overkill for a 2D linear problem), FEniCSx (conda-only, heavier), SU2/Elmer (not installed, compile or conda needed).
- 2026-11-10 — **3D path** — keep `make_geometry`/`mesh`/`run_cfd` as separate tools with file-based hand-off so a 3D solver (OpenFOAM or FEniCSx) can replace `run_cfd` later; gmsh already does 3D — not built in V1.
- 2026-11-10 — **1D screen lives in this repo, not thermal-mcp-server** — the server uses circular-tube constants (Nu = 4.36, f = 64/Re) with no fin efficiency; straight-fin model is added here first and upstreamed after the demo.
- 2026-11-10 — **1D Nu: Shah & London four-wall Nu_H1 polynomial + fin efficiency** — standard engineering practice and the constants are certain; the three-heated-wall table would be more faithful to the adiabatic cover and is a candidate upgrade once checked against the book — rejected for now: three-wall table from memory.
- 2026-11-10 — **1D stays fully developed (no entrance or minor-loss terms)** — parity plot should compare like with like; entrance effects are a stated limitation of both.
- 2026-11-10 — **Requirements: 1 kW uniform over 60×60 mm, PG25 at 40 °C inlet, ΔP ≤ 30 kPa core-only, copper k = 390, 2 mm base** — Riccardo's call.
- 2026-11-10 — **Manufacturing limits: ≥ 0.15 mm features, aspect ratios ≤ 10, fin height ≤ 4 mm (skived copper)** — Riccardo accepted defaults.
- 2026-11-10 — **Objectives: R_th vs ΔP (hero plot); pumping power reported only** — Riccardo's call.
- 2026-11-10 — **R_th = (avg base T − T_inlet)/Q, incl. base + caloric, excl. TIM and R_jc** — isolates what the plate design controls.
- 2026-11-10 — **License MIT** — matches thermal-mcp-server.
- 2026-11-10 — **Model use: Fable as main model, no separate advisor** — Riccardo's call; revisit at each usage check (delegate coding to Opus subagents if Fable usage passes ~50% before the Monday reset).
- 2026-11-10 — **PG25 properties from CoolProp `INCOMP::MPG[0.26]` (Melinder data) at 40 °C** — reproducible in code and citable; cross-checked against DOWFROST LC 25 typical data (μ = 1.15 mPa·s and k = 0.485 W/m·K at 50 °C vs CoolProp 1.18 and 0.487 at 47 °C) — rejected: hand-typed datasheet values; note inhibited commercial fluids run about 1–2% denser and Dow lists a higher cp (4.13 vs 3.98 kJ/kg·K), not resolved.
- 2026-11-10 — **SPEC.md frozen** — Riccardo approved.
