# Build log

Raw material for the write-up (README, portfolio page, LinkedIn post). One entry per working session:
what was tried, what the numbers were, what surprised us. Keep figures and commands reproducible.

## 2026-11-10 (Tue) — spec, toolchain, first solve
- Audit of the existing 1D model found circular-tube constants (Nu = 4.36, f = 64/Re) on rectangular
  channels and no fin efficiency, so the straight-fin screen is written fresh in this repo.
- Hand estimate before any code: channel-core ΔP is 0.6–5 kPa for 0.2–0.5 mm channels, Re 100–400.
  The 30 kPa budget only binds near 0.15 mm channels; no turbulence model is needed.
- Formulation: 2D conjugate cross-section, fully developed. Chosen because it keeps all four design
  variables and has exact reference solutions.
- Toolchain: gmsh + scikit-fem, pip-only. Validation against Shah & London at four aspect ratios:
  worst error 0.007% on f·Re and Nu_H1, 1.6 s total.
- First conjugate solves (PG25, 40 °C, 1 kW, 60×60 mm), 0.1–0.3 s each. The 1D column is a hand
  estimate (Shah & London polynomials + fin efficiency), not yet the packaged screen:

| Design (w × t × h mm, LPM) | Re | ΔP CFD / 1D (kPa) | R_th CFD / 1D (mK/W) |
|---|---|---|---|
| 0.5 × 0.5 × 2.0, 1.5 | 242 | 1.99 / 1.99 | 25.3 / 26.4 |
| 0.15 × 0.15 × 1.5, 1.5 | 110 | 26.45 / 26.46 | 9.31 / 9.24 |
| 0.3 × 0.3 × 3.0, 2.0 | 146 | 4.41 / 4.41 | 10.94 / 10.80 |
| 1.0 × 1.0 × 4.0, 1.0 | 161 | 0.17 / 0.17 | 46.7 / 48.9 |

- Observation: ΔP agrees to four figures because both models are fully developed with the same f·Re.
  The informative parity axis is R_th (−1% to +5% so far); ΔP parity only becomes interesting with
  entrance effects.
