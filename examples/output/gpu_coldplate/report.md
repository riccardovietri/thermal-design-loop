# Cold-plate design study: gpu_coldplate

**Recommended design:** 0.44 mm channels, 0.47 mm fins, 4.0 mm tall (65 channels) at 1.5 LPM. CFD gives R_th = 15.7 mK/W at ΔP = 1.2 kPa, so the chip interface averages 55.7 °C (peak 60.7 °C) for 1000 W and 40 °C inlet.

## Why this design
- It is the knee of the CFD-verified front: the point farthest from the straight line between the two ends of the front in log-log space.
- The lowest-resistance design within budget (0.15 mm channels, 0.15 mm fins, 1.5 mm tall, 1.5 LPM) reaches 9.3 mK/W at 26.4 kPa: 41% lower resistance for 22.2× the pressure drop.
- At the recommended point the resistance splits into base conduction 1.4, fin and convection 9.2, and coolant heating 5.0 mK/W (1D breakdown).

![Pareto](pareto.png)

## 1D screen vs CFD
The 1D screen evaluated 1008 designs (530 feasible) in 3 ms. 10 designs spread along its Pareto front were solved with the 2D conjugate model in 1.6 s.

- R_th: 1D differs from CFD by 0.9% on average, 2.9% at worst (mean signed +0.0%).
- ΔP: identical to within 0.03%. This is expected, not a result: both models assume fully developed flow and share the same f·Re.

![Parity](parity.png)

| Channel (mm) | Fin (mm) | Height (mm) | Flow (LPM) | Re | ΔP CFD (kPa) | R_th 1D (mK/W) | R_th CFD (mK/W) | 1D − CFD | Fin eff. |
|---|---|---|---|---|---|---|---|---|---|
| 1.00 | 0.47 | 4.0 | 1.0 | 121 | 0.12 | 40.19 | 39.05 | +2.9% | 0.92 |
| 0.76 | 0.47 | 4.0 | 1.0 | 106 | 0.22 | 29.12 | 28.72 | +1.4% | 0.89 |
| 0.58 | 0.47 | 4.0 | 1.0 | 93 | 0.41 | 22.19 | 22.14 | +0.2% | 0.86 |
| 0.44 | 0.47 | 4.0 | 1.0 | 84 | 0.79 | 18.11 | 18.21 | -0.6% | 0.82 |
| 0.44 **(recommended)** | 0.47 | 4.0 | 1.5 | 126 | 1.19 | 15.61 | 15.71 | -0.6% | 0.82 |
| 0.34 | 0.32 | 3.0 | 1.5 | 119 | 2.56 | 13.17 | 13.27 | -0.8% | 0.81 |
| 0.34 | 0.32 | 3.0 | 2.0 | 159 | 3.41 | 11.92 | 12.02 | -0.8% | 0.81 |
| 0.26 | 0.32 | 2.5 | 2.0 | 170 | 8.10 | 10.54 | 10.60 | -0.5% | 0.82 |
| 0.20 | 0.15 | 1.5 | 1.5 | 123 | 13.85 | 10.51 | 10.55 | -0.4% | 0.82 |
| 0.15 | 0.15 | 1.5 | 1.5 | 110 | 26.45 | 9.24 | 9.31 | -0.8% | 0.78 |

## Confidence checks

| Check | Value | Limit | Result |
|---|---|---|---|
| Mesh independence, R_th (production vs 2× finer) | 0.001% | < 1.0% | pass |
| Mesh independence, ΔP (production vs 2× finer) | 0.000% | < 1.0% | pass |
| Linear solver residual, worst run | 9.1e-12 | < 1e-08 | pass |
| Energy balance: flux recovered at chip interface vs applied, worst run | 0.000% | < 1.0% | pass |
| Laminar assumption: highest Reynolds number | 170 | < 2000 | pass |
| Mesh quality: lowest element quality (SICN) | 0.76 | > 0.3 | pass |

Mesh study on `w0444_t0468_h4000_q150`:

| Cells across half channel | Unknowns | R_th (mK/W) | ΔP (kPa) |
|---|---|---|---|
| 4 | 4879 | 15.7090 | 1.1894 |
| 8 | 17951 | 15.7096 | 1.1894 |
| 16 | 68255 | 15.7098 | 1.1894 |

![Fields](fields_recommended.png)

## Inputs
- Heat load 1000 W, uniform over 60 × 60 mm; 2.0 mm base, k = 390 W/m·K; ΔP budget 30 kPa (channel core).
- PG25 at 40 °C (CoolProp): ρ = 1010.6 kg/m³, cp = 3958 J/kg·K, k = 0.481 W/m·K, μ = 1.394 mPa·s.
- R_th = (average chip-interface temperature − inlet temperature) / heat load. TIM and junction-to-case are excluded.

## Design basis
- Feature limits (≥ 0.15 mm, aspect ratio ≤ 10) sit inside published skived-copper practice: 0.1–0.8 mm fins with 0.1 mm minimum gap (Boyd fabrication guide); 100–200 µm fins and gaps, 1–2 mm tall, in direct-to-chip cold plates (ASME mini-channel study).
- Base thickness is a fixed assumption. Published designs span 0.1 mm (manifold microcoolers) to 4.5 mm (a commercial 55 mm plate); no typical value was found.

## What this model does not capture
- Entrance effects and inlet/outlet losses (both models are fully developed).
- Flow maldistribution between channels, manifolds and fittings.
- Temperature-dependent properties, non-uniform heat flux and heat spreading from a die smaller than the plate.
- Anything three-dimensional. Each CFD run is a 2D cross-section of one fin and channel.

Every CFD number traces to `runs/gpu_coldplate/<design_id>/` (inputs, mesh, results, fields).
