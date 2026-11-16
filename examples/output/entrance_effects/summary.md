# Entrance effects in the straight-fin study (issue #12)

**Answer.** For the recommended design (0.42 mm channels, 0.44 mm fins, 4.0 mm tall, 69 channels, Re 119, x* = L/(Dh·Re·Pr) = 0.058), developing flow lowers R_th from 14.93 to 14.44 mK/W (-3.3%) and raises the channel-core ΔP from 1.294 to 1.320 kPa (+2.0%). Inlet and outlet losses, which are outside the channel-core definition, add 19 Pa more (1.339 kPa in all). **The 15 mK/W objective is still met**, with a margin of 3.7% instead of 0.5%. The resolved check (tall-fin limit) does not confirm the 1D increments to within the 10% gate (see below), but every estimate moves R_th down and ΔP up by a few percent: with the resolved ratios in place of the 1D ones the recommended design reads 14.55 mK/W at 1.312 kPa (core) and 1.317 kPa plenum to plenum, so the conclusion does not depend on which is right.

**Front and recommendation.** The order of the 8 CFD-verified designs by R_th is unchanged; the wide channels gain most (x* is smallest there), so the front flattens. 0.60 / 0.44 / 4.0 mm now meets the 20 mK/W threshold (18.58 mK/W), so the marginal-return walk starts one point earlier: 0.60 / 0.44 / 4.0 → 0.50 / 0.44 / 4.0 → 0.42 / 0.44 / 4.0 mm, and stops at the same design (fully developed: 0.50 / 0.44 / 4.0 → 0.42 / 0.44 / 4.0 mm). With inlet/outlet losses added to ΔP the walk ends at the same design. The recommendation would not change.

Values in the tables are the Level 1 CFD result plus the 1D increment (CFD-anchored): the 1D screen with the new option, minus the same screen fully developed. Geometry is channel width / fin thickness / fin height in mm.

## Sources and formulas
All read in full text; see `docs/literature.md` and `ASSUMPTIONS.md` rows P2 and P9–P12.
- **Hydrodynamic entrance:** Muzychka & Yovanovich (2002), Part I. Apparent friction (mean over the length, including the momentum change of the developing profile, so no separate K(∞) term): f_app·Re_√A = [(f·Re_√A)² + (3.44/√L+)²]^½, L+ = L/(√A·Re_√A), f·Re_√A = 12 / (√ε (1+ε) [1 − 192ε/π⁵ tanh(π/2ε)]). Stated accuracy ±10% on f_app·Re for most duct shapes.
- **Thermal and simultaneous entrance, uniform wall flux:** Muzychka & Yovanovich (2002), Part II, Eqs. 50–52 with their Table 2 (C1 = 3.66, C3 = 0.501, γ = 1/10, f(Pr) = 0.886/[1 + (1.909 Pr^(1/6))^(9/2)]^(2/9), m = 2.27 + 1.65 Pr^(1/3)), local form (C2 = C4 = 1), z* = z/(√A·Re_√A·Pr). Valid for 0.1 < Pr < ∞; stated accuracy ±12% (Graetz) and ±15% (combined entry). R_th uses the mean of (wall − bulk) temperature over the length, so the effective h is the harmonic length-mean of the local h. Both corrections are applied as ratios to the model's own developed limit, times the Shah & London f·Re and Nu_H1 already in the screen, so they vanish for a long channel.
- **Inlet and outlet losses:** laminar abrupt-contraction and abrupt-expansion coefficients of Kays (1950) for flat and square ducts, as polynomial fits in σ interpolated linearly in aspect ratio by dos Reis & Altemani (COBEM 1999, Eq. 5, Table 1); ΔP_io = (Kc + Ke)·ρU²/2 with U the channel velocity and σ = w/(w + t) (plenum as tall as the fins, no other area change). Kays's paper itself was not obtained. The flat-duct expansion fit equals the momentum balance with a parabolic exit profile (1 − 2.4σ + σ², tested).

## R_th and ΔP per design

| Design | x* | R_th fully developed | thermal entrance | hydrodynamic entrance | inlet/outlet | combined |
|---|---|---|---|---|---|---|
| 1.00 / 0.44 / 4.0 | 0.018 | 35.98 | 28.63 | 35.98 | 35.98 | 28.55 |
| 0.84 / 0.44 / 4.0 | 0.023 | 29.09 | 24.74 | 29.09 | 29.09 | 24.69 |
| 0.71 / 0.44 / 4.0 | 0.029 | 23.69 | 21.21 | 23.69 | 23.69 | 21.18 |
| 0.60 / 0.44 / 4.0 | 0.036 | 20.05 | 18.59 | 20.05 | 20.05 | 18.58 |
| 0.50 / 0.44 / 4.0 | 0.046 | 17.12 | 16.28 | 17.12 | 17.12 | 16.27 |
| 0.42 / 0.44 / 4.0 **(recommended)** | 0.058 | 14.93 | 14.45 | 14.93 | 14.93 | 14.44 |
| 0.36 / 0.44 / 3.5 | 0.065 | 13.78 | 13.41 | 13.78 | 13.78 | 13.40 |
| 0.30 / 0.26 / 2.5 | 0.081 | 12.58 | 12.30 | 12.58 | 12.58 | 12.30 |

R_th in mK/W. Thermal entrance: Graetz model (developed velocity, developing temperature). Combined: velocity and temperature developing together from a uniform inlet. The hydrodynamic and inlet/outlet effects do not change R_th in this model.

| Design | Re | ΔP fully developed | thermal entrance | hydrodynamic entrance | inlet/outlet | combined, core | combined with inlet/outlet |
|---|---|---|---|---|---|---|---|
| 1.00 / 0.44 / 4.0 | 177 | 182 | 182 | 196 | 187 | 196 | 202 |
| 0.84 / 0.44 / 4.0 | 163 | 264 | 264 | 280 | 271 | 280 | 287 |
| 0.71 / 0.44 / 4.0 | 148 | 382 | 382 | 399 | 391 | 399 | 408 |
| 0.60 / 0.44 / 4.0 | 138 | 573 | 573 | 593 | 585 | 593 | 604 |
| 0.50 / 0.44 / 4.0 | 128 | 856 | 856 | 879 | 871 | 879 | 893 |
| 0.42 / 0.44 / 4.0 **(recommended)** | 119 | 1294 | 1294 | 1320 | 1313 | 1320 | 1339 |
| 0.36 / 0.44 / 3.5 | 125 | 2277 | 2277 | 2317 | 2309 | 2317 | 2349 |
| 0.30 / 0.26 / 2.5 | 121 | 3794 | 3794 | 3850 | 3830 | 3850 | 3886 |

ΔP in Pa. Inlet/outlet column: fully developed core plus (Kc + Ke)·ρU²/2. Coefficients and factors:

| Design | σ | Kc | Ke | h factor, thermal | h factor, combined | f_app / f | L/(Dh·Re) |
|---|---|---|---|---|---|---|---|
| 1.00 / 0.44 / 4.0 | 0.693 | 0.709 | -0.252 | 1.351 | 1.356 | 1.0808 | 0.212 |
| 0.84 / 0.44 / 4.0 | 0.655 | 0.715 | -0.198 | 1.260 | 1.264 | 1.0608 | 0.265 |
| 0.71 / 0.44 / 4.0 | 0.615 | 0.723 | -0.141 | 1.189 | 1.192 | 1.0452 | 0.337 |
| 0.60 / 0.44 / 4.0 | 0.573 | 0.732 | -0.082 | 1.140 | 1.142 | 1.0346 | 0.418 |
| 0.50 / 0.44 / 4.0 | 0.531 | 0.742 | -0.019 | 1.103 | 1.104 | 1.0262 | 0.527 |
| 0.42 / 0.44 / 4.0 | 0.488 | 0.752 | 0.047 | 1.075 | 1.076 | 1.0200 | 0.661 |
| 0.36 / 0.44 / 3.5 | 0.445 | 0.766 | 0.112 | 1.066 | 1.067 | 1.0177 | 0.742 |
| 0.30 / 0.26 / 2.5 | 0.537 | 0.737 | -0.026 | 1.058 | 1.059 | 1.0148 | 0.928 |

![Straight-fin front with and without entrance effects](front_entrance.png)

*R_th against ΔP for the 8 CFD-verified straight-fin designs of the Level 1 study at 1.5 LPM PG25: fully developed (Level 1 CFD, blue) and with developing flow in the channel core (copper, arrows join each design), and with inlet and outlet losses added to ΔP (open diamonds). Corrections are the 1D increments added to the CFD values. Budget 5 kPa (channel core), threshold 20 mK/W, objective 15 mK/W. Circled: the Level 1 recommended design (0.42 / 0.44 / 4.0 mm). Region of interest as in the Level 1 report (from 0.3 kPa): 1.00 / 0.44 / 4.0, 0.84 / 0.44 / 4.0 mm (below 0.3 kPa) are off the plot and in the tables.*

## Resolved check (tall-fin limit)

Three solves per design with the plan-view solver (`level2/planar.py`, the solver behind `level2.planview.run_developing`), uniform flux on the fin walls, inlet temperature uniform:
1. **Combined entry:** one straight half channel over the full 60 mm, uniform velocity 1 mm upstream of the fin leading edge with slip walls there (a uniform profile imposed right at the wall corner makes the inlet pressure singular, and the pressure drop then grows with mesh refinement; the slip lead-in removes that).
2. **Thermal entrance only:** the same channel with a parabolic inlet profile at x = 0.
3. **Plenum to plenum (ΔP only):** a 2 mm inlet plenum, half a fin of thickness t with square ends over 60 mm, a 10 mm outlet plenum; plenums as tall as the fins (σ = w/(w + t)). The pressure difference changes by under 0.1% between 5 and 10 mm after the array.

Length-mean h = q''/(mean wall T − mean bulk T over 0–60 mm); ΔP from the inlet line to x = 60 mm. Developed reference: the streamwise-periodic cell. The 1D column evaluates the new option for the same channel: parallel plates of gap w, aspect ratio 0.01, as the source models parallel plates.

**Tall-fin limit.** The plan view has fins of infinite height: no floor, no cover, both channel walls heated. So it checks development between two fins, not the real duct with an aspect ratio of about 0.1 and an adiabatic cover. On the developed values it reads 5.7% low on R_th and 7.1% low on ΔP for the 0.42 / 0.44 / 4.0 mm design (ASSUMPTIONS V1). Corner regions of the real duct develop differently and are not checked here; the 1D option applied to the real duct uses that duct's aspect ratio, so the two are compared only in the parallel-plate limit.

| Design | Re (2w) | x* (2w) | h increase, resolved | h increase, 1D | ΔP increase, resolved | ΔP increase, 1D | ΔR_th resolved (mK/W) | ΔR_th 1D (mK/W) | energy balance |
|---|---|---|---|---|---|---|---|---|---|
| 1.00 / 0.44 / 4.0 | 218 | 0.012 | +25.2% | +13.7% | +5.19% | +7.26% | -4.553 | -2.725 | +0.037% |
| 0.50 / 0.44 / 4.0 | 143 | 0.037 | +7.8% | +4.7% | +1.73% | +2.44% | -0.572 | -0.358 | +0.016% |
| 0.42 / 0.44 / 4.0 | 131 | 0.047 | +5.9% | +3.7% | +1.34% | +1.89% | -0.345 | -0.220 | +0.013% |
| 0.36 / 0.44 / 3.5 | 138 | 0.053 | +5.2% | +3.3% | +1.19% | +1.68% | -0.271 | -0.174 | +0.011% |
| 0.30 / 0.26 / 2.5 | 135 | 0.065 | +4.2% | +2.7% | +0.98% | +1.38% | -0.182 | -0.118 | +0.009% |

ΔR_th: change of the plan-view network R_th from the developed value. Thermal entrance only (parabolic inlet):

| Design | h increase, resolved | h increase, 1D (Graetz) | ΔR_th resolved (mK/W) | ΔR_th 1D (mK/W) | energy balance |
|---|---|---|---|---|---|
| 1.00 / 0.44 / 4.0 | +25.0% | +13.4% | -4.526 | -2.666 | -0.050% |
| 0.50 / 0.44 / 4.0 | +7.7% | +4.6% | -0.570 | -0.351 | -0.014% |
| 0.42 / 0.44 / 4.0 | +5.9% | +3.6% | -0.343 | -0.215 | -0.013% |
| 0.36 / 0.44 / 3.5 | +5.2% | +3.2% | -0.270 | -0.171 | -0.011% |
| 0.30 / 0.26 / 2.5 | +4.2% | +2.7% | -0.181 | -0.116 | -0.010% |

Plenum to plenum, in units of the channel dynamic pressure ρU²/2:

| Design | σ | ρU²/2 (Pa) | total excess over developed, resolved | of which developing channel | of which inlet + outlet | 1D: f_app excess | 1D: Kc + Ke (flat duct) | 1D total |
|---|---|---|---|---|---|---|---|---|
| 1.00 / 0.44 / 4.0 | 0.693 | 11.4 | 0.59 | 0.69 | -0.10 | 0.96 | 0.43 | 1.39 |
| 0.50 / 0.44 / 4.0 | 0.531 | 19.5 | 0.83 | 0.70 | +0.14 | 0.98 | 0.70 | 1.68 |
| 0.42 / 0.44 / 4.0 | 0.488 | 23.1 | 0.91 | 0.70 | +0.21 | 0.98 | 0.78 | 1.76 |
| 0.36 / 0.44 / 3.5 | 0.445 | 36.2 | 0.96 | 0.70 | +0.26 | 0.99 | 0.86 | 1.84 |
| 0.30 / 0.26 / 2.5 | 0.537 | 48.7 | 0.83 | 0.70 | +0.13 | 0.99 | 0.69 | 1.68 |

Mesh check, recommended design, 8 → 16 cells across the half channel: combined entry h +0.008%, ΔP -0.0011% (entrance increments 5.89% → 5.90% on h and 1.340% → 1.339% on ΔP); thermal entrance h +0.0002%; plenum ΔP +0.0047%.

**Gates.** mesh change at 2x refinement, recommended design: 0.01% (limit 2%), pass; energy balance, all resolved solves: 0.050% (limit 1%), pass; 1D vs resolved, entrance contribution to R_th (recommended): -36.3% (limit ±10%), **FAIL**; 1D vs resolved, entrance contribution to ΔP (recommended): +40.8% (limit ±10%), **FAIL**.

**Diagnosis of the disagreement (not tuned).**
- *Heat transfer.* In the parallel-plate limit the 1D option predicts a 3.7% rise in h against 5.9% resolved, so it understates the R_th reduction by 36%. The cause is the developed limit the correction is scaled by: the source's Eq. 31 at aspect ratio 0.01 gives a developed Nu 19% above the exact parallel-plate value (8.235 on 2w), outside its stated ±10%, and the blend with an over-large developed value hides part of the entrance gain. Blending the source's own entrance asymptotes with the exact developed value instead gives 6.0% (resolved 5.9%); the same holds on all five designs. For the real duct (aspect ratio about 0.1) the source's developed Nu is within about 2% of the Shah & London Nu_H1 used by the screen, so this particular error is small there; but no resolved check of the real duct exists, and the 1D real-duct factor (7.6%) is larger than either plate-limit value.
- *Pressure drop.* The source's n = 2 blend has a long-duct tail equivalent to an incremental pressure-drop number 2·3.44²·(√A/Dh)/(f·Re_√A), about 0.98 for parallel plates, while the resolved developing channel gives 0.69–0.70 on all five designs. So the 1D option overstates the hydrodynamic increment by about 41%. The source's ±10% accuracy is on f_app·Re as a whole, which over 60 mm is dominated by the developed part; the increment itself is not covered by it.
- *Inlet and outlet.* Plenum to plenum, the resolved excess over developed flow is 0.91 ρU²/2 for the recommended design, against 1.76 from the 1D option (f_app excess plus Kc + Ke). The resolved inlet + outlet part (+0.21) is below the Kays coefficients (0.78) by a nearly constant 0.53–0.59 across σ = 0.44–0.69, close to the developing-channel excess. That pattern fits laminar Kc values that already contain part of the entrance pressure defect, so adding them to f_app counts it twice; Kays's paper was not available to confirm how the laminar curves were derived, so this stays a hypothesis.
- None of this changes the answer: all estimates lower R_th and raise ΔP by a few percent, the objective is met with either set of ratios, and ΔP stays far inside the 5 kPa budget.

## Not covered
- Flow maldistribution between channels, manifold, port and fitting losses (ASSUMPTIONS P5, P6).
- Axial conduction in the copper base, which moves heat toward the inlet where h is higher. With it the effective h lies between the harmonic length-mean used here and the arithmetic mean, so the R_th reductions here are the smaller estimate.
- Temperature-dependent viscosity (P3), and the 3D development of the real duct (floor, cover, corners).

Reproduce: `uv run python -m level2.entrance` (about 5 minutes); `--plot-only` redraws the figure and this file from `results.json`.
