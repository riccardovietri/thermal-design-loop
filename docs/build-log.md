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

## 2026-11-10 (Tue, later) — field plot and a Level 2 prototype
- `tools/postprocess.py` plots axial velocity and temperature on the cross-section
  (`docs/figures/fields_w0300_t0300_h3000_q200.png`). Fin root-to-tip drop is about 2.3 K at 0.3 × 3 mm.
- Design ranges checked against published practice: Boyd lists copper skived fins at 0.1–0.8 mm thick with a
  0.1 mm minimum gap; an ASME direct-to-chip paper reports 100–200 µm fins and gaps, 1–2 mm tall
  (aspect ratio 10–20). Our 0.15 mm minimum and aspect ratio ≤ 10 sit inside that. No source found for base thickness.
- Prototype `prototypes/level2_wavy_channel.py` (not part of V1): 2D planar steady Navier-Stokes, Taylor-Hood,
  Newton with Reynolds continuation, 15 s for eight solves.
  - Straight channel reproduces plane Poiseuille ΔP to three decimals at Re 20–200.
  - Sine-wave walls (0.3 mm gap, 0.3 mm amplitude, 2 mm wavelength): ΔP is 2.05× the straight-channel value at
    Re 20 and 2.31× at Re 200. Most of that is geometric (the gap normal to the flow narrows on the slopes);
    inertia adds about 13% between Re 20 and 200. No recirculation appears up to Re 200.
  - Surprise: in 2D at microchannel Reynolds numbers, wavy fins are a pressure-drop penalty with no new flow
    structure. The published enhancement (Sui et al. 2010) comes from Dean vortices, which are 3D.
    Wavy fins are therefore a Level 3 case, not a Level 2 one.

## 2026-11-10 (Tue, late) — V1 study runs end to end
- `run_study.py examples/gpu_coldplate.yaml`: 1008 designs screened in 3 ms, 10 solved with CFD in 1.6 s, about 3 s in total.
- 1D vs CFD on R_th along the front: mean |difference| 0.9%, max 2.9% (the widest, lowest-aspect-ratio channel).
- The honest reading: for straight fins in fully developed laminar flow, a correct 1D model is enough.
  The earlier 1D model (circular-tube constants, no fin efficiency) would not have been. CFD here is a check, not a discovery tool.
- Resistance floor: at 2 LPM the coolant heating term is 3.8 mK/W and the base 1.4 mK/W, so no fin design gets below about 5.2 mK/W.
  The best design in budget reaches 9.3 mK/W; flow rate and base, not fins, set what is left.
- Base thickness: literature spans 0.1–4.5 mm with no typical value; kept at 2 mm as a labelled assumption.

## 2026-11-10 (Tue, night) — requirements sense check, Level 2 solver, first slot-jet case
- ΔP budget check against published numbers: an 80 kW rack example budgets 124 kPa in total with 34 kPa in hoses and fittings;
  a 4 kW cold plate is quoted at 55 kPa for the full loop including fittings and quick-disconnects; one measured plate
  reaches 7.6 mK/W at 4 kPa. A 30 kPa budget for the channel core alone is high; 5–10 kPa is more representative.
  Flow guidance of 1.2–1.5 LPM per kW matches the 1–2 LPM range used here.
- The knee is real, not a plotting artefact: left of it, 10× more ΔP buys about 60% lower R_th (0.12 → 1.19 kPa, 39 → 15.7 mK/W);
  right of it, 22× more buys 41% (1.19 → 26.4 kPa, 15.7 → 9.3 mK/W). The bend comes from the coolant-heating and base terms,
  which fins cannot reduce.
- Level 2 solver (`level2/planar.py`): steady Navier-Stokes (Taylor-Hood, Newton with Reynolds continuation) plus energy.
  Validation on thermally developing flow between parallel plates (`validation/level2_thermal_entrance.py`):
  pressure drop exact to 4 decimals, fully developed Nu 8.243 vs 8.235 (+0.10%), local Nu within 0.64% of an independent
  finite-difference marching solution from x* = 0.002 to 0.1, energy balance −0.03%.
- First slot-jet case (0.3 mm slots, 0.6 mm gap, 6 mm pitch, 1.5 LPM): Re 30, ΔP 15 Pa, stagnation h 14.5 kW/m²K,
  mean 4.0 kW/m²K, R_th 71 mK/W. The field plot shows the stagnation zone, wall jet and a recirculation cell beside the jet
  (`docs/figures/slot_jet_example/fields.png`). Thermally it is far worse than fins: no area enhancement, and almost no pressure used.
- Scaling insight: slot Reynolds number depends only on the number of slots (Re = ρQ/μNL), and the wall-jet region behaves like a
  thin channel heated from one side, so heat transfer is set by the gap. Narrow gaps are where slot jets become competitive.

## 2026-11-10 (Tue, night) — slot-jet sweep vs straight fins
- 24 cases (slot 0.15/0.3 mm, gap 0.05–0.4 mm, pitch 2–8 mm, 1.5 LPM), 3–209 s each, all converged, energy balance within 0.01%.
  Results in `examples/output/level2_slot_jet/`.
- Thermal resistance is set almost entirely by the gap: about 57, 37, 24 and 16–17 mK/W at 0.4, 0.2, 0.1 and 0.05 mm.
  Slot width and pitch change it by a few percent.
- Pressure drop rises about 4× per doubling of pitch at fixed gap, with almost no thermal benefit, so the tightest pitch wins.
- Against the Level 1 fin front: slot jets are better below about 1 kPa (24 mK/W at 0.11 kPa vs 39 mK/W for fins at 0.12 kPa;
  16.4 at 0.85 kPa vs 18.2 at 0.79 kPa); fins are better above it (12 mK/W at 3.4 kPa, where the best jet is still 16).
  Caveats: the fin points at low ΔP are at 1.0 LPM and the jets at 1.5 LPM; every competitive jet design needs a gap of 0.05–0.1 mm,
  below the 0.15 mm feature limit applied to fins; pitches under 2 mm and gaps under 0.05 mm were not explored.

## 2026-11-10 (Tue, late night) — requirements amendment: rated flow, 6 kPa budget, 15 mK/W target
- Riccardo's sense check (SPEC Amendment 1): rated flow 1.5 LPM; core ΔP ≤ 6 kPa at rated flow; R_th ≤ 15 mK/W; recommend the
  lowest-ΔP CFD-verified design that meets both, else the lowest R_th within budget; the knee becomes secondary because it depends
  on where the front's end points are. Minimum slot-jet gap and slot width 0.1 mm. Plots restricted to the region of interest.
- Level 1, rated flow only, refined grid 12 × 8 × 7 = 672 designs (273 feasible). The 1D front within 6 kPa has only 9 points,
  so 9 designs (not 10) were solved with CFD.
- Recommended (rule `meets_target_lowest_dp`): 0.42 mm channels, 0.44 mm fins, 4.0 mm tall, 69 channels → 14.93 mK/W at 1.29 kPa,
  interface average 54.9 °C, peak 59.9 °C. The margin to the target is 0.5%; the 1D screen gives 14.79.
- Knee (secondary): 0.50 / 0.44 / 4.0 mm → 17.12 mK/W at 0.86 kPa. Lowest R_th within budget: 0.25 / 0.26 / 2.5 mm → 11.35 mK/W at 5.75 kPa.
- 1D vs CFD on R_th over the 9 designs: mean |difference| 1.1%, max 2.9%, mean signed +0.3%. ΔP within 0.03%. All confidence checks pass.
- Flow sensitivity of the recommended design (same geometry): 1.0 LPM → 0.86 kPa, 17.43 mK/W, 57.4 °C average;
  1.5 LPM → 1.29 kPa, 14.93 mK/W, 54.9 °C; 2.0 LPM → 1.73 kPa, 13.68 mK/W, 53.7 °C. The design misses the target below rated flow.
- Slot jets, 18 cases at 1.5 LPM (slots 0.1/0.2 mm, gaps 0.1/0.15/0.2 mm, pitches 1/2/4 mm), all converged in 1–85 s,
  energy balance within 0.002%. Best R_th 23.76 mK/W at 0.49 kPa (0.1 mm slot, 0.1 mm gap, 4 mm pitch). No case meets 15 mK/W;
  the best is 58% above it. R_th is set by the gap: about 24–26, 30–33 and 36–39 mK/W at 0.1, 0.15 and 0.2 mm.
- Jets vs fins at equal flow: jets win below about 0.35 kPa (24.6 mK/W at 0.03 kPa, where fins do not reach); fins match them at
  0.38 kPa (23.7 mK/W) and beat them above it, reaching the target at 1.3 kPa. With a flat base and ≥ 0.1 mm gaps, slot jets do not
  compete for this requirement; the earlier competitive cases needed 0.05 mm gaps.
- Slot-jet mesh study (0.2 mm slot, 0.1 mm gap, 2 mm pitch): 10/20/40 cells → R_th 24.3773/24.3780/24.3784 mK/W,
  ΔP 104.93/105.27/105.53 Pa; 20 → 40 changes R_th by +0.002% and ΔP by +0.25%.

## 2026-11-11 (Wed) — plan-view Level 2: fin planforms against the rated-flow requirements
- Riccardo asked for more ambitious fin topologies than straight fins (SPEC Amendment 2). Method: steady laminar 2D Navier-Stokes +
  energy in the plane of the fin layer (tall-fin limit), streamwise-periodic cell, developed dp/dx and h, put into the Level 1
  resistance network. New code: `level2/planforms.py` (gmsh builders), `level2/planview.py` (solve + network + figures),
  `level2/planview_sweep.py`, `validation/level2_planview_straight.py`; `level2/planar.py` gained the periodic solvers.
- Validation (gates): G1 straight channel f·Re = 24.0000 (exact to 1e-6) and Nu on 2w = 8.2353 (+0.004% vs 8.235);
  G2 energy balance 2e-15 on the straight cell; G3 Level 1 recommended design 0.422 / 0.444 / 4.0 mm: plan view ΔP 1.20 kPa vs
  Level 1 CFD 1.29 kPa (−7.1%), R_th 14.08 vs 14.93 mK/W (−5.7%). Both low, as expected without floor and cover.
- Why a periodic cell: a developing strip from a uniform inlet shows how slowly the thermal field settles at Pr 11.5. Periods
  needed to bring h within 2% of the periodic value: pins 3–4 (4–5 mm), offset strips 5 (10 mm), sinusoidal and zigzag 10 (20 mm),
  straight 20–25 mm. The last two developing periods agree within 0.2% on h and dp/dx and sit within 1% of the periodic cell
  (diamond +1.0%, others ≤ 0.5%). For straight fins over 60 mm the length-mean h is 8% above the developed value, so the
  developed values are conservative, as in Level 1.
- Mesh study (10 vs 20 cells across the reference gap, one design per topology): dp/dx ≤ 0.43%, h ≤ 0.55%. Energy balance ≤ 0.75%
  on every run (sharp-cornered diamonds worst). Thermal undershoot ≤ 0.3% of the wall-to-bulk difference, so no oscillation problem.
- Sweep at 1.5 LPM, 63 designs (12 straight, 10 sinusoidal, 10 zigzag, 10 offset strip, 11 circular pins, 10 diamond pins),
  all inside the limits, all converged, 231 s of solver time (Stage 2 checks 347 s).
- Best R_th within 6 kPa: circular pins 9.79 mK/W at 5.62 kPa (d 0.4, S_T 0.6, S_L 0.7, H 4); straight 10.48 at 4.35 (w 0.25);
  offset strips 11.36 at 5.66; zigzag 11.45 at 3.06; sinusoidal 11.54 at 3.24; diamond pins 10.38 at 5.30. Every topology meets 15 mK/W.
- At equal ΔP against the straight plan-view front (log-log interpolation): circular pins 1–7% lower R_th (−3.1% at 0.80 kPa,
  −3.8% at 1.86, −7.4% at 3.21, −3.3% at 5.62); zigzag +3%, sinusoidal +5%, diamond pins +2 to +8%, offset strips +2 to +12%.
  Nothing beats straight fins by more than 7.4%, so no "too good" alarm; the pin advantage is inside the plan view's own
  uncertainty (the tall-fin bias was measured on straight fins only).
- Why the others lose: in 2D, wavy channels at Re_Dh 90–150 have no Dean vortices and almost no separation, so they add length and
  pressure drop for Nu slightly below 8.235 (7.5); offset strips are held back by the 0.15 mm gap between rows, which forces in-row
  channels of 0.7 mm or more, half the fin density of straight fins at the same gap; pins win on h (8–16 kW/m²K) despite fin
  efficiencies of 0.38–0.55 at H = 4 mm.
- Recommendation (`tools.pareto.recommend`, rule `meets_target_lowest_dp`): circular pins d 0.4, S_T 1.0, S_L 0.9 mm, H 4 mm →
  14.90 mK/W at 0.80 kPa (Re_Dh 216, Re_d,max 50). Margin to target 0.7%; with the −5.7% tall-fin bias seen on straight fins it
  would miss. Straight fins in the same model: w 0.4 / t 0.4 → 13.36 mK/W at 1.31 kPa.

## 2026-11-11 (Wed, late) — requirements amendment 3: threshold and objective, 5 kPa, marginal-return rule
- Riccardo's changes (SPEC Amendment 3): threshold R_th ≤ 20 mK/W (must meet) and objective ≤ 15 mK/W (want) replace the single
  15 mK/W target; channel-core ΔP ≤ 5 kPa (one third of a 15 kPa plate allocation); recommendation by marginal return with cut-off
  0.25 replaces "lowest ΔP meeting the target" and the knee; capacity and pumping ratios reported; H / gap ≤ 10 applied to the plan view.
  Basis for the numbers in docs/literature.md (second pass, full-text OCP and Ortega et al. sources).
- New code: `tools.pareto.recommend_marginal_return` (pure function; returns index, rule, steps with returns, front, and whether the
  threshold and objective are met) and `tools.pareto.capacity_ratios`; `recommend` and `knee_index` kept but unused in reports.
  `tools.report.draw_requirements` draws the budget, threshold and objective lines with the two windows, used by all three comparison plots.
- Level 1 at 5 kPa: 672 designs screened, 258 feasible; the 1D front within budget has 8 points, so 8 designs were solved with CFD.
  1D vs CFD on R_th: mean |difference| 1.1%, max 2.9%, mean signed +0.4%. All confidence checks pass.
- Rule steps on the Level 1 CFD front (eligible front: 0.50, 0.42, 0.36, 0.30 mm channels):
  from 0.50 / 0.44 / 4.0 mm (17.12 mK/W at 0.86 kPa): returns 0.331 to 0.42 mm, 0.222 to 0.36 mm, 0.207 to 0.30 mm → move to 0.42 mm;
  from 0.42 / 0.44 / 4.0 mm (14.93 mK/W at 1.29 kPa): returns 0.141 to 0.36 mm, 0.159 to 0.30 mm → stop.
  Recommended: 0.42 / 0.44 / 4.0 mm, 14.93 mK/W at 1.29 kPa (the same design as under Amendment 1), meets the objective,
  interface 54.9 °C average. Capacity at the 20 K rise 1.34 kW, 1.04 kW per kPa, 30.9 kW per pumping W.
- Reference designs: lowest ΔP meeting the threshold 0.50 / 0.44 / 4.0 mm (17.12 mK/W, 0.86 kPa); lowest ΔP meeting the objective is the
  recommended design; lowest R_th within 5 kPa 0.30 / 0.26 / 2.5 mm (12.58 mK/W, 3.79 kPa). The 0.60 mm design (20.05 mK/W at 0.57 kPa)
  misses the threshold by 0.3%; had it qualified, its return to 0.50 mm is 0.39, so the walk would end at the same design.
- Slot jets (plot redrawn, no re-solve): best 23.76 mK/W at 0.49 kPa; no case meets the 20 mK/W threshold.
- Plan view, H / min gap ≤ 10 (gap = true minimum distance between solids): of the 63 original designs, 36 break it (straight 6,
  sinusoidal 4, zigzag 4, offset strip 7, circular pins 7, diamond pins 8). Offset strips (3 eligible), circular pins (4) and diamond
  pins (2) fell below 5, so 15 designs were added at H = 2–4 mm with gaps ≥ H/10 (`--extend` solved only those, all converged, 1–11 s
  each, energy balance ≤ 0.49%). Eligible now: 6, 6, 6, 8, 8, 8 of 78 converged.
- Per-topology picks (marginal return, plan-view numbers): straight w 0.4 / t 0.4 / H 4 → 13.36 mK/W at 1.31 kPa (one step from w 0.5,
  return 0.283); sinusoidal w 0.4, amp 0.3, λ 4 → 13.58 at 1.54 (single-point eligible front); zigzag w 0.4, amp 0.3, λ 4 → 13.44 at 1.46;
  offset strip l 2.0 / t 0.4 / p_t 1.6 → 18.36 at 0.50 (stops at its first point: every return 0.19–0.23; misses the objective);
  circular pins d 0.4 / S_T 0.8 / S_L 0.9 → 12.65 at 1.42 (two steps); diamond pins s 0.4 / S_T 1.0 / S_L 0.7 → 12.28 at 2.51.
- Cross-topology pick: circular pins d 0.4, S_T 0.8, S_L 0.9, H 4 → 12.65 mK/W at 1.42 kPa, meets the objective; path
  offset strip (0.50 kPa) → straight w 0.5 (0.75 kPa, e 0.404) → pins S_T 1.0 / S_L 0.9 (0.80 kPa, e 0.822) → S_T 0.9 / S_L 0.9
  (1.03 kPa, e 0.309) → S_T 0.8 / S_L 0.9 (1.42 kPa, e 0.264) → stop (best further 0.230). The 0.822 step joins two points 6% apart in ΔP;
  the nearest-qualifying rule would have taken that step at any return above 0.25, so it does not change the pick.
- Circular pins against the eligible straight front at equal ΔP: −3.1% (0.80 kPa), −3.7% (1.03), −2.6% (1.37), −4.1% (1.42), −5.8% (1.86).
  The eligible straight front ends at 3.10 kPa (w 0.3 / t 0.3 / H 3), so there is no comparison above that. Best within 5 kPa:
  pins 10.61 mK/W (d 0.3, H 3, 4.19 kPa) vs straight 11.63 (3.10 kPa). Wavy fins +2.4 to +4.3%, diamond pins +2 to +7%, offset strips
  +2 to +12%. The pin advantage remains inside the plan view's own 6% tall-fin bias.
