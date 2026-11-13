# Work plan: next issues

Scoped instructions for whoever picks up the open issues, human or agent. Read `AGENTS.md` first; its rules apply to everything here.
This file says what to do, in what order, and when to stop. The issues hold the same acceptance criteria in short form.

## Ground rules for this phase
- **One owner per issue.** Comment "taking this" on the issue before starting and say when you stop. Do not work on an issue someone else has claimed.
- **One issue, one branch, one pull request.** Branch from an up-to-date `main`. Merge only when CI passes. Do not bundle issues.
- **The frozen study does not change.** `SPEC.md` (with Amendments 1–3), the requirement numbers, the selection rule and the committed
  results under `examples/output/` stay as they are unless the issue says otherwise. New work adds results beside them.
- **No number without a source or a derivation.** A correlation or property value goes into `ASSUMPTIONS.md` with a source that was
  actually read. Do not use constants from memory.
- **Never loosen a gate, a test or a confidence limit to pass.** If a gate fails, stop and report the numbers.
- **Physics-judgment calls go to Riccardo.** That means: thermal boundary conditions, what counts as the reference result, any change to
  requirements, and anything that would change the recommended design. Plumbing is yours to decide; log it in `DECISIONS.md`.
- **Report honestly.** Each pull request states what was verified, what was not, and anything that looked physically wrong.
- **Machine limits.** Apple M3 Pro, 11 cores, 18 GB RAM. Everything must run locally. Docker is acceptable for a solver with the reason logged.

## Order
1. Issue #12, entrance effects. Small, uses the existing solvers, and addresses the thinnest margin in the study.
2. Issue #10, 3D straight fins against circular pins. The main piece. Two stages with a hard gate between them.
3. Issue #11, validation against a measured plate. Can run in parallel with either; starts as a literature search.
4. Issue #14, command-line wrappers. Any time.
5. Issue #13, selection-rule clause. **Blocked: needs Riccardo's decision. Do not start.**

## Definitions to keep fixed
- **R_th** = (area-average temperature of the plate base at the chip interface − coolant inlet temperature) / heat load. Includes base
  conduction, fins, convection and coolant heating. Excludes interface material and the chip.
- **ΔP** = channel-core pressure drop over the 60 mm flow length. No manifolds, ports or fittings unless a result says so explicitly.
- **Operating point:** 1 kW uniform over 60 × 60 mm, 25% propylene glycol at 40 °C inlet with constant properties from `tools/fluids.py`,
  1.5 LPM total, copper k = 390 W/m·K, 2 mm base.
- **Comparisons between designs are at equal total flow and equal footprint**, then read at equal ΔP. Never at equal Reynolds number.

## Issue #12: entrance effects in the straight-fin study
**Question.** How much do developing flow and inlet/outlet losses move R_th and ΔP for the recommended design, and does it still meet the
15 mK/W objective?

**Do.**
1. Find a published developing-flow treatment for laminar flow in rectangular ducts (hydrodynamic and thermal, uniform-flux-type boundary)
   and for contraction and expansion losses. Read the source; record it in `docs/literature.md` and `ASSUMPTIONS.md`.
2. Add it to the 1D screen as an option that is off by default, so existing results reproduce unchanged. Tests against values from the source.
3. Check it with a resolved solve on the recommended design. The plan-view developing strip in `level2/planview.py` already does the
   in-plane part; state clearly that it is the tall-fin limit.
4. Report R_th and ΔP for the recommended design and the two neighbours on the front, with and without entrance effects, in a short
   section added to the example report or a new `examples/output/entrance_effects/summary.md`.

**Done when.** Tests pass; the with/without table exists; `ASSUMPTIONS.md` row P2 is updated; the write-up says whether the objective is still met.
**Do not.** Change the default study, the recommendation or the committed Level 1 results.

## Issue #10: 3D check of straight fins against circular pins
**Question.** The plan-view model puts circular pins 3–6% below straight fins at equal ΔP. The same model reads 5.7% low on R_th and 7.1% low
on ΔP for straight fins because it has no floor or cover. Does the pin lead survive in 3D?

### Stage 1: toolchain and validation gate
1. Choose a 3D conjugate solver that runs on the machine above. Candidates: OpenFOAM in Docker, FEniCSx. The existing scikit-fem direct
   solves will not scale to 3D here. Log the options and the choice in `DECISIONS.md`.
2. Build a 3D model of the **straight-fin** unit cell: one fin pitch wide with symmetry planes, full fin height with the floor, the
   adiabatic cover, the 2 mm base with uniform flux underneath, conjugate conduction in fin and base. Use a streamwise-periodic
   (fully developed) formulation so it is comparable with Level 1; a short developing run is optional.
3. Run the Level 1 recommended design (0.42 mm channel, 0.44 mm fin, 4.0 mm tall).

**Gate.** R_th and ΔP within 2% of the Level 1 cross-section CFD (14.93 mK/W, 1.29 kPa), mesh change under 2% at one refinement, energy
balance within 1%. **If the gate fails, stop and report.** Do not proceed to pins on an unvalidated model.

### Stage 2: pins
1. Circular-pin unit cell with the same boundaries: floor, cover, base, conduction along the pin. Streamwise-periodic, lateral symmetry.
2. Designs: the cross-topology pick (0.4 mm diameter, 0.8 × 0.9 mm pitch, 4 mm tall) and its two neighbours on the pin front in
   `examples/output/level2_planview/results.json`; for straight fins, the plan-view pick (0.4 / 0.4 / 4.0 mm) and the Level 1 pick.
3. Steadiness: at pin Reynolds numbers up to about 140 the wake may shed. Run a transient or otherwise show the steady solution is
   stable. If it sheds, report time-averaged values and say so.
4. Mesh independence and energy balance for at least one pin case.

**Deliver.** `level3/` code with a validation script; results under `examples/output/level3_pins_vs_straight/`; one figure (R_th against
ΔP for both, 3D and plan view) in the shared plot style; a table at two or three equal-ΔP levels; a paragraph answering the question;
`ASSUMPTIONS.md`, `DECISIONS.md`, `README.md` and `PROGRESS.md` updated. CI should run the validation only if it is fast; otherwise
commit its output and document the command.
**Do not.** Add fin shapes, add manifolds, or change the plan-view results. If 3D contradicts the plan view, report it; do not tune either model.

## Issue #11: validation against a measured cold plate
**Question.** How far is the model from a real plate?

**Do.**
1. Find a published single-phase test of a straight-fin copper cold plate with **complete** data: fin thickness, channel width, fin
   height, flow length, plate or finned width, base thickness, coolant and inlet temperature, heat load and heated area, where the
   temperatures were measured, and R_th (and ideally ΔP) against flow. Ortega et al. (ASME InterPACK 2022, in `docs/literature.md`)
   lacks width and base thickness in its text; check its figures and cited CFD setup before discarding it.
2. Use scholarly APIs and read the full text. Mark every number as read from text, read from a figure, or inferred.
3. Run Level 0 and Level 1 at the test conditions. Compare R_th over the measured flow range. Compare ΔP only if manifold losses can be
   separated; otherwise say it cannot be compared.

**Done when.** A script under `validation/` with the reference data and citation, a model-against-measurement plot, and a stated error.
If no complete data set exists, a short note listing what was searched and what each candidate lacked is an acceptable result.
**Do not.** Fill missing geometry with guesses and present the result as validation. A sensitivity band over a stated assumption is fine if labelled.

## Issue #14: command-line entry points
Thin `argparse` wrappers for `run_cfd`, `check_confidence`, `pareto` and `report`, JSON in and out, one test each. No behaviour changes.

## How to hand back
In the pull request: what was done, the evidence for each acceptance point, what was not done, deviations from this plan with the reason,
and anything that looked physically suspicious. Update `PROGRESS.md` ("Resume here") in the same pull request.
