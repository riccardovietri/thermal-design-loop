# PROGRESS

## Resume here
Next: M2 — 1D straight-fin screen as `tools/screen_1d.py` with tests, then M3 batch loop.
Branch `feat/m1-single-design` holds M1; `main` has the frozen spec.

## Done
- 2026-11-10 (Tue): interview complete; SPEC frozen; repo created as a uv project.
- 2026-11-10: toolchain proven. `uv run python validation/duct_fully_developed.py` → f·Re and Nu_H1
  within 0.007% of Shah & London at α = 1, 0.5, 0.25, 0.125.
- 2026-11-10: M1 done a day early. `tools/mesh.py` (gmsh unit cell + quality) and `tools/run_cfd.py`
  (conjugate solve) run one design in 0.1–0.3 s; `uv run pytest` → 4 passed. Numbers in `docs/build-log.md`.

## Open items
- Level 2 (2D developing-flow Navier-Stokes) approved by Riccardo as a conditional stretch after V1 runs end to end;
  prototype in `prototypes/`. Showcase case not yet chosen (wavy fins ruled out: needs 3D).
- ΔP parity is exact by construction (both models fully developed); decide whether the hero parity plot
  is R_th only, or whether to add the stretch streamwise/entrance case.
- Energy-balance check is currently the pinned-node reaction (near machine zero by construction);
  M4 should add an independent interface heat-flux integral.
- `tools/run_cfd.py` has no CLI yet (mesh has one); add with M3.
- Dow lists cp about 4% above CoolProp for PG25; unresolved, low impact.
- thermal-mcp-server not yet a dependency (add in M2 if its resistance-network code is reused).
- Portfolio write-up: keep `docs/build-log.md` current each session; it feeds the README and site page.

## Notes for a fresh session
- Run everything with `uv run`. In Claude Code's sandbox, `uv` and `git` writes need the sandbox disabled.
- Fable usage: 7% used after M1 on 11/10.
- GitHub: private repo riccardovietri/coldplate-design-loop; `main` pushed once with Riccardo's OK, work goes on branches.
