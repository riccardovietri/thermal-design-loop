"""State of the detached Level 3 pin jobs (mesh gate, transient) and their latest monitored values.

    uv run python -m level3.status                      # running / finished / failed, latest R_th, ΔP, time or iteration
    uv run python -m level3.pins_study --summary-only   # post-process finished jobs into examples/output/level3_pins_vs_straight/

Jobs are listed in runs/level3/long/jobs.json (written by `python -m level3.pins_study --launch-long`); each runs in its own
Docker container (`docker ps -a --filter name=l3pins_`), logs in its case directory.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

LONG = Path("runs/level3/long")


def _tail(path: Path, n_bytes: int = 400_000) -> str:
    if not path.exists():
        return ""
    with open(path, "rb") as fh:
        fh.seek(0, 2)
        size = fh.tell()
        fh.seek(max(0, size - n_bytes))
        return fh.read().decode(errors="replace")


def _last_line(path: Path) -> list[str] | None:
    if not path.exists():
        return None
    lines = [ln for ln in _tail(path, 20_000).splitlines() if ln.strip() and not ln.startswith("#")]
    return lines[-1].split() if lines else None


def latest(case: Path, log: Path) -> dict | None:
    """Latest R_th and ΔP from the tail of a solver log and the base-temperature monitor."""
    from .pins import pin_cell, r_th_pins
    from .pins_study import setup
    text = _tail(log)
    times = re.findall(r"^Time = (\S+)", text, re.M)
    grads = re.findall(r"[Pp]ressure gradient = (\S+)", text)
    bulks = re.findall(r"L3MON bulkT (\S+)", text)
    if not (times and grads and bulks):
        return None
    s = json.loads((case / "setup.json").read_text())
    _, req, fluid = setup()
    cell = pin_cell(s["params_mm"], req, fluid, full_pitch=s["full_pitch"])
    base = None
    for name in ("baseTt", "baseT"):
        files = sorted((case / "postProcessing/solid" / name).glob("*/surfaceFieldValue.dat"))
        if files:
            base = _last_line(files[-1])
            break
    out = {"time_or_iteration": times[-1], "dp_pa": round(cell.dp_core_pa(float(grads[-1])), 3)}
    if base:
        out["r_th_mk_per_w"] = round(1e3 * r_th_pins(cell, float(base[1]), float(bulks[-1])), 5)
    probe = sorted((case / "postProcessing/wake/fluid").glob("*/U"))
    if probe:
        row = _last_line(probe[-1])
        if row:
            vals = [float(v) for v in " ".join(row).replace("(", " ").replace(")", " ").split()]
            out["wake_uy_m_s"] = [f"{vals[1 + 3 * k + 1]:.2e}" for k in range(3)]
    return out


def job_state(name: str, spec: dict) -> dict:
    case = Path(spec["case"])
    r = subprocess.run(["docker", "inspect", "-f", "{{.State.Status}} {{.State.ExitCode}}", spec["container"]],
                       capture_output=True, text=True)
    status, code = (r.stdout.split() + ["?", "?"])[:2] if r.returncode == 0 else ("missing", "?")
    log = case / "log.solver"
    phase = None
    if name == "transient" and not (case / "transient.started").exists():
        log, phase = case / "log.solver.steady", "steady start"
    ended = bool(re.search(r"^End\s*$", _tail(log, 5_000), re.M)) if log.exists() else False
    if status == "running":
        state = "running"
    elif status == "exited" and code == "0" and ended and phase is None:
        state = "finished"
    elif status == "missing" and ended and phase is None:
        state = "finished"  # container removed after a complete run
    else:
        state = "failed"
    out = {"state": state, "docker": f"{status} (exit {code})", "case": str(case), "log": str(log)}
    if phase:
        out["phase"] = phase
    try:
        out["latest"] = latest(case, log)
    except Exception as exc:  # status must never crash on a half-written file
        out["latest_error"] = repr(exc)
    return out


def main() -> int:
    jf = LONG / "jobs.json"
    if not jf.exists():
        print(f"no detached jobs: {jf} not found")
        return 1
    jobs = json.loads(jf.read_text())
    for name, spec in jobs.items():
        st = job_state(name, spec)
        print(f"{name}: {st['state']} [{st['docker']}]" + (f", phase: {st['phase']}" if "phase" in st else ""))
        print(f"    {spec['what']}; started {spec['started']}; case {st['case']}")
        if st.get("latest"):
            print("    latest: " + ", ".join(f"{k} {v}" for k, v in st["latest"].items()))
        if st.get("latest_error"):
            print("    could not read monitors: " + st["latest_error"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
