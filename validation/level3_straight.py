"""Level 3 stage-one gates (issue #10): the 3D periodic straight-fin cell (OpenFOAM) against the Level 1 cross-section CFD.

    uv run python validation/level3_straight.py          # check the committed results (seconds, no solver)
    uv run python validation/level3_straight.py --run    # re-solve both designs at three mesh levels (Docker, about 10 min)

Exits non-zero if a gate fails on the Level 1 recommended design, or if the committed 3D results were computed against a
Level 1 reference that no longer matches examples/output/gpu_coldplate/results.json.
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from level3.straight import LEVEL1, OUT, gates, main as run_main  # noqa: E402


def check(results: Path = OUT / "results.json", level1_path: Path = LEVEL1) -> int:
    data = json.loads(Path(results).read_text())
    level1 = json.loads(Path(level1_path).read_text())
    rc = 0
    if data["recommended"] != level1["recommended"]:
        print(f"STALE: 3D results are for {data['recommended']}, Level 1 now recommends {level1['recommended']}")
        rc = 1
    for did, e in data["designs"].items():
        ref = next(d["cfd"] for d in level1["designs"] if d["cfd"]["design_id"] == did)
        for k in ("r_th_k_per_w", "dp_pa"):
            if abs(e["level1_cfd"][k] / ref[k] - 1) > 1e-9:
                print(f"STALE: {did} {k} reference {e['level1_cfd'][k]} != Level 1 {ref[k]}")
                rc = 1
        g = gates({int(k): v for k, v in e["levels"].items()}, ref)  # recomputed from the saved runs, not read back
        for k, v in g.items():
            print(f"{did} {k}: {v['value']} (limit {v['limit']}) {'pass' if v['pass'] else 'FAIL'}")
            if did == level1["recommended"] and not v["pass"]:
                rc = 1
    print("all gates pass on the recommended design" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", action="store_true", help="re-solve with OpenFOAM in Docker before checking")
    a = ap.parse_args()
    if a.run and run_main([]) != 0:
        sys.exit(1)
    sys.exit(check())
