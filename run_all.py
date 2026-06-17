#!/usr/bin/env python3
"""
run_all.py -- run the whole tutorial end to end.

Executes scripts 00-05 in order with sensible defaults so you can confirm the
environment works and see every output produced in one go.

    python run_all.py                  # full run (includes GIF animations)
    python run_all.py --quick          # skip the slow animation scripts (01, 05 gifs)

Outputs land in ./outputs/.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable


def run(rel_cmd):
    cmd = [PY, *rel_cmd]
    print("\n" + ">" * 4, " ".join(str(c) for c in cmd))
    subprocess.run(cmd, cwd=HERE, check=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="skip slow animations")
    ap.add_argument("--field", default="pressure", help="field for DMD / reconstruction")
    args = ap.parse_args()

    run(["scripts/00_inspect_data.py"])
    run(["scripts/02_rms_fields.py"])
    run(["scripts/03_point_psd.py"])
    run(["scripts/04_dmd_analysis.py", "--field", args.field,
         "--x-min-mm", "40", "--x-max-mm", "100", "--y-max-mm", "2.5",
         "--display-y-max-mm", "2.0"])
    if not args.quick:
        run(["scripts/01_animate_field.py", "--step", "1", "--display-y-max-mm", "2.5"])
        run(["scripts/05_reconstruct_from_modes.py", "--field", args.field, "--pairs", "3",
             "--x-min-mm", "40", "--x-max-mm", "100", "--y-max-mm", "2.5",
             "--display-y-max-mm", "2.0"])
    print("\n[all done] see ./outputs/")


if __name__ == "__main__":
    main()
