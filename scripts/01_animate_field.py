#!/usr/bin/env python3
"""
01_animate_field.py: visualize the flow as a movie.

Two GIFs per field:
    * <field>_raw.gif     : the field itself (steady background plus waves)
    * <field>_fluct.gif   : the fluctuation only (time-mean removed), which is
                            where you actually *see* the travelling instability
                            or forced wave packet moving downstream.

The fluctuation movie is the honest picture of the unsteady physics, and it is
the thing the DMD reconstruction (script 05) will be compared against.

Usage:
    python scripts/01_animate_field.py                       # both fields
    python scripts/01_animate_field.py --field pressure --fps 20
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dmdkit import dataset, fields, viz

TUT = Path(__file__).resolve().parents[1]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=TUT / "data")
    ap.add_argument("--out", type=Path, default=TUT / "outputs" / "01_animation")
    ap.add_argument("--field", default="all", help="'all', or a single field key")
    ap.add_argument("--fps", type=int, default=15)
    ap.add_argument("--step", type=int, default=1, help="use every Nth snapshot")
    ap.add_argument("--display-y-max-mm", type=float, default=None)
    args = ap.parse_args()
    viz.configure()
    args.out.mkdir(parents=True, exist_ok=True)

    meta = dataset.load_metadata(args.data)
    keys = meta.field_keys() if args.field == "all" else [args.field]

    for key in keys:
        f = dataset.load_field(args.data, key, mmap=False)
        times = f.times_s()
        print(f"[{key}] animating {f.n_time} snapshots ({f.name}, {f.units}) ...")

        p_raw = viz.animate_field(
            f.data, f.x_mm, f.y_mm, args.out / f"{key}_raw",
            title=f"{f.name}", units=f.units, cmap="inferno",
            symmetric=False, fps=args.fps, step=args.step, times_s=times,
            display_y_max_mm=args.display_y_max_mm)

        fluct = fields.fluctuation_cube(f.data, detrend="mean")
        p_fluct = viz.animate_field(
            fluct, f.x_mm, f.y_mm, args.out / f"{key}_fluct",
            title=f"{f.symbol}' fluctuation", units=f.units, cmap="RdBu_r",
            symmetric=True, fps=args.fps, step=args.step, times_s=times,
            display_y_max_mm=args.display_y_max_mm)
        print(f"   wrote {p_raw.name} and {p_fluct.name}")

    print(f"[done] {args.out}")


if __name__ == "__main__":
    main()
