#!/usr/bin/env python3
"""
02_rms_fields.py: where does the flow fluctuate?

The RMS field is the root-mean-square of the fluctuation at every grid point,
taken after we remove a least-squares line in time so that a slow residual
drift does not get counted as unsteadiness. Bright regions mean strong unsteady
activity, and for a forced or instability flow this lights up the wave packet
inside the shear or boundary layer.

Outputs per field:
    figures/rms_<field>.png      RMS map (with the time-mean alongside)
    tables/rms_summary.csv       min, mean, max, and p99 of each RMS field

Usage:
    python scripts/02_rms_fields.py
    python scripts/02_rms_fields.py --detrend mean
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dmdkit import dataset, fields, viz, io_utils

TUT = Path(__file__).resolve().parents[1]
EXAMPLE = TUT / "examples" / "HypersonicFlowOverPlate"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=EXAMPLE / "data")
    ap.add_argument("--out", type=Path, default=EXAMPLE / "results" / "02_rms")
    ap.add_argument("--field", default="all")
    ap.add_argument("--detrend", default="linear", choices=["linear", "mean"])
    args = ap.parse_args()
    viz.configure()
    dirs = io_utils.ensure_dirs(args.out)

    meta = dataset.load_metadata(args.data)
    keys = meta.field_keys() if args.field == "all" else [args.field]

    rows = []
    for key in keys:
        f = dataset.load_field(args.data, key)
        rms = fields.rms_field(f.data, detrend=args.detrend)
        stats = fields.field_stats(rms)
        print(f"[{key}] {f.symbol}'_rms ({args.detrend}-detrended): "
              f"mean {stats['mean']:.3e}, max {stats['max']:.3e} {f.units}")
        viz.plot_field(rms, f.x_mm, f.y_mm,
                       title=rf"$\,{f.symbol}'_{{\rm rms}}$  ({f.name}, {args.detrend} detrend)",
                       units=f.units, out_path=dirs["figures"] / f"rms_{key}.png",
                       cmap="magma", vmin=0.0)
        rows.append({"field": key, "name": f.name, "units": f.units,
                     "detrend": args.detrend, **stats})

    io_utils.write_csv(dirs["tables"] / "rms_summary.csv", rows)
    print(f"[done] {args.out}")


if __name__ == "__main__":
    main()
