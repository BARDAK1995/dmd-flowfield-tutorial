#!/usr/bin/env python3
"""
00_inspect_data.py -- understand how the data is structured and dimensionalized.

Run this first.  It answers, for any case:

    * What is in the dataset?  (fields, shapes, units)
    * How does a snapshot index become a physical time and a frequency?
        dt = snapshot spacing [s];  fs = 1/dt;  Nyquist = fs/2;  df = fs/N.
        => this is *the* reason DMD/PSD frequencies come out in kHz.
    * Is the number-density field truly dimensional?  (cross-check vs. n_inf).

It also saves the time-mean field and one representative fluctuation snapshot for
each variable so you can see the steady background and the unsteady part.

Usage:
    python scripts/00_inspect_data.py                 # uses ../data and ../outputs
    python scripts/00_inspect_data.py --data DIR --out DIR
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from dmdkit import dataset, fields, viz, io_utils

TUT = Path(__file__).resolve().parents[1]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=TUT / "data")
    ap.add_argument("--out", type=Path, default=TUT / "outputs" / "00_inspect")
    args = ap.parse_args()
    viz.configure()
    args.out.mkdir(parents=True, exist_ok=True)

    meta = dataset.load_metadata(args.data)
    g, t = meta.raw["grid"], meta.raw["time"]

    print("=" * 70)
    print(meta.raw["title"])
    print("=" * 70)
    print(f"grid        : {g['ny']} x {g['nx']} cells   "
          f"x in [{g['x_min_mm']:.2f},{g['x_max_mm']:.2f}] mm, "
          f"y in [{g['y_min_mm']:.2f},{g['y_max_mm']:.2f}] mm")
    print(f"snapshots   : {t['n_snapshots']}   "
          f"steps {t['step_first']}..{t['step_first']+(t['n_snapshots']-1)*t['step_interval']} "
          f"every {t['step_interval']}")
    print()
    print("TIME -> FREQUENCY (why results are in kHz)")
    dt = meta.dt_seconds
    record_s = (t["n_snapshots"] - 1) * dt
    df = meta.fs_hz / t["n_snapshots"]
    print(f"  snapshot dt        = {dt:.4e} s")
    print(f"  sampling freq fs   = {meta.fs_hz/1e3:,.1f} kHz")
    print(f"  Nyquist  fs/2      = {meta.nyquist_hz/1e3:,.1f} kHz   (max resolvable frequency)")
    print(f"  record length      = {record_s*1e6:.2f} us")
    print(f"  PSD resolution df  = {df/1e3:.2f} kHz   (= fs / N)")
    if meta.raw.get("forcing", {}).get("frequency_hz"):
        ff = meta.raw["forcing"]["frequency_hz"]
        print(f"  >> wall actuator forces at {ff/1e3:.0f} kHz "
              f"({'below' if ff < meta.nyquist_hz else 'ABOVE'} Nyquist -- "
              f"{'resolved' if ff < meta.nyquist_hz else 'ALIASED!'})")
    print()

    report = {"grid": g, "time": t, "fs_khz": meta.fs_hz / 1e3,
              "nyquist_khz": meta.nyquist_hz / 1e3, "df_khz": df / 1e3,
              "record_us": record_s * 1e6, "fields": {}}

    for key in meta.field_keys():
        f = dataset.load_field(args.data, key)
        mean = fields.mean_field(f.data)
        rms = fields.rms_field(f.data)
        stats = fields.field_stats(mean)
        print(f"field '{key}'  [{f.name}, {f.symbol}, {f.units}]  shape={f.data.shape}")
        print(f"   mean: min {stats['min']:.3e}  mean {stats['mean']:.3e}  max {stats['max']:.3e}")
        report["fields"][key] = {"name": f.name, "units": f.units, "shape": list(f.data.shape),
                                 "mean_stats": stats}

        viz.plot_field(mean, f.x_mm, f.y_mm, title=f"Time-mean {f.name}",
                       units=f.units, out_path=args.out / f"mean_{key}.png",
                       cmap="inferno")
        # one representative fluctuation snapshot (mid-record, mean removed)
        fluct = f.data[f.n_time // 2] - mean
        viz.plot_field(fluct, f.x_mm, f.y_mm,
                       title=f"Representative {f.symbol}' fluctuation snapshot",
                       units=f.units, out_path=args.out / f"fluct_{key}.png",
                       cmap="RdBu_r", symmetric=True)

    # dimensional cross-check on number density, if present
    if "number_density" in meta.field_keys():
        f = dataset.load_field(args.data, "number_density")
        fs = meta.raw["suggested_probes"]["freestream"]
        ref = meta.raw["freestream"]["number_density_m3"]
        chk = dataset.verify_freestream_number_density(f, fs["x_mm"], fs["y_mm"], ref)
        print()
        print(f"dimensional check: freestream n at ({chk['probe_x_mm']:.1f},"
              f"{chk['probe_y_mm']:.2f}) mm = {chk['mean_value']:.3e} m^-3, "
              f"reference n_inf = {chk['reference']:.3e}  "
              f"(rel. err {100*chk['relative_error']:+.1f}%)")
        report["number_density_freestream_check"] = chk

    io_utils.write_json(args.out / "data_report.json", report)
    print(f"\n[done] figures + data_report.json in {args.out}")


if __name__ == "__main__":
    main()
