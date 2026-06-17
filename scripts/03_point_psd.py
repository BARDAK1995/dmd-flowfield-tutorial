#!/usr/bin/env python3
"""
03_point_psd.py: spectra at a boundary point versus a freestream point.

Picks two probes (taken from the metadata, or overridden on the command line):
    * boundary  : just downstream of the actuator, inside the shear/boundary layer
    * freestream: top of the domain, in the undisturbed stream

and compares their power-spectral densities. The contrast is the whole point:

    boundary   gives a sharp TONE at the instability/forcing frequency
    freestream gives a flat broadband FLOOR (here, DSMC statistical scatter)

For each field it overlays the two PSDs, marks the known forcing frequency,
and reports the peak frequency along with the band and total RMS recovered
from the spectrum.

Usage:
    python scripts/03_point_psd.py
    python scripts/03_point_psd.py --boundary 60.8 0.5 --freestream 30 2.8 --method welch
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from dmdkit import dataset, psd, viz, io_utils

TUT = Path(__file__).resolve().parents[1]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=TUT / "data")
    ap.add_argument("--out", type=Path, default=TUT / "outputs" / "03_psd")
    ap.add_argument("--field", default="all")
    ap.add_argument("--method", default="periodogram", choices=["periodogram", "welch"])
    ap.add_argument("--boundary", type=float, nargs=2, default=None, metavar=("X_MM", "Y_MM"))
    ap.add_argument("--freestream", type=float, nargs=2, default=None, metavar=("X_MM", "Y_MM"))
    ap.add_argument("--fmax-khz", type=float, default=None, help="upper limit of the x-axis")
    args = ap.parse_args()
    viz.configure()
    dirs = io_utils.ensure_dirs(args.out)

    meta = dataset.load_metadata(args.data)
    probes = meta.raw.get("suggested_probes", {})
    b = args.boundary or (probes["boundary"]["x_mm"], probes["boundary"]["y_mm"])
    fsm = args.freestream or (probes["freestream"]["x_mm"], probes["freestream"]["y_mm"])
    forcing_khz = (meta.raw.get("forcing", {}).get("frequency_hz") or 0.0) / 1e3
    fmax_khz = args.fmax_khz or meta.nyquist_hz / 1e3
    keys = meta.field_keys() if args.field == "all" else [args.field]

    rows, summary = [], []
    for key in keys:
        f = dataset.load_field(args.data, key)
        curves = []
        for label, (px, py), color, ls in (
            ("boundary", b, "#c1121f", "-"),
            ("freestream", fsm, "#1f77b4", "-"),
        ):
            s = f.point_series(px, py)
            freq_hz, pxx = psd.point_psd(s["series"], f.dt_seconds, method=args.method)
            f_khz = psd.to_khz(freq_hz)
            curves.append({"f_khz": f_khz, "psd": pxx,
                           "label": f"{label}  ({s['x_actual_mm']:.1f}, {s['y_actual_mm']:.2f}) mm",
                           "color": color, "linestyle": ls})
            peak = psd.peak_frequency_hz(freq_hz, pxx, f_min_hz=10e3)
            total_rms = psd.rms_from_psd(freq_hz, pxx)
            print(f"[{key}/{label}] peak {peak/1e3:7.1f} kHz   total RMS {total_rms:.3e} {f.units}")
            summary.append({"field": key, "probe": label,
                            "x_mm": s["x_actual_mm"], "y_mm": s["y_actual_mm"],
                            "peak_khz": peak / 1e3, "total_rms": total_rms,
                            "mean_value": float(np.mean(s["series"]))})
            for fk, pv in zip(f_khz, pxx):
                rows.append({"field": key, "probe": label,
                             "frequency_khz": f"{fk:.6f}", "psd": f"{pv:.8e}"})

        viz.plot_psd(curves, dirs["figures"] / f"psd_{key}.png",
                     units_symbol=f.symbol, units=f.units,
                     title=f"PSD of {f.name}: boundary vs. freestream",
                     mark_khz=[forcing_khz] if forcing_khz else [],
                     xlim_khz=(0.0, fmax_khz))

    io_utils.write_csv(dirs["tables"] / "psd_data.csv", rows)
    io_utils.write_csv(dirs["tables"] / "psd_summary.csv", summary)
    print(f"[done] {args.out}")


if __name__ == "__main__":
    main()
