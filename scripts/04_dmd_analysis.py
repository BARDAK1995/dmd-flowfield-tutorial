#!/usr/bin/env python3
"""
04_dmd_analysis.py: Dynamic Mode Decomposition, reported in physical units.

Runs the full DMD pipeline on one field over the (already converged) window and
writes a tidy, self-describing bundle:

    figures/svd_distribution.png          singular-value decay (how many modes matter)
    figures/field_overview.png            time-mean plus one detrended snapshot
    figures/dmd_frequency_spectrum.png    normalized amplitude vs. frequency [kHz]
    modes/mode_XX.png                     real-part and amplitude map of each top mode
    modes/mode_XX.npz                     raw 2-D mode fields (re-plot without re-running)
    tables/dmd_mode_summary.csv           f[kHz], growth[1/s], amplitude, lambda[mm], c[m/s], ...
    tables/svd_singular_values.csv
    analysis_metadata.json                every setting, for exact reproduction
    README.md                             human-readable digest

Each mode is reported dimensionally: frequency in kHz (from the eigenvalue and
the real dt), streamwise wavelength in mm (spatial FFT of the mode), and phase
speed c = f*lambda in m/s.

Usage:
    python scripts/04_dmd_analysis.py --field pressure
    python scripts/04_dmd_analysis.py --field number_density --svd-rank 30 --modes-to-plot 6 \
        --x-min-mm 40 --x-max-mm 100 --y-max-mm 2.5
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from dmdkit import dataset, fields, dmd, viz, io_utils

TUT = Path(__file__).resolve().parents[1]
EXAMPLE = TUT / "examples" / "HypersonicFlowOverPlate"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=EXAMPLE / "data")
    ap.add_argument("--out", type=Path, default=None, help="defaults to ../outputs/04_dmd/<field>")
    ap.add_argument("--field", default="pressure")
    ap.add_argument("--x-min-mm", type=float, default=None)
    ap.add_argument("--x-max-mm", type=float, default=None)
    ap.add_argument("--y-min-mm", type=float, default=None)
    ap.add_argument("--y-max-mm", type=float, default=None)
    ap.add_argument("--detrend", default="mean", choices=["mean", "moving_average", "none"])
    ap.add_argument("--moving-average-window", type=int, default=300)
    ap.add_argument("--svd-rank", type=int, default=30)
    ap.add_argument("--modes-to-plot", type=int, default=6)
    ap.add_argument("--display-y-max-mm", type=float, default=None)
    args = ap.parse_args()
    viz.configure()

    out = args.out or (EXAMPLE / "results" / "04_dmd" / args.field)
    dirs = io_utils.ensure_dirs(out)
    meta = dataset.load_metadata(args.data)

    f = dataset.load_field(args.data, args.field)
    if any(v is not None for v in (args.x_min_mm, args.x_max_mm, args.y_min_mm, args.y_max_mm)):
        f = f.crop(args.x_min_mm, args.x_max_mm, args.y_min_mm, args.y_max_mm)
    print(f"[dmd] field '{args.field}' ({f.name}) cube {f.data.shape}, dt={f.dt_seconds:.3e}s")

    # Fit the DMD model.
    result = dmd.run_dmd(
        f.data, f.dt_seconds, f.x_mm, f.y_mm, units=f.units, symbol=f.symbol,
        detrend=args.detrend, moving_average_window=args.moving_average_window,
        svd_rank=args.svd_rank)
    print(f"[dmd] {result.n_positive} positive-frequency modes; "
          f"dominant at {result.mode_info(1)['frequency_khz']:.1f} kHz")

    # Overview figures.
    viz.plot_svd(result.singular_values, dirs["figures"] / "svd_distribution.png",
                 rank_marker=args.svd_rank)

    mean = fields.mean_field(f.data)
    detr = fields.fluctuation_cube(f.data, detrend=("linear" if args.detrend == "moving_average" else "mean"))
    viz.plot_field(mean, f.x_mm, f.y_mm, title=f"Time-mean {f.name}",
                   units=f.units, out_path=dirs["figures"] / "field_overview_mean.png",
                   cmap="inferno")
    viz.plot_field(detr[detr.shape[0] // 2], f.x_mm, f.y_mm,
                   title=f"Representative detrended {f.symbol}' snapshot",
                   units=f.units, out_path=dirs["figures"] / "field_overview_fluct.png",
                   cmap="RdBu_r", symmetric=True)

    # Frequency spectrum.
    pos = result.pos_order
    freqs_khz = result.freq_hz[pos] / 1e3
    amps = np.abs(result.amplitudes[pos])
    amps_norm = amps / amps.max()
    viz.plot_spectrum(freqs_khz, amps_norm, dirs["figures"] / "dmd_frequency_spectrum.png",
                      title=f"Positive-frequency DMD spectrum: {f.name}",
                      highlight_n=args.modes_to_plot,
                      mark_khz=[(meta.raw.get("forcing", {}).get("frequency_hz") or 0) / 1e3]
                      if meta.raw.get("forcing", {}).get("frequency_hz") else [])

    # Export each mode.
    rows = []
    for rank in range(1, result.n_positive + 1):
        info = result.mode_info(rank)
        rows.append({
            "positive_mode_rank": rank,
            "full_index": info["full_index"],
            "frequency_khz": f"{info['frequency_khz']:.6f}",
            "growth_rate_per_s": f"{info['growth_rate_per_s']:.6e}",
            "amplitude_abs": f"{info['amplitude_abs']:.6e}",
            "normalized_amplitude": f"{info['normalized_amplitude']:.6f}",
            "coefficient_rms": f"{info['coefficient_rms']:.6e}",
            "wavelength_mm": f"{info['wavelength_mm']:.6f}",
            "phase_speed_mps": f"{info['phase_speed_mps']:.6f}",
            "eigenvalue_real": f"{info['eigenvalue_real']:.8f}",
            "eigenvalue_imag": f"{info['eigenvalue_imag']:.8f}",
            "mode_real_abs_max_units": f"{info['mode_real_abs_max']:.6e}",
            "mode_amplitude_max_units": f"{info['mode_amplitude_max']:.6e}",
        })
        if rank <= args.modes_to_plot:
            real_map = result.mode_real_map(rank)
            amp_map = result.mode_amplitude_map(rank)
            viz.plot_mode(real_map, amp_map, f.x_mm, f.y_mm, rank=rank,
                          frequency_khz=info["frequency_khz"], wavelength_mm=info["wavelength_mm"],
                          phase_speed_mps=info["phase_speed_mps"], units=f.units,
                          variable_name=f.name, out_path=dirs["modes"] / f"mode_{rank:02d}.png",
                          display_y_max_mm=args.display_y_max_mm)
            io_utils.save_mode_npz(
                dirs["modes"] / f"mode_{rank:02d}.npz", x_mm=f.x_mm, y_mm=f.y_mm,
                frequency_hz=info["frequency_hz"], growth_rate_per_s=info["growth_rate_per_s"],
                units=f.units, real_map=real_map, amplitude_map=amp_map)

    io_utils.write_csv(dirs["tables"] / "dmd_mode_summary.csv", rows)
    io_utils.write_svd_csv(dirs["tables"] / "svd_singular_values.csv", result.singular_values)

    cond = float(result.singular_values[0] / result.singular_values[-1])
    metadata = {
        "field": args.field, "variable_name": f.name, "units": f.units,
        "cube_shape": list(f.data.shape),
        "crop_mm": {"x": [float(f.x_mm[0]), float(f.x_mm[-1])],
                    "y": [float(f.y_mm[0]), float(f.y_mm[-1])]},
        "dt_seconds": f.dt_seconds, "sampling_frequency_hz": f.fs_hz,
        "step_first": f.step_first, "step_interval": f.step_interval,
        "detrend": args.detrend, "moving_average_window": args.moving_average_window,
        "svd_rank": args.svd_rank, "exact": True, "opt": True,
        "forward_backward": True, "sorted_eigs": "abs",
        "n_positive_modes": result.n_positive,
        "snapshot_matrix_condition_number": cond,
        "dominant_mode": result.mode_info(1),
    }
    io_utils.write_json(out / "analysis_metadata.json", metadata)

    top = result.mode_info(1)
    io_utils.write_readme(out / "README.md", f"DMD analysis: {f.name}", [
        f"- Field: `{args.field}` ({f.name}, {f.units})",
        f"- Window: steps {f.step_first}..{f.step_first + (f.n_time-1)*f.step_interval}, "
        f"{f.n_time} snapshots, dt = {f.dt_seconds:.3e} s (fs = {f.fs_hz/1e3:.0f} kHz).",
        f"- Spatial crop: x in [{f.x_mm[0]:.1f}, {f.x_mm[-1]:.1f}] mm, "
        f"y in [{f.y_mm[0]:.2f}, {f.y_mm[-1]:.2f}] mm.",
        f"- Detrend: {args.detrend}; DMD: forward-backward exact, rank {args.svd_rank}, "
        f"optimal amplitudes, eigenvalue-sorted by magnitude.",
        f"- Snapshot-matrix condition number: {cond:.3e}.",
        "",
        "## Method",
        "Physical frequency f = Im(log lambda)/(2*pi*dt); growth = Re(log lambda)/dt. "
        "Positive-frequency modes are ranked by |amplitude|. Each mode's streamwise "
        "wavelength comes from a spatial FFT of its real part, and phase speed c = f*lambda.",
        "",
        "## Leading modes",
        *[f"- Mode {r['positive_mode_rank']:02d}: {float(r['frequency_khz']):.1f} kHz, "
          f"norm. amp {float(r['normalized_amplitude']):.3f}, "
          f"lambda {float(r['wavelength_mm']):.2f} mm, c {float(r['phase_speed_mps']):.0f} m/s"
          for r in rows[:min(8, len(rows))]],
        "",
        f"Dominant mode: {top['frequency_khz']:.1f} kHz "
        f"(phase speed {top['phase_speed_mps']:.0f} m/s).",
    ])
    print(f"[done] {out}")


if __name__ == "__main__":
    main()
