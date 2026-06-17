#!/usr/bin/env python3
"""
05_reconstruct_from_modes.py rebuilds the flow from a few DMD modes.

The promise of DMD is that a handful of modes should capture the coherent
unsteady physics. Here we prove it by reconstructing the field from the first
few positive-frequency modes and watching it reproduce the travelling wave.

Key idea (conjugate pairs):
    Real data gives DMD eigenvalues in complex-conjugate pairs (a +f mode and its
    -f twin). To rebuild a real field you must keep both members of each pair.
    So "reconstruct from the first N modes" here means the first N
    positive-frequency modes plus their conjugates, which gives 2N complex modes
    that combine into a real field. The toolkit handles the pairing automatically.

    field(t) ~= mean  +  sum_{k in selected pairs} phi_k * b_k * lambda_k^t

Outputs:
    figures/reconstruction_error_vs_pairs.png   relative L2 error vs #pairs
    figures/snapshot_compare.png                original vs 1/2/3-pair reconstruction
    <field>_reconstruction_<N>pairs.gif         original fluctuation vs reconstruction
    <field>_fullfield_<N>pairs.gif              mean plus reconstruction (looks like the flow)
    tables/reconstruction_error.csv

Usage:
    python scripts/05_reconstruct_from_modes.py --field pressure --pairs 3
    python scripts/05_reconstruct_from_modes.py --field number_density --pairs 2 \
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=TUT / "data")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--field", default="pressure")
    ap.add_argument("--pairs", type=int, default=3,
                    help="number of conjugate mode-pairs in the showcase reconstruction (gives 2*pairs modes)")
    ap.add_argument("--ranks", type=int, nargs="*", default=None, metavar="R",
                    help="reconstruct from these EXACT positive-frequency mode ranks "
                         "(1-based, e.g. --ranks 1 3 5). Overrides --pairs for the "
                         "showcase, and conjugate twins are added automatically.")
    ap.add_argument("--max-pairs", type=int, default=6, help="largest number of pairs in the error sweep")
    ap.add_argument("--x-min-mm", type=float, default=None)
    ap.add_argument("--x-max-mm", type=float, default=None)
    ap.add_argument("--y-min-mm", type=float, default=None)
    ap.add_argument("--y-max-mm", type=float, default=None)
    ap.add_argument("--svd-rank", type=int, default=30)
    ap.add_argument("--detrend", default="mean", choices=["mean", "moving_average", "none"])
    ap.add_argument("--fps", type=int, default=15)
    ap.add_argument("--anim-max-frames", type=int, default=160,
                    help="cap GIF length (frames are subsampled to this many)")
    ap.add_argument("--display-y-max-mm", type=float, default=None)
    args = ap.parse_args()
    viz.configure()

    out = args.out or (TUT / "outputs" / "05_reconstruction" / args.field)
    dirs = io_utils.ensure_dirs(out)

    f = dataset.load_field(args.data, args.field)
    if any(v is not None for v in (args.x_min_mm, args.x_max_mm, args.y_min_mm, args.y_max_mm)):
        f = f.crop(args.x_min_mm, args.x_max_mm, args.y_min_mm, args.y_max_mm)

    result = dmd.run_dmd(f.data, f.dt_seconds, f.x_mm, f.y_mm, units=f.units,
                         symbol=f.symbol, detrend=args.detrend, svd_rank=args.svd_rank)
    print(f"[recon] field '{args.field}'; {result.n_positive} positive-frequency modes available")

    # reference fluctuation, which is what the reconstructed fluctuation should match
    fluct = fields.fluctuation_cube(
        f.data, detrend=("linear" if args.detrend == "moving_average" else "mean"))

    # --- the coherent DMD model (ALL retained modes) --------------------- #
    # A few modes can only ever rebuild the coherent part of the flow. In DSMC
    # the raw fluctuation is coherent_wave plus large incoherent statistical noise,
    # so the error against raw data plateaus at the noise floor. The meaningful
    # question is how fast the first N pairs converge to the full coherent model.
    max_pairs = min(args.max_pairs, result.n_positive)
    full_model = result.reconstruct(list(range(1, result.n_positive + 1)), add_mean=False)
    coherent_fraction = float(np.linalg.norm(full_model) / np.linalg.norm(fluct))
    print(f"[recon] coherent fraction ||DMD model|| / ||fluctuation|| = "
          f"{coherent_fraction:.2f}  (the rest is incoherent DSMC noise)")

    # --- error sweep from 1 to max_pairs --------------------------------- #
    err_rows, err_full, err_raw = [], [], []
    for n in range(1, max_pairs + 1):
        recon = result.reconstruct(list(range(1, n + 1)), add_mean=False)
        e_full = dmd.reconstruction_error(full_model, recon)
        e_raw = dmd.reconstruction_error(fluct, recon)
        err_full.append(e_full)
        err_raw.append(e_raw)
        freqs = [result.mode_info(r)["frequency_khz"] for r in range(1, n + 1)]
        err_rows.append({"n_pairs": n, "n_modes": 2 * n,
                         "error_vs_full_DMD_model": f"{e_full:.6f}",
                         "error_vs_raw_fluctuation": f"{e_raw:.6f}",
                         "frequencies_khz": ";".join(f"{x:.1f}" for x in freqs)})
        print(f"   {n} pair(s) (={2*n} modes, f={[f'{x:.0f}' for x in freqs]} kHz): "
              f"err vs full model {e_full:.3f}, vs raw {e_raw:.3f}")

    io_utils.write_csv(dirs["tables"] / "reconstruction_error.csv", err_rows)

    viz.configure()
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.6, 5.0), constrained_layout=True)
    xs = list(range(1, max_pairs + 1))
    ax.plot(xs, err_full, "o-", color="#1f4e79", label="vs. full DMD model (coherent)")
    ax.plot(xs, err_raw, "s--", color="#c1121f",
            label="vs. raw fluctuation (incl. DSMC noise)")
    ax.axhline(np.sqrt(max(0.0, 1.0 - coherent_fraction ** 2)), color="#999999", ls=":",
               label=f"noise floor (coherent frac {coherent_fraction:.2f})")
    ax.set_xlabel("number of conjugate mode-pairs")
    ax.set_ylabel("relative L2 error")
    ax.set_title(f"Reconstruction convergence: {f.name}")
    ax.grid(True, alpha=0.25)
    ax.set_ylim(0.0, 1.05)
    ax.legend()
    fig.savefig(dirs["figures"] / "reconstruction_error_vs_pairs.png", bbox_inches="tight")
    plt.close(fig)

    # --- snapshot comparison of raw, full model, and 1/2/3-pair reconstruction --- #
    k = fluct.shape[0] // 2
    panels = [("raw fluctuation (noisy)", fluct[k]),
              ("full DMD model (coherent)", full_model[k])]
    for n in (1, 2, min(3, max_pairs)):
        rc = result.reconstruct(list(range(1, n + 1)), add_mean=False,
                                time_indices=np.array([k]))[0]
        panels.append((f"{n} pair(s) = {2*n} modes", rc))
    m = float(np.quantile(np.abs(fluct[k]), 0.99))
    fig, axes = plt.subplots(len(panels), 1, figsize=(12, 2.1 * len(panels)),
                             sharex=True, constrained_layout=True)
    ext = (float(f.x_mm[0]), float(f.x_mm[-1]), float(f.y_mm[0]), float(f.y_mm[-1]))
    for ax, (lab, fld) in zip(axes, panels):
        im = ax.imshow(fld, origin="lower", aspect="auto", extent=ext, cmap="RdBu_r",
                       vmin=-m, vmax=m, interpolation="nearest")
        ax.set_ylabel("y [mm]")
        ax.set_title(f"{f.symbol}':  {lab}")
        if args.display_y_max_mm:
            ax.set_ylim(float(f.y_mm[0]), args.display_y_max_mm)
        fig.colorbar(im, ax=ax, pad=0.01).set_label(f.units)
    axes[-1].set_xlabel("x [mm]")
    fig.savefig(dirs["figures"] / "snapshot_compare.png", bbox_inches="tight")
    plt.close(fig)

    # --- showcase animations --------------------------------------------- #
    # The showcase reconstructs from either an explicit set of ranks (--ranks,
    # any subset such as 1 3 5) or the first N pairs (--pairs).
    if args.ranks:
        sel_ranks = [r for r in args.ranks if 1 <= r <= result.n_positive]
        tag = "ranks_" + "-".join(str(r) for r in sel_ranks)
    else:
        n = min(args.pairs, max_pairs)
        sel_ranks = list(range(1, n + 1))
        tag = f"{n}pairs"
    sel_freqs = [result.mode_info(r)["frequency_khz"] for r in sel_ranks]
    print(f"[recon] showcase reconstruction from ranks {sel_ranks} "
          f"= {2*len(sel_ranks)} modes at {[f'{x:.0f}' for x in sel_freqs]} kHz")

    times = f.times_s()
    astep = max(1, fluct.shape[0] // args.anim_max_frames)
    recon_fluct = result.reconstruct(sel_ranks, add_mean=False)
    viz.animate_pair(
        fluct, recon_fluct, f.x_mm, f.y_mm,
        out / f"{args.field}_reconstruction_{tag}.gif",
        title_left=f"original {f.symbol}' fluctuation",
        title_right=f"DMD reconstruction (modes {sel_ranks} = {2*len(sel_ranks)} modes)",
        suptitle=f"{f.name}: original vs. DMD reconstruction from ranks {sel_ranks}",
        units=f.units, cmap="RdBu_r", symmetric=True, fps=args.fps, step=astep,
        times_s=times, display_y_max_mm=args.display_y_max_mm)

    recon_full = result.reconstruct(sel_ranks, add_mean=True)
    viz.animate_field(
        recon_full, f.x_mm, f.y_mm, out / f"{args.field}_fullfield_{tag}.gif",
        title=f"{f.name}: full-field reconstruction (modes {sel_ranks})", units=f.units,
        cmap="inferno", symmetric=False, fps=args.fps, step=astep, times_s=times,
        display_y_max_mm=args.display_y_max_mm)

    print(f"[done] {out}")


if __name__ == "__main__":
    main()
