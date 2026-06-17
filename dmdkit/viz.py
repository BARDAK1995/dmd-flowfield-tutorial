"""
viz.py: plotting and animation helpers.

Everything renders through a non-interactive backend so the scripts run headless
(on a cluster, in CI, over SSH). Font sizes obey journal minimums:
labels >= 14 pt, ticks >= 12 pt, titles >= 16 pt.

Animations can fall back to GIFs written through Matplotlib's Pillow writer, which
needs no external binaries, so a collaborator can reproduce them anywhere.
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter, writers

# Prefer MP4 (ffmpeg), which is roughly 20x smaller than GIF for wide colormap movies.
# Fall back to a size-limited GIF if ffmpeg is unavailable.
_HAS_FFMPEG = "ffmpeg" in writers.list()


def _save_anim(anim, out_path: Path, fps: int, dpi: int):
    """Save an animation as an .mp4 (ffmpeg) when possible, otherwise a lean .gif."""
    out_path = Path(out_path)
    if _HAS_FFMPEG:
        mp4 = out_path.with_suffix(".mp4")
        anim.save(str(mp4), writer=FFMpegWriter(fps=fps, bitrate=3000), dpi=dpi)
        return mp4
    gif = out_path.with_suffix(".gif")
    anim.save(str(gif), writer=PillowWriter(fps=fps), dpi=min(dpi, 70))
    return gif


def configure() -> None:
    plt.rcParams.update({
        "figure.dpi": 130, "savefig.dpi": 200,
        "font.family": "DejaVu Sans", "mathtext.fontset": "stix",
        "axes.titlesize": 17, "axes.labelsize": 15,
        "xtick.labelsize": 13, "ytick.labelsize": 13,
        "legend.fontsize": 13, "figure.facecolor": "white", "axes.facecolor": "white",
    })


# --------------------------------------------------------------------------- #
# Static maps
# --------------------------------------------------------------------------- #
def plot_field(field2d: np.ndarray, x_mm, y_mm, *, title: str, units: str,
               out_path: Path, cmap: str = "inferno", symmetric: bool = False,
               vmin=None, vmax=None) -> None:
    configure()
    ext = (float(x_mm[0]), float(x_mm[-1]), float(y_mm[0]), float(y_mm[-1]))
    fig, ax = plt.subplots(figsize=(12.5, 3.4), constrained_layout=True)
    if symmetric:
        m = float(np.max(np.abs(field2d)))
        vmin, vmax = -m, m
    im = ax.imshow(field2d, origin="lower", aspect="auto", extent=ext,
                   cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_title(title)
    cb = fig.colorbar(im, ax=ax, pad=0.012)
    cb.set_label(units)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_svd(singular_values: np.ndarray, out_path: Path, *,
             title: str = "Singular-value spectrum of detrended snapshots",
             rank_marker: int | None = None) -> None:
    configure()
    fig, ax = plt.subplots(figsize=(8.0, 4.6), constrained_layout=True)
    idx = np.arange(1, singular_values.size + 1)
    ax.semilogy(idx, singular_values, "o-", color="#1f4e79", markersize=3.5)
    if rank_marker:
        ax.axvline(rank_marker + 0.5, color="#c1121f", linestyle="--",
                   label=f"retained rank = {rank_marker}")
        ax.legend()
    ax.set_xlabel("Singular-value index")
    ax.set_ylabel("Singular value")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.2)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_spectrum(freq_khz: np.ndarray, normalized_amp: np.ndarray, out_path: Path, *,
                  title: str, highlight_n: int = 0, mark_khz: Sequence[float] = ()) -> None:
    configure()
    fig, ax = plt.subplots(figsize=(9.6, 4.6), constrained_layout=True)
    bars = ax.bar(freq_khz, normalized_amp, width=max(6.0, 0.01 * float(np.max(freq_khz))),
                  color="#6c757d", edgecolor="black", linewidth=0.7)
    for k in range(min(highlight_n, len(bars))):
        bars[k].set_facecolor("#c1121f")
        bars[k].set_hatch("//")
    for fk in mark_khz:
        ax.axvline(fk, color="black", linestyle=":", linewidth=1.3, alpha=0.7)
    ax.set_xlabel("Frequency [kHz]")
    ax.set_ylabel("Normalized DMD amplitude")
    ax.set_title(title)
    ax.set_xlim(left=0.0)
    ax.set_ylim(0.0, 1.08)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_psd(curves: List[dict], out_path: Path, *, units_symbol: str, units: str,
             title: str, mark_khz: Sequence[float] = (), xlim_khz=None) -> None:
    """
    Overlay several PSD curves on one semilog-y axis.
    Each curve is a dict with keys: f_khz, psd, label, and optionally color,
    linestyle, and linewidth.
    """
    configure()
    fig, ax = plt.subplots(figsize=(10.0, 5.0), constrained_layout=True)
    eps = 1e-300
    for c in curves:
        ax.semilogy(c["f_khz"], np.maximum(c["psd"], eps), label=c["label"],
                    color=c.get("color"), linestyle=c.get("linestyle", "-"),
                    linewidth=c.get("linewidth", 1.9))
    for fk in mark_khz:
        ax.axvline(fk, color="black", linestyle=":", linewidth=1.3, alpha=0.6)
    if xlim_khz:
        ax.set_xlim(*xlim_khz)
    ax.set_xlabel("Frequency [kHz]")
    ax.set_ylabel(rf"PSD of ${units_symbol}'$  [{units}$^2$/Hz]")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.2)
    ax.legend()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_mode(real_map: np.ndarray, amp_map: np.ndarray, x_mm, y_mm, *,
              rank: int, frequency_khz: float, wavelength_mm: float,
              phase_speed_mps: float, units: str, variable_name: str,
              out_path: Path, display_y_max_mm: float | None = None) -> None:
    """Two-panel mode figure: phase-aligned real part (top) and amplitude (bottom)."""
    configure()
    ext = (float(x_mm[0]), float(x_mm[-1]), float(y_mm[0]), float(y_mm[-1]))
    fig, axes = plt.subplots(2, 1, figsize=(12, 5.0), sharex=True, constrained_layout=True)

    vmax = float(np.max(np.abs(real_map)))
    im0 = axes[0].imshow(real_map, origin="lower", aspect="auto", extent=ext,
                         cmap="RdBu_r", vmin=-vmax, vmax=vmax, interpolation="nearest")
    axes[0].set_ylabel("y [mm]")
    axes[0].set_title(f"Mode {rank:02d}: {variable_name}  (real part)")
    txt = f"f = {frequency_khz:.1f} kHz"
    if np.isfinite(phase_speed_mps):
        txt += f"\nc = {phase_speed_mps:.0f} m/s"
    axes[0].text(0.01, 0.95, txt, transform=axes[0].transAxes, va="top", ha="left",
                 fontsize=14, bbox={"facecolor": "white", "alpha": 0.7, "edgecolor": "none"})
    cb0 = fig.colorbar(im0, ax=axes[0], pad=0.01)
    cb0.set_label(units)

    im1 = axes[1].imshow(amp_map, origin="lower", aspect="auto", extent=ext,
                         cmap="turbo", vmin=0.0, vmax=float(np.max(amp_map)),
                         interpolation="nearest")
    axes[1].set_xlabel("x [mm]")
    axes[1].set_ylabel("y [mm]")
    if np.isfinite(wavelength_mm):
        axes[1].text(0.01, 0.95, rf"$\lambda$ = {wavelength_mm:.2f} mm",
                     transform=axes[1].transAxes, va="top", ha="left", color="white",
                     fontsize=14, bbox={"facecolor": "black", "alpha": 0.35, "edgecolor": "none"})
    axes[1].set_title(f"Mode {rank:02d}: {variable_name}  (amplitude)")
    cb1 = fig.colorbar(im1, ax=axes[1], pad=0.01)
    cb1.set_label(units)

    if display_y_max_mm:
        for ax in axes:
            ax.set_ylim(float(y_mm[0]), display_y_max_mm)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Animations (GIF via Pillow when ffmpeg is not available)
# --------------------------------------------------------------------------- #
def _frame_limits(cube: np.ndarray, symmetric: bool, robust: float = 0.995):
    if symmetric:
        m = float(np.quantile(np.abs(cube), robust))
        return -m, m
    return float(np.quantile(cube, 1 - robust)), float(np.quantile(cube, robust))


def animate_field(cube: np.ndarray, x_mm, y_mm, out_path: Path, *,
                  title: str, units: str, cmap: str = "inferno",
                  symmetric: bool = False, fps: int = 15, step: int = 1,
                  dpi: int = 95, times_s: np.ndarray | None = None,
                  display_y_max_mm: float | None = None):
    """Animate a snapshot cube (n_time, n_y, n_x). Returns the written path."""
    configure()
    frames = range(0, cube.shape[0], step)
    vmin, vmax = _frame_limits(cube, symmetric)
    ext = (float(x_mm[0]), float(x_mm[-1]), float(y_mm[0]), float(y_mm[-1]))
    fig, ax = plt.subplots(figsize=(11.0, 3.1), constrained_layout=True)
    im = ax.imshow(cube[0], origin="lower", aspect="auto", extent=ext, cmap=cmap,
                   vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    cb = fig.colorbar(im, ax=ax, pad=0.012)
    cb.set_label(units)
    ttl = ax.set_title(title)
    if display_y_max_mm:
        ax.set_ylim(float(y_mm[0]), display_y_max_mm)

    def update(k):
        im.set_data(cube[k])
        if times_s is not None:
            ttl.set_text(f"{title}   t = {times_s[k] * 1e6:.2f} us")
        return im, ttl

    anim = FuncAnimation(fig, update, frames=frames, blit=False)
    path = _save_anim(anim, out_path, fps, dpi)
    plt.close(fig)
    return path


def animate_pair(cube_left: np.ndarray, cube_right: np.ndarray, x_mm, y_mm, out_path: Path, *,
                 title_left: str, title_right: str, suptitle: str, units: str,
                 cmap: str = "RdBu_r", symmetric: bool = True, fps: int = 15, step: int = 1,
                 dpi: int = 95, times_s: np.ndarray | None = None,
                 display_y_max_mm: float | None = None):
    """Stacked movie, for example the original fluctuation against its DMD reconstruction."""
    configure()
    frames = range(0, min(cube_left.shape[0], cube_right.shape[0]), step)
    vmin, vmax = _frame_limits(np.concatenate([cube_left, cube_right]), symmetric)
    ext = (float(x_mm[0]), float(x_mm[-1]), float(y_mm[0]), float(y_mm[-1]))
    fig, axes = plt.subplots(2, 1, figsize=(11.0, 5.2), sharex=True, constrained_layout=True)
    ims = []
    for ax, cube, t in zip(axes, (cube_left, cube_right), (title_left, title_right)):
        im = ax.imshow(cube[0], origin="lower", aspect="auto", extent=ext, cmap=cmap,
                       vmin=vmin, vmax=vmax, interpolation="nearest")
        ax.set_ylabel("y [mm]")
        ax.set_title(t)
        fig.colorbar(im, ax=ax, pad=0.01).set_label(units)
        if display_y_max_mm:
            ax.set_ylim(float(y_mm[0]), display_y_max_mm)
        ims.append(im)
    axes[-1].set_xlabel("x [mm]")
    st = fig.suptitle(suptitle)

    def update(k):
        ims[0].set_data(cube_left[k])
        ims[1].set_data(cube_right[k])
        if times_s is not None:
            st.set_text(f"{suptitle}   t = {times_s[k] * 1e6:.2f} us")
        return (*ims, st)

    anim = FuncAnimation(fig, update, frames=frames, blit=False)
    path = _save_anim(anim, out_path, fps, dpi)
    plt.close(fig)
    return path
