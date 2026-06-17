"""
dmd.py: Dynamic Mode Decomposition, reported in physical units.

DMD factorises a sequence of snapshots into a set of spatial **modes**, each one
evolving in time as a single complex exponential.  Mode k has:

    * a complex eigenvalue        lambda_k        (one step of the map)
    * a frequency                 f_k = Im(log lambda_k) / (2*pi*dt)   [Hz]
    * a temporal growth rate      g_k = Re(log lambda_k) / dt          [1/s]
    * a spatial structure         phi_k                                 (a field)
    * an amplitude                b_k                                   (how much
                                                                         it weighs)

So a DMD is "an FFT in time that also hands you the 2-D spatial picture of each
frequency".  Everything below is *dimensional*: feed it the real snapshot
timestep dt and you get frequencies in Hz/kHz and phase speeds in m/s.

Pipeline (matching the reference workflow, but generalised):

    detrend snapshots  ->  SVD (singular-value spectrum)  ->  fit DMD
        ->  physical frequencies / growth rates
        ->  rank the positive-frequency modes by amplitude
        ->  per-mode real-part map, amplitude map, wavelength, phase speed
        ->  reconstruct the field from a chosen subset of modes

The reconstruction step keeps **conjugate pairs** together so the rebuilt field
stays real.  Reconstructing from "the first N modes" really means the first N
positive-frequency modes *and* their negative-frequency twins.
"""
from __future__ import annotations

from dataclasses import dataclass, field as _dc_field
from typing import Dict, List, Sequence

import numpy as np
import scipy.linalg
from pydmd import DMD


# --------------------------------------------------------------------------- #
# Detrending
# --------------------------------------------------------------------------- #
def moving_average(cube: np.ndarray, window: int) -> np.ndarray:
    """Centred moving average along time (cumsum implementation)."""
    cube = np.asarray(cube, dtype=np.float64)
    if window <= 1:
        return cube.copy()
    if window >= cube.shape[0]:
        return np.repeat(cube.mean(axis=0, keepdims=True), cube.shape[0], axis=0)
    csum = np.cumsum(cube, axis=0, dtype=np.float64)
    ma = (csum[window:] - csum[:-window]) / window
    pad0 = np.repeat(ma[:1], window // 2, axis=0)
    pad1 = np.repeat(ma[-1:], window // 2 + window % 2, axis=0)
    return np.concatenate([pad0, ma, pad1], axis=0)


def detrend_snapshots(cube: np.ndarray, mode: str = "mean", window: int = 300):
    """
    Remove a slowly varying background before DMD.

    mode = "mean"           : subtract the time-mean field (recommended for a
                              statistically steady window, the tutorial default).
    mode = "moving_average" : subtract a centred moving average (reference
                              pipeline default; use for records with residual drift).
    mode = "none"           : no detrending.

    Returns (detrended_cube, trend).  ``trend`` is (n_y, n_x) for "mean" or
    (n_time, n_y, n_x) for "moving_average"; add it back to reconstruct full fields.
    """
    cube = np.asarray(cube, dtype=np.float64)
    if mode == "mean":
        trend = cube.mean(axis=0)
        return cube - trend[None], trend
    if mode == "moving_average":
        trend = moving_average(cube, window)
        return cube - trend, trend
    if mode == "none":
        return cube, np.zeros(cube.shape[1:], dtype=np.float64)
    raise ValueError("mode must be 'mean', 'moving_average' or 'none'")


def singular_values(cube: np.ndarray) -> np.ndarray:
    """Singular values of the (space x time) snapshot matrix."""
    cube = np.asarray(cube, dtype=np.float64)
    snapshot_matrix = cube.reshape(cube.shape[0], -1).T   # (n_space, n_time)
    return scipy.linalg.svdvals(snapshot_matrix)


# --------------------------------------------------------------------------- #
# Result container
# --------------------------------------------------------------------------- #
@dataclass
class DmdResult:
    dmd: DMD
    eigs: np.ndarray            # (r,)
    modes: np.ndarray           # (n_space, r)
    amplitudes: np.ndarray      # (r,)
    dynamics: np.ndarray        # (r, n_time)  == b_k * lambda_k^t
    freq_hz: np.ndarray         # (r,)
    growth_per_s: np.ndarray    # (r,)
    pos_order: np.ndarray       # indices into the above, positive-f, |amp|-sorted
    singular_values: np.ndarray
    dt_seconds: float
    x_mm: np.ndarray
    y_mm: np.ndarray
    units: str
    symbol: str
    trend: np.ndarray = _dc_field(repr=False, default=None)   # type: ignore[assignment]
    n_time: int = 0

    @property
    def shape_yx(self):
        return self.y_mm.size, self.x_mm.size

    @property
    def n_positive(self) -> int:
        return int(self.pos_order.size)

    # per-(positive)-mode accessors
    def full_index(self, rank: int) -> int:
        """Full DMD index of the rank-th positive-frequency mode (1-based)."""
        return int(self.pos_order[rank - 1])

    def mode_info(self, rank: int) -> Dict[str, float]:
        i = self.full_index(rank)
        amp_norm = float(np.abs(self.amplitudes[i]) /
                         np.max(np.abs(self.amplitudes[self.pos_order])))
        coeff = _coeff_scale(self.dynamics[i])
        real_map = self.mode_real_map(rank)
        amp_map = self.mode_amplitude_map(rank)
        lam_mm = estimate_wavelength_mm(real_map, self.x_mm)
        f_hz = float(self.freq_hz[i])
        return {
            "rank": rank, "full_index": i,
            "frequency_hz": f_hz, "frequency_khz": f_hz / 1e3,
            "growth_rate_per_s": float(self.growth_per_s[i]),
            "amplitude_abs": float(np.abs(self.amplitudes[i])),
            "normalized_amplitude": amp_norm,
            "coefficient_rms": coeff,
            "wavelength_mm": lam_mm,
            "phase_speed_mps": phase_speed_mps(f_hz, lam_mm),
            "eigenvalue_real": float(self.eigs[i].real),
            "eigenvalue_imag": float(self.eigs[i].imag),
            "mode_real_abs_max": float(np.max(np.abs(real_map))),
            "mode_amplitude_max": float(np.max(amp_map)),
        }

    def mode_real_map(self, rank: int) -> np.ndarray:
        """Phase-aligned real part of the mode, scaled to physical units."""
        i = self.full_index(rank)
        v = _phase_align(self.modes[:, i])
        coeff = _coeff_scale(self.dynamics[i])
        return v.reshape(self.shape_yx).real * coeff

    def mode_amplitude_map(self, rank: int) -> np.ndarray:
        """Magnitude (envelope) of the mode, scaled to physical units."""
        i = self.full_index(rank)
        coeff = _coeff_scale(self.dynamics[i])
        return np.abs(self.modes[:, i].reshape(self.shape_yx)) * coeff

    # reconstruction
    def conjugate_partner(self, full_index: int) -> int:
        """Index of the mode whose eigenvalue is the complex conjugate."""
        target = np.conj(self.eigs[full_index])
        d = np.abs(self.eigs - target)
        d[full_index] = np.inf
        j = int(np.argmin(d))
        # if there is no genuine partner (a real eigenvalue), the mode is its own conjugate
        if np.abs(self.eigs[full_index].imag) < 1e-12:
            return full_index
        return j

    def selected_full_indices(self, ranks: Sequence[int]) -> List[int]:
        """Expand positive-frequency ranks into a conjugate-closed index set."""
        sel = set()
        for r in ranks:
            i = self.full_index(r)
            sel.add(i)
            sel.add(self.conjugate_partner(i))
        return sorted(sel)

    def reconstruct(self, ranks: Sequence[int], *, add_mean: bool = True,
                    time_indices: np.ndarray | None = None) -> np.ndarray:
        """
        Rebuild the (n_time, n_y, n_x) field from the chosen positive-frequency
        modes (their conjugate twins are included automatically so the result is
        real).

        ranks       : 1-based positive-frequency mode ranks, e.g. [1, 2, 3].
        add_mean    : if True, add the detrended-out background back for the full
                      field; if False, return only the reconstructed fluctuation.
        time_indices: subset of snapshots to rebuild (default: all).
        """
        ti = np.arange(self.n_time) if time_indices is None else np.asarray(time_indices)
        sel = self.selected_full_indices(ranks)
        recon = (self.modes[:, sel] @ self.dynamics[np.ix_(sel, ti)]).real  # (n_space, nt)
        ny, nx = self.shape_yx
        cube = recon.T.reshape(ti.size, ny, nx)
        if add_mean and self.trend is not None:
            if self.trend.ndim == 2:
                cube = cube + self.trend[None]
            else:
                cube = cube + self.trend[ti]
        return cube


# --------------------------------------------------------------------------- #
# Internal helpers
# --------------------------------------------------------------------------- #
def _phase_align(mode_vec: np.ndarray) -> np.ndarray:
    ref = int(np.argmax(np.abs(mode_vec)))
    return mode_vec * np.exp(-1j * np.angle(mode_vec[ref]))


def _coeff_scale(dynamics_row: np.ndarray) -> float:
    """RMS of a mode's time coefficient, which gives the mode a physical amplitude."""
    return float(np.sqrt(np.mean(np.abs(dynamics_row) ** 2)))


def eigs_to_frequency_hz(eigs: np.ndarray, dt_seconds: float) -> np.ndarray:
    return np.imag(np.log(eigs)) / (2.0 * np.pi * dt_seconds)


def eigs_to_growth_per_s(eigs: np.ndarray, dt_seconds: float) -> np.ndarray:
    return np.real(np.log(eigs)) / dt_seconds


# --------------------------------------------------------------------------- #
# The fit
# --------------------------------------------------------------------------- #
def run_dmd(cube: np.ndarray, dt_seconds: float, x_mm: np.ndarray, y_mm: np.ndarray, *,
            units: str = "", symbol: str = "",
            detrend: str = "mean", moving_average_window: int = 300,
            svd_rank: int = 30, exact: bool = True, opt: bool = True,
            forward_backward: bool = True, sorted_eigs: str = "abs",
            freq_min_hz: float = 1.0) -> DmdResult:
    """
    Detrend, compute the singular-value spectrum, fit a DMD, and package the
    physically-dimensioned result.

    Defaults reproduce the reference pipeline's DMD family (forward-backward
    exact DMD, optimal amplitudes, rank 30), but with ``detrend='mean'``, which is
    the natural choice for a statistically steady analysis window.
    """
    cube = np.asarray(cube, dtype=np.float64)
    n_time = cube.shape[0]

    detrended, trend = detrend_snapshots(cube, mode=detrend, window=moving_average_window)
    sv = singular_values(detrended)

    snapshots = [detrended[k] for k in range(n_time)]
    dmd = DMD(svd_rank=svd_rank, exact=exact, opt=opt,
              forward_backward=forward_backward, sorted_eigs=sorted_eigs)
    dmd.fit(snapshots)

    eigs = np.asarray(dmd.eigs)
    freq_hz = eigs_to_frequency_hz(eigs, dt_seconds)
    growth = eigs_to_growth_per_s(eigs, dt_seconds)
    amps = np.abs(np.asarray(dmd.amplitudes))

    positive = np.flatnonzero(freq_hz > freq_min_hz)
    if positive.size == 0:
        raise RuntimeError("No positive-frequency DMD modes found.")
    pos_order = positive[np.argsort(amps[positive])[::-1]]

    return DmdResult(
        dmd=dmd, eigs=eigs, modes=np.asarray(dmd.modes),
        amplitudes=np.asarray(dmd.amplitudes), dynamics=np.asarray(dmd.dynamics),
        freq_hz=freq_hz, growth_per_s=growth, pos_order=pos_order,
        singular_values=sv, dt_seconds=dt_seconds, x_mm=x_mm, y_mm=y_mm,
        units=units, symbol=symbol, trend=trend, n_time=n_time,
    )


# --------------------------------------------------------------------------- #
# Spatial wavelength / phase speed (dimensional, no stability theory)
# --------------------------------------------------------------------------- #
def estimate_wavelength_mm(real_map: np.ndarray, x_mm: np.ndarray,
                           y_index: int | None = None) -> float:
    """
    Dominant streamwise wavelength of a mode's real part, in mm.

    The mode is collapsed to a 1-D streamwise signal (averaged over y, or sampled
    at a chosen row), then a windowed real-FFT in x finds the peak wavenumber.
    No boundary-layer line or edge criterion, just the spatial period.
    """
    real_map = np.asarray(real_map, dtype=np.float64)
    signal = real_map.mean(axis=0) if y_index is None else real_map[y_index, :]
    signal = signal - signal.mean()
    if np.allclose(signal, 0.0) or signal.size < 4:
        return float("nan")
    dx_mm = float(np.mean(np.diff(x_mm)))
    spec = np.fft.rfft(signal * np.hanning(signal.size))
    k_per_mm = np.fft.rfftfreq(signal.size, d=dx_mm)
    pos = np.flatnonzero(k_per_mm > 0.0)
    if pos.size == 0:
        return float("nan")
    k_peak = k_per_mm[pos][int(np.argmax(np.abs(spec[pos])))]
    return 1.0 / k_peak if k_peak > 0 else float("nan")


def phase_speed_mps(frequency_hz: float, wavelength_mm: float) -> float:
    """c = f * lambda, returned in m/s (lambda given in mm)."""
    if not np.isfinite(wavelength_mm):
        return float("nan")
    return float(frequency_hz * wavelength_mm * 1.0e-3)


def reconstruction_error(reference_cube: np.ndarray, recon_cube: np.ndarray) -> float:
    """Relative L2 error ||ref - recon|| / ||ref|| over the whole cube."""
    ref = np.asarray(reference_cube, dtype=np.float64)
    rec = np.asarray(recon_cube, dtype=np.float64)
    denom = float(np.linalg.norm(ref))
    return float(np.linalg.norm(ref - rec) / denom) if denom > 0 else float("nan")
