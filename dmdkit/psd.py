"""
psd.py -- point power-spectral density (PSD).

Given a time series sampled at fixed dt, where is its energy in frequency?  This
is how you tell a *tone* (a sharp peak -- e.g. a forced or resonant instability)
from *broadband noise* (a flat floor -- e.g. DSMC statistical scatter in the
free stream).

Conventions (matching the reference pipeline):

* sampling frequency  fs = 1 / dt
* one-sided PSD with **density** scaling  -> units are [signal]^2 / Hz
* a Hann window and a **linear** detrend are applied inside the estimator
* two estimators are offered:
      - "periodogram" : full-resolution single-segment estimate (df = fs / N)
      - "welch"       : segment-averaged estimate (smoother, lower variance)

The integral of the one-sided PSD over frequency returns the variance, so
``rms_from_psd`` recovers the time-domain RMS:  sqrt( integral PSD df ).

Frequencies are returned in Hz; helpers convert to kHz for plotting/reporting.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy.signal import periodogram, welch


def point_psd(series: np.ndarray, dt_seconds: float, *,
              method: str = "periodogram", window: str = "hann",
              detrend: str = "linear", nperseg: int | None = None
              ) -> Tuple[np.ndarray, np.ndarray]:
    """
    One-sided PSD of a 1-D time series.

    Returns (frequency_hz, psd).  ``psd`` has units [series]^2 / Hz.

    method = "periodogram" : single segment, finest frequency resolution.
    method = "welch"       : averaged segments (set ``nperseg``, default N//8),
                             trades resolution for a smoother, lower-variance PSD.
    """
    series = np.asarray(series, dtype=np.float64)
    fs = 1.0 / float(dt_seconds)
    if method == "periodogram":
        f_hz, pxx = periodogram(series, fs=fs, window=window, detrend=detrend,
                                scaling="density", return_onesided=True)
    elif method == "welch":
        if nperseg is None:
            nperseg = max(8, series.size // 8)
        f_hz, pxx = welch(series, fs=fs, window=window, detrend=detrend,
                          nperseg=nperseg, scaling="density", return_onesided=True)
    else:
        raise ValueError("method must be 'periodogram' or 'welch'")
    return f_hz, pxx


def rms_from_psd(freq_hz: np.ndarray, psd: np.ndarray) -> float:
    """
    Recover the time-domain RMS from a one-sided density PSD:

        variance = integral_0^{fNyq} PSD(f) df,   RMS = sqrt(variance).
    """
    return float(np.sqrt(np.trapz(np.asarray(psd), np.asarray(freq_hz))))


def band_rms_from_psd(freq_hz: np.ndarray, psd: np.ndarray,
                      f_lo_hz: float, f_hi_hz: float) -> float:
    """RMS contained in a frequency band [f_lo, f_hi] (Hz)."""
    f = np.asarray(freq_hz)
    mask = (f >= f_lo_hz) & (f <= f_hi_hz)
    if mask.sum() < 2:
        return 0.0
    return float(np.sqrt(np.trapz(np.asarray(psd)[mask], f[mask])))


def peak_frequency_hz(freq_hz: np.ndarray, psd: np.ndarray,
                      f_min_hz: float = 0.0) -> float:
    """Frequency of the largest PSD value above ``f_min_hz`` (Hz)."""
    f = np.asarray(freq_hz)
    p = np.asarray(psd)
    mask = f > f_min_hz
    if not np.any(mask):
        return float("nan")
    idx = np.flatnonzero(mask)[int(np.argmax(p[mask]))]
    return float(f[idx])


def to_khz(freq_hz: np.ndarray) -> np.ndarray:
    return np.asarray(freq_hz) / 1.0e3
