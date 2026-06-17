"""
fields.py: time-mean, RMS fields, and fluctuation cubes.

These answer "what does the flow look like, and where does it fluctuate."
Two conventions matter, so let's be explicit about them:

* The *time-mean* field is a simple average over the time axis.

* The *RMS field* is the root-mean-square of the fluctuation **after removing a
  least-squares straight line in time** at every grid point.  Removing a linear
  trend (rather than just the mean) stops a slow residual drift, such as a
  transient that has not fully converged yet, from masquerading as
  turbulence or instability energy.  This matches the reference pipeline's RMS
  definition.

All functions take an (n_time, n_y, n_x) cube and return (n_y, n_x) maps, so they
compose directly with :class:`dmdkit.dataset.Field`.
"""
from __future__ import annotations

import numpy as np


# --------------------------------------------------------------------------- #
# Detrending primitives
# --------------------------------------------------------------------------- #
def _centered_time(n: int) -> np.ndarray:
    """Symmetric time index 0..n-1 recentred about its mean (for stable fits)."""
    return np.arange(n, dtype=np.float64) - 0.5 * (n - 1)


def linear_detrend(cube: np.ndarray) -> np.ndarray:
    """
    Remove the least-squares line along axis 0 (time) at every (y, x) point.

    Vectorised closed form: slope = <t, x> / <t, t>, then trend = mean + t * slope.
    """
    cube = np.asarray(cube, dtype=np.float64)
    n = cube.shape[0]
    if n <= 1:
        return cube - cube.mean(axis=0, keepdims=True)
    t = _centered_time(n)
    denom = float(np.dot(t, t))
    mean = cube.mean(axis=0, keepdims=True)
    slope = np.tensordot(t, cube, axes=(0, 0)) / denom
    trend = mean + t.reshape((n,) + (1,) * (cube.ndim - 1)) * slope
    return cube - trend


# --------------------------------------------------------------------------- #
# Public diagnostics
# --------------------------------------------------------------------------- #
def mean_field(cube: np.ndarray) -> np.ndarray:
    """Time-mean field, shape (n_y, n_x)."""
    return np.asarray(cube, dtype=np.float64).mean(axis=0)


def rms_field(cube: np.ndarray, detrend: str = "linear") -> np.ndarray:
    """
    RMS fluctuation field, shape (n_y, n_x).

    detrend = "linear" : subtract a per-point least-squares line in time (default,
                         matches the reference pipeline).
    detrend = "mean"   : subtract only the per-point time mean.
    """
    cube = np.asarray(cube, dtype=np.float64)
    if detrend == "linear":
        resid = linear_detrend(cube)
    elif detrend == "mean":
        resid = cube - cube.mean(axis=0, keepdims=True)
    else:
        raise ValueError("detrend must be 'linear' or 'mean'")
    return np.sqrt(np.mean(resid ** 2, axis=0))


def fluctuation_cube(cube: np.ndarray, detrend: str = "mean") -> np.ndarray:
    """
    Fluctuation field at every snapshot (cube minus its trend), useful for
    animating the unsteady part of the flow without the large steady background.
    """
    cube = np.asarray(cube, dtype=np.float64)
    if detrend == "linear":
        return linear_detrend(cube)
    if detrend == "mean":
        return cube - cube.mean(axis=0, keepdims=True)
    if detrend == "none":
        return cube
    raise ValueError("detrend must be 'linear', 'mean' or 'none'")


def field_stats(field2d: np.ndarray) -> dict:
    """Summary statistics for a 2-D map (used in tables / annotations)."""
    f = np.asarray(field2d, dtype=np.float64)
    return {
        "min": float(f.min()), "mean": float(f.mean()), "max": float(f.max()),
        "p95": float(np.quantile(f, 0.95)), "p99": float(np.quantile(f, 0.99)),
    }
