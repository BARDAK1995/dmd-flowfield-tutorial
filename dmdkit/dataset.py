"""
dataset.py -- load 2-D snapshot data and its physical metadata.

The whole tutorial is built around one idea: a *case* is a stack of 2-D field
snapshots plus a small JSON file that tells you the grid and the timestep.

    field array  : shape (n_time, n_y, n_x)         a .npy memmap-friendly cube
    metadata     : dataset_metadata.json            grid (mm), timestep (s), units

Given those two things, every physical quantity (time in seconds, frequency in
kHz, a probe location in mm, a wavelength in mm) follows mechanically.  No
boundary-layer edge, no Reynolds number, no stability theory -- purely
dimensional bookkeeping.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field as _dc_field
from pathlib import Path
from typing import Dict, Tuple

import numpy as np


# --------------------------------------------------------------------------- #
# Metadata container
# --------------------------------------------------------------------------- #
@dataclass
class CaseMetadata:
    """Parsed view of ``dataset_metadata.json``."""

    raw: Dict
    data_dir: Path

    # grid (uniform cell centres, millimetres)
    x_mm: np.ndarray = _dc_field(repr=False, default=None)  # type: ignore[assignment]
    y_mm: np.ndarray = _dc_field(repr=False, default=None)  # type: ignore[assignment]

    # time axis
    dt_seconds: float = 0.0
    step_first: int = 0
    step_interval: int = 0

    def __post_init__(self) -> None:
        grid = self.raw["grid"]
        self.x_mm = uniform_cell_centers(grid["x_min_mm"], grid["x_max_mm"], int(grid["nx"]))
        self.y_mm = uniform_cell_centers(grid["y_min_mm"], grid["y_max_mm"], int(grid["ny"]))
        t = self.raw["time"]
        self.dt_seconds = float(t["dt_seconds"])
        self.step_first = int(t["step_first"])
        self.step_interval = int(t["step_interval"])

    # -- convenience accessors --------------------------------------------- #
    @property
    def fs_hz(self) -> float:
        """Snapshot sampling frequency (Hz)."""
        return 1.0 / self.dt_seconds

    @property
    def nyquist_hz(self) -> float:
        return 0.5 * self.fs_hz

    def field_spec(self, key: str) -> Dict[str, str]:
        return self.raw["fields"][key]

    def field_keys(self):
        return list(self.raw["fields"].keys())


# --------------------------------------------------------------------------- #
# Field container
# --------------------------------------------------------------------------- #
@dataclass
class Field:
    """A loaded field cube with its coordinates and labels."""

    data: np.ndarray            # (n_time, n_y, n_x)
    x_mm: np.ndarray            # (n_x,)
    y_mm: np.ndarray            # (n_y,)
    dt_seconds: float
    step_first: int
    step_interval: int
    name: str                   # e.g. "Pressure"
    symbol: str                 # e.g. "p"
    units: str                  # e.g. "Pa"
    meta: CaseMetadata

    # -- shapes / axes ----------------------------------------------------- #
    @property
    def n_time(self) -> int:
        return self.data.shape[0]

    @property
    def shape_yx(self) -> Tuple[int, int]:
        return self.data.shape[1], self.data.shape[2]

    @property
    def fs_hz(self) -> float:
        return 1.0 / self.dt_seconds

    def times_s(self) -> np.ndarray:
        """Physical time of each snapshot, seconds (t=0 at first snapshot)."""
        return np.arange(self.n_time, dtype=np.float64) * self.dt_seconds

    def steps(self) -> np.ndarray:
        """Solver timestep number of each snapshot."""
        return self.step_first + np.arange(self.n_time) * self.step_interval

    def extent_mm(self):
        """matplotlib imshow extent: (x0, x1, y0, y1) in mm."""
        return (float(self.x_mm[0]), float(self.x_mm[-1]),
                float(self.y_mm[0]), float(self.y_mm[-1]))

    # -- point / window selection ----------------------------------------- #
    def nearest_ij(self, x_mm: float, y_mm: float) -> Tuple[int, int]:
        ix = int(np.argmin(np.abs(self.x_mm - x_mm)))
        iy = int(np.argmin(np.abs(self.y_mm - y_mm)))
        return iy, ix

    def point_series(self, x_mm: float, y_mm: float) -> Dict[str, object]:
        """Time series at the grid point nearest to (x_mm, y_mm)."""
        iy, ix = self.nearest_ij(x_mm, y_mm)
        return {
            "series": np.asarray(self.data[:, iy, ix], dtype=np.float64),
            "x_requested_mm": float(x_mm), "y_requested_mm": float(y_mm),
            "x_actual_mm": float(self.x_mm[ix]), "y_actual_mm": float(self.y_mm[iy]),
            "ix": ix, "iy": iy,
        }

    def crop(self, x_min_mm=None, x_max_mm=None, y_min_mm=None, y_max_mm=None) -> "Field":
        """Return a spatially cropped copy (used to focus the DMD region)."""
        xmask = np.ones_like(self.x_mm, dtype=bool)
        ymask = np.ones_like(self.y_mm, dtype=bool)
        if x_min_mm is not None:
            xmask &= self.x_mm >= x_min_mm
        if x_max_mm is not None:
            xmask &= self.x_mm <= x_max_mm
        if y_min_mm is not None:
            ymask &= self.y_mm >= y_min_mm
        if y_max_mm is not None:
            ymask &= self.y_mm <= y_max_mm
        xi = np.flatnonzero(xmask)
        yi = np.flatnonzero(ymask)
        if xi.size == 0 or yi.size == 0:
            raise ValueError("crop() removed the whole domain; check the mm bounds")
        sub = np.asarray(self.data[:, yi[0]:yi[-1] + 1, xi[0]:xi[-1] + 1], dtype=np.float64)
        return Field(sub, self.x_mm[xi], self.y_mm[yi], self.dt_seconds,
                     self.step_first, self.step_interval,
                     self.name, self.symbol, self.units, self.meta)

    def time_window(self, start: int = 0, stop: int | None = None) -> "Field":
        """Return a copy restricted to snapshot indices [start, stop)."""
        stop = self.n_time if stop is None else stop
        sub = np.asarray(self.data[start:stop], dtype=np.float64)
        return Field(sub, self.x_mm, self.y_mm, self.dt_seconds,
                     self.step_first + start * self.step_interval, self.step_interval,
                     self.name, self.symbol, self.units, self.meta)


# --------------------------------------------------------------------------- #
# Loading helpers
# --------------------------------------------------------------------------- #
def uniform_cell_centers(lo_mm: float, hi_mm: float, n: int) -> np.ndarray:
    """Cell-centre coordinates for ``n`` uniform cells spanning [lo, hi]."""
    edges = np.linspace(float(lo_mm), float(hi_mm), n + 1)
    return 0.5 * (edges[:-1] + edges[1:])


def load_metadata(data_dir: str | Path) -> CaseMetadata:
    data_dir = Path(data_dir)
    meta_path = data_dir / "dataset_metadata.json"
    raw = json.loads(meta_path.read_text(encoding="utf-8"))
    return CaseMetadata(raw=raw, data_dir=data_dir)


def load_field(data_dir: str | Path, key: str, mmap: bool = True) -> Field:
    """
    Load one field cube by its metadata key (e.g. ``"number_density"`` or
    ``"pressure"``).  ``mmap=True`` keeps it memory-mapped so large arrays load
    instantly; pass ``mmap=False`` to pull it fully into RAM.
    """
    meta = load_metadata(data_dir)
    spec = meta.field_spec(key)
    arr = np.load(Path(data_dir) / spec["file"], mmap_mode="r" if mmap else None)
    if not mmap:
        arr = np.asarray(arr, dtype=np.float64)
    return Field(
        data=arr, x_mm=meta.x_mm, y_mm=meta.y_mm,
        dt_seconds=meta.dt_seconds, step_first=meta.step_first,
        step_interval=meta.step_interval,
        name=spec["name"], symbol=spec["symbol"], units=spec["units"], meta=meta,
    )


# --------------------------------------------------------------------------- #
# Dimensional conversions (DSMC-specific, but optional & self-checking)
# --------------------------------------------------------------------------- #
def number_density_to_mass_density(n_m3: np.ndarray, molecular_mass_kg: float) -> np.ndarray:
    """rho [kg/m^3] = n [1/m^3] * m_molecule [kg]."""
    return n_m3 * float(molecular_mass_kg)


def counts_to_number_density(counts: np.ndarray, fnum: float, cell_volume_m3: float) -> np.ndarray:
    """
    Convert raw DSMC cell particle *counts* to number density:

        n = counts * FNUM / V_cell

    FNUM is the number of real molecules represented by one simulated molecule.
    (The tutorial ships number density directly, but this documents the chain
    that produced it -- see tools/prepare_tutorial_data.py.)
    """
    return np.asarray(counts, dtype=np.float64) * float(fnum) / float(cell_volume_m3)


def verify_freestream_number_density(field: Field, x_mm: float, y_mm: float,
                                     reference_m3: float) -> Dict[str, float]:
    """
    Sanity check that the dimensional number-density field reads the documented
    freestream value at a point that sits in the undisturbed stream.
    """
    s = field.point_series(x_mm, y_mm)
    mean = float(np.mean(s["series"]))
    return {
        "probe_x_mm": s["x_actual_mm"], "probe_y_mm": s["y_actual_mm"],
        "mean_value": mean, "reference": float(reference_m3),
        "relative_error": (mean - reference_m3) / reference_m3 if reference_m3 else float("nan"),
    }
