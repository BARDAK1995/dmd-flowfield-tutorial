"""
dataset.py loads 2-D snapshot data and its physical metadata.

The whole tutorial is built around one idea: a *case* is a stack of 2-D field
snapshots plus a tiny text file that tells you the timestep and the domain.

    field arrays : <name>_<unit>.npy, shape (n_time, n_y, n_x)   memmap-friendly
    fieldinputs  : fieldinputs.dat   dt (s), x range, y range, unit

The toolkit reads the resolution and snapshot count from the array shapes and the
field name and units from the filename, so fieldinputs.dat only carries the
timestep, the x and y ranges, and the unit. From those, every physical quantity
(time in seconds, frequency in kHz, a probe location in mm, a wavelength in mm)
follows mechanically. There is no boundary-layer edge, no Reynolds number, and no
stability theory, just plain dimensional bookkeeping. A richer
dataset_metadata.json sidecar is also supported when present.
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

    # grid (uniform cell centres, in millimetres)
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

    # ==== convenience accessors ========================================== #
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

    # ==== shapes / axes ================================================== #
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

    # ==== point / window selection ====================================== #
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


# Map a filename unit token and a field key to display units and a math symbol.
_UNIT_TOKENS = {"m3": "m^-3", "Pa": "Pa", "K": "K", "mps": "m/s", "ms": "m/s"}
_SYMBOLS = {"number_density": "n", "density": "rho", "pressure": "p",
            "temperature": "T", "u_velocity": "u", "v_velocity": "v",
            "u": "u", "v": "v", "temp": "T"}
_LENGTH_TO_MM = {"mm": 1.0, "cm": 10.0, "m": 1000.0, "um": 1.0e-3, "µm": 1.0e-3,
                 "micron": 1.0e-3}


def _field_from_filename(stem: str) -> Dict[str, str]:
    """
    Turn an array filename into a field spec. The convention is
    ``<name>_<unit>.npy``, for example ``pressure_Pa.npy`` or
    ``number_density_m3.npy``. The part after the last underscore is read as the
    unit if it is a known token; otherwise the whole name is the key.
    """
    parts = stem.split("_")
    if len(parts) >= 2 and parts[-1] in _UNIT_TOKENS:
        key = "_".join(parts[:-1])
        units = _UNIT_TOKENS[parts[-1]]
    else:
        key, units = stem, ""
    name = key.replace("_", " ").strip()
    name = name[:1].upper() + name[1:] if name else key
    symbol = _SYMBOLS.get(key, key[:1] if key else "q")
    return {"file": f"{stem}.npy", "name": name, "symbol": symbol, "units": units}


def _read_data_lines(path: Path):
    """Non-comment, non-blank lines of a .dat file (comments start with #)."""
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.split("#", 1)[0].strip()
        if s:
            out.append(s)
    return out


def metadata_from_fieldinputs(data_dir: Path) -> Dict:
    """
    Build the metadata dictionary from a simple ``fieldinputs.dat`` plus the
    .npy files in the folder. The .dat supplies the three physical inputs that
    the arrays cannot carry: the timestep, the x and y domain ranges, and the
    length unit. The grid resolution and snapshot count come from the array shapes.

    Four values, one per line (comments start with #):
        1. dt in seconds
        2. the x range as "x_min x_max" (for example 40 100)
        3. the y range as "y_min y_max" (for example 0 2.5)
        4. the length unit (for example mm)
    """
    lines = _read_data_lines(data_dir / "fieldinputs.dat")
    if len(lines) < 4:
        raise ValueError("fieldinputs.dat needs four values: dt, 'x_min x_max', "
                         "'y_min y_max', and the unit")
    dt = float(lines[0])
    x_lo, x_hi = (float(v) for v in lines[1].split()[:2])
    y_lo, y_hi = (float(v) for v in lines[2].split()[:2])
    unit = lines[3].split()[0]
    to_mm = _LENGTH_TO_MM.get(unit.lower(), 1.0)

    field_files = sorted(p for p in data_dir.glob("*.npy"))
    if not field_files:
        raise ValueError(f"no .npy field arrays found in {data_dir}")
    sample = np.load(field_files[0], mmap_mode="r")
    nt, ny, nx = (int(v) for v in sample.shape)
    fields = {}
    for p in field_files:
        spec = _field_from_filename(p.stem)
        # the field key (what --field expects) is the filename without the unit token
        key = p.stem.rsplit("_", 1)[0] if p.stem.split("_")[-1] in _UNIT_TOKENS else p.stem
        fields[key] = spec

    return {
        "title": f"Dataset described by fieldinputs.dat ({data_dir.name})",
        "grid": {"x_min_mm": x_lo * to_mm, "x_max_mm": x_hi * to_mm,
                 "y_min_mm": y_lo * to_mm, "y_max_mm": y_hi * to_mm,
                 "nx": nx, "ny": ny, "note": "uniform cell centres; domain from fieldinputs.dat"},
        "time": {"n_snapshots": nt, "step_first": 0, "step_interval": 1,
                 "dt_seconds": dt, "sampling_frequency_hz": 1.0 / dt, "nyquist_hz": 0.5 / dt,
                 "domain_unit": unit},
        "fields": fields,
    }


def load_metadata(data_dir: str | Path) -> CaseMetadata:
    """
    Load the physical description of a dataset. Two inputs are supported:

    * ``fieldinputs.dat`` (preferred for simple cases): four commented lines
      giving the timestep, the x range, the y range, and the length unit. The
      fields and their units are read from the .npy filenames.
    * ``dataset_metadata.json``: the richer sidecar, used when present and there
      is no fieldinputs.dat. It can carry extra context like forcing and probes.
    """
    data_dir = Path(data_dir)
    if (data_dir / "fieldinputs.dat").exists():
        raw = metadata_from_fieldinputs(data_dir)
    else:
        raw = json.loads((data_dir / "dataset_metadata.json").read_text(encoding="utf-8"))
    return CaseMetadata(raw=raw, data_dir=data_dir)


def load_field(data_dir: str | Path, key: str, mmap: bool = True) -> Field:
    """
    Load one field cube by its metadata key (e.g. ``"number_density"`` or
    ``"pressure"``). With ``mmap=True`` it stays memory-mapped, so large arrays
    load instantly. Pass ``mmap=False`` to pull it fully into RAM.
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
    that produced it; see tools/prepare_tutorial_data.py.)
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
