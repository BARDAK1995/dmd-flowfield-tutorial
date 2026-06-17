"""
io_utils.py, structured and self-describing outputs.

Every analysis writes the same tidy bundle, so a collaborator always knows where
to look:

    <output>/figures/   PNG plots
    <output>/tables/    CSV summaries (one row per mode, one row per singular value)
    <output>/modes/     mode_XX.npz   (the raw 2-D mode fields, for re-plotting)
    <output>/analysis_metadata.json    every knob used, for exact reproduction
    <output>/README.md                 a human-readable digest of the run
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np


def ensure_dirs(output_dir: Path) -> Dict[str, Path]:
    out = Path(output_dir)
    dirs = {"root": out, "figures": out / "figures",
            "tables": out / "tables", "modes": out / "modes"}
    for p in dirs.values():
        p.mkdir(parents=True, exist_ok=True)
    return dirs


def write_csv(path: Path, rows: Sequence[Dict[str, object]],
              fieldnames: List[str] | None = None) -> None:
    rows = list(rows)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with Path(path).open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def write_svd_csv(path: Path, singular_values: np.ndarray) -> None:
    with Path(path).open("w", newline="", encoding="utf-8") as fh:
        fh.write("index,singular_value\n")
        for i, v in enumerate(singular_values, start=1):
            fh.write(f"{i},{v:.16e}\n")


def save_mode_npz(path: Path, *, x_mm, y_mm, frequency_hz, growth_rate_per_s,
                  units: str, real_map: np.ndarray, amplitude_map: np.ndarray) -> None:
    np.savez_compressed(
        path, x_mm=np.asarray(x_mm), y_mm=np.asarray(y_mm),
        frequency_hz=float(frequency_hz), growth_rate_per_s=float(growth_rate_per_s),
        units=units, real_map_units=real_map, amplitude_map_units=amplitude_map,
    )


def write_json(path: Path, payload: Dict[str, object]) -> None:
    Path(path).write_text(json.dumps(payload, indent=2, default=_json_default),
                          encoding="utf-8")


def _json_default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return str(obj)


def write_readme(path: Path, title: str, lines: Sequence[str]) -> None:
    Path(path).write_text("\n".join([f"# {title}", "", *lines]) + "\n", encoding="utf-8")
