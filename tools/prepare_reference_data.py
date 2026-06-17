#!/usr/bin/env python3
"""
prepare_reference_data.py builds the tutorial dataset from the original
second-mode reference case.

This is the provenance script for the HypersonicFlowOverPlate example. It takes
the last N snapshots of the unforced Mach 6 flat-plate reference run (the case
the second-mode DMD analysis in the paper was built on), calibrates the fields
to physical units, and writes the two arrays plus a metadata sidecar.

What it does
------------
1. Loads the reference N (particle counts) and P (raw pressure) snapshot cubes.
   Their stored axis order is (time, x, y), so we transpose to (time, y, x).
2. Keeps the last 1000 snapshots (the converged, statistically steady tail).
3. Calibrates to physical units against the validated freestream of this case
   (number density 1.2e24 m^-3 and pressure 828 Pa, from the solver output and
   inputs). The raw counts and raw pressure are each linearly proportional to
   the physical quantity, so a single freestream-anchored factor per field is
   exact:
       number density  n = counts * (n_inf / freestream_counts)
       pressure        p = P_raw  * (p_inf / freestream_P_raw)
4. Writes number_density_m3.npy, pressure_Pa.npy (float32), and
   dataset_metadata.json into the example data directory.

The grid is x = 40 to 100 mm by 1501 points and y = 0 to 2.5 mm by 63 points,
and the snapshot spacing is dt = 1e-7 s (so fs = 10 MHz, Nyquist 5 MHz). The
case is unforced, so the dominant DMD mode is the natural second mode near
250 kHz rather than a driven tone.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TUTORIAL_DIR = HERE.parent
DEFAULT_SOURCE = TUTORIAL_DIR.parent / "2ndMode_refDSMC" / "Datas" / "expensive_Ref"

# Validated freestream of this case (solver output 2d_out.tec and inputs).
N_INF_M3 = 1.2e24       # freestream number density, m^-3
P_INF_PA = 828.0        # freestream pressure, Pa
U_INF_MPS = 860.6
T_INF_K = 50.0
MACH = 5.97
N2_MASS_KG = 4.65e-26
GAS_R = 296.89
GAMMA = 1.4
KB = 1.380649e-23

# Physical extent of the stored region (streamwise offset already applied).
X_MIN_MM, X_MAX_MM = 40.0, 100.0
Y_MIN_MM, Y_MAX_MM = 0.0, 2.5
DT_SECONDS = 1.0e-7


def freestream_factor(counts_or_raw_top, t_top, target):
    """Anchor factor so the freestream (top rows, T near T_inf) maps to target."""
    mask = (t_top > 48.0) & (t_top < 55.0)
    ref = float(counts_or_raw_top[mask].mean())
    return target / ref, ref


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, default=DEFAULT_SOURCE,
                    help="reference case folder (holds temp/N_snapshots.npy etc.)")
    ap.add_argument("--out", type=Path,
                    default=TUTORIAL_DIR / "examples" / "HypersonicFlowOverPlate" / "data")
    ap.add_argument("--n-snapshots", type=int, default=1000,
                    help="number of trailing snapshots to keep")
    args = ap.parse_args()

    n_path = args.source / "temp" / "N_snapshots.npy"
    p_path = args.source / "temp" / "P_snapshots.npy"
    t_path = args.source / "T_snapshots.npy"

    N = np.load(n_path, mmap_mode="r")   # (nt, x, y) counts
    P = np.load(p_path, mmap_mode="r")   # (nt, x, y) raw pressure
    T = np.load(t_path, mmap_mode="r")   # (nt, x, y) Kelvin
    nt = N.shape[0]
    start = nt - args.n_snapshots

    # calibration anchors from the freestream (top rows, upstream x, T ~ T_inf)
    t_top = np.asarray(T[start:, 0:400, 58:63], dtype=np.float64)
    fac_n, fs_counts = freestream_factor(np.asarray(N[start:, 0:400, 58:63], dtype=np.float64),
                                         t_top, N_INF_M3)
    fac_p, fs_praw = freestream_factor(np.asarray(P[start:, 0:400, 58:63], dtype=np.float64),
                                       t_top, P_INF_PA)

    print(f"[source] {args.source}")
    print(f"[time]   {nt} snapshots stored; keeping last {args.n_snapshots} "
          f"(index {start}:{nt}), dt={DT_SECONDS:.1e}s, fs={1/DT_SECONDS:.1e}Hz")
    print(f"[calib]  freestream counts={fs_counts:.1f} -> n factor={fac_n:.4e} "
          f"(freestream n = {fs_counts*fac_n:.3e})")
    print(f"[calib]  freestream P_raw={fs_praw:.2f} -> p factor={fac_p:.4e} "
          f"(freestream p = {fs_praw*fac_p:.1f} Pa)")

    # load tail, transpose (time, x, y) -> (time, y, x), calibrate, cast to float32
    n_out = (np.asarray(N[start:], dtype=np.float64).transpose(0, 2, 1) * fac_n).astype(np.float32)
    p_out = (np.asarray(P[start:], dtype=np.float64).transpose(0, 2, 1) * fac_p).astype(np.float32)
    ny, nx = n_out.shape[1], n_out.shape[2]

    args.out.mkdir(parents=True, exist_ok=True)
    np.save(args.out / "number_density_m3.npy", n_out)
    np.save(args.out / "pressure_Pa.npy", p_out)

    # freestream sanity check on the written arrays (top of domain, upstream x)
    iy_fs, ix_fs = ny - 3, nx // 6
    n_fs = float(n_out[:, iy_fs, ix_fs].mean())
    p_fs = float(p_out[:, iy_fs, ix_fs].mean())
    print(f"[check]  written freestream-ish n ~ {n_fs:.3e} m^-3, p ~ {p_fs:.1f} Pa")

    meta = {
        "title": "Hypersonic flat-plate DSMC: second-mode reference case",
        "description": ("Last 1000 snapshots of the unforced Mach 6 N2 flat-plate "
                        "reference run, the case the paper's second-mode DMD analysis "
                        "was built on. Fields calibrated to physical units against the "
                        "validated freestream. Two fields: number density and pressure."),
        "provenance": {
            "source_case": "2ndMode_refDSMC/Datas/expensive_Ref (unforced reference)",
            "original_n_snapshots": int(nt),
            "kept_last": int(args.n_snapshots),
            "stored_axis_order": "(time, x, y), transposed here to (time, y, x)",
            "calibration": {
                "number_density": f"n = counts * {fac_n:.6e} (anchored to n_inf = {N_INF_M3:.3e})",
                "pressure": f"p = P_raw * {fac_p:.6e} (anchored to p_inf = {P_INF_PA} Pa)",
                "freestream_counts": fs_counts, "freestream_P_raw": fs_praw,
            },
        },
        "grid": {"x_min_mm": X_MIN_MM, "x_max_mm": X_MAX_MM,
                 "y_min_mm": Y_MIN_MM, "y_max_mm": Y_MAX_MM,
                 "nx": int(nx), "ny": int(ny), "note": "uniform cell centres in mm"},
        "time": {"n_snapshots": int(n_out.shape[0]), "step_first": 0, "step_interval": 1,
                 "dt_seconds": DT_SECONDS, "sampling_frequency_hz": 1.0 / DT_SECONDS,
                 "nyquist_hz": 0.5 / DT_SECONDS,
                 "note": "step_first/step_interval are nominal snapshot indices; dt is physical"},
        "fields": {
            "number_density": {"file": "number_density_m3.npy", "name": "Number density",
                               "symbol": "n", "units": "m^-3", "dtype": "float32"},
            "pressure": {"file": "pressure_Pa.npy", "name": "Pressure",
                         "symbol": "p", "units": "Pa", "dtype": "float32"},
        },
        "forcing": {"type": "none (unforced reference; natural second-mode instability)"},
        "freestream": {"gas": "N2", "mach": MACH, "U_mps": U_INF_MPS, "T_K": T_INF_K,
                       "number_density_m3": N_INF_M3, "pressure_Pa": P_INF_PA,
                       "molecular_mass_kg": N2_MASS_KG, "gas_constant_J_kgK": GAS_R,
                       "gamma": GAMMA},
        "suggested_probes": {
            "boundary": {"x_mm": 70.0, "y_mm": 0.5,
                         "note": "inside the boundary layer, mid-plate"},
            "freestream": {"x_mm": 50.0, "y_mm": 2.4,
                           "note": "top of the domain, outside the boundary layer"},
        },
    }
    (args.out / "dataset_metadata.json").write_text(json.dumps(meta, indent=2, default=float),
                                                    encoding="utf-8")
    mb = lambda p: p.stat().st_size / 1e6
    print(f"[write]  number_density_m3.npy  {mb(args.out/'number_density_m3.npy'):.1f} MB  shape {n_out.shape}")
    print(f"[write]  pressure_Pa.npy        {mb(args.out/'pressure_Pa.npy'):.1f} MB")
    print(f"[write]  dataset_metadata.json")
    print("[done]")


if __name__ == "__main__":
    main()
