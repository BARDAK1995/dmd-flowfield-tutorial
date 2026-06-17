#!/usr/bin/env python3
"""
prepare_tutorial_data.py builds the portable tutorial dataset.

This is the provenance script. It shows exactly how the shipped data/ folder
was produced from a full DSMC case, and it lets anyone regenerate or re-trim it.

What it does:

1. Reads the source case's DSMC input/output files (data.dan, code.dan,
   snap.data, surfaceMove.data, Prot0.dat) to recover the grid, timestep,
   particle weight (FNUM), molecular mass, and forcing/freestream conditions.
2. Loads the full N (particle-count) and P (pressure) snapshot cubes.
3. Selects the converged steady tail (a step window) and optionally subsamples
   in time. The full arrays are gigabytes, so the tutorial only needs a portable
   excerpt that still resolves the physics.
4. Converts particle counts to number density: n = counts * FNUM / V_cell.
5. Writes data/number_density_m3.npy, data/pressure_Pa.npy (float32) and
   data/dataset_metadata.json describing the grid, timestep and units.

Defaults target the 250 kHz, Vn=100 m/s wall-actuator case, last 50k steps,
time-subsampled by 5x (so 201 snapshots, dt = 2.5e-7 s, fs = 4 MHz). The forcing
frequency (250 kHz) sits well below Nyquist (2 MHz), so the excerpt is faithful.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TUTORIAL_DIR = HERE.parent
DEFAULT_SOURCE = (TUTORIAL_DIR.parent / "AVS_CASES" /
                  "flatplate120x3_600k_800k_AVS250kHz_Vn100ms")
DEFAULT_N2_MASS_KG = 4.65e-26


# --------------------------------------------------------------------------- #
# Forgiving text parsers for the DSMC input/output files
# --------------------------------------------------------------------------- #
def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="ignore") if path.exists() else ""


def parse_snap_data(path: Path) -> dict:
    lines = [ln.strip() for ln in _read(path).splitlines() if ln.strip()]
    x0, y0 = (float(v) for v in lines[0].split()[:2])
    x1, y1 = (float(v) for v in lines[1].split()[:2])
    s0, s1 = (int(v) for v in lines[2].split()[:2])
    return {"x0_m": x0, "x1_m": x1, "y0_m": y0, "y1_m": y1,
            "step_marker_start": s0, "step_marker_end": s1}


def parse_data_dan(path: Path) -> dict:
    out = {}
    txt = _read(path).splitlines()
    for i, ln in enumerate(txt):
        p = ln.strip().split()
        if not p:
            continue
        if p[0] == "TAU":
            out["tau_seconds"] = float(p[1])
        elif p[0] == "PFnum":
            out["fnum"] = float(p[1])
        elif p[0] == "U" and len(p) >= 2:
            try:
                out["U_mps"] = float(p[1])
            except ValueError:
                pass
        elif p[0] == "NDens":
            out["number_density_m3"] = float(p[1])
        elif p[0] == "T(K)" and len(p) >= 2 and p[1] == "Ttrn":
            out["T_K"] = float(txt[i + 1].strip().split()[0])
    return out


def parse_code_dan(path: Path) -> dict:
    out = {}
    for ln in _read(path).splitlines():
        p = ln.strip().split()
        if len(p) >= 2 and p[0] in ("MACP", "STEP", "MACS"):
            try:
                out[p[0]] = int(float(p[1]))
            except ValueError:
                pass
    return out


def parse_surface_move(path: Path) -> dict:
    if not path.exists():
        return {}
    lines = [ln.strip() for ln in _read(path).splitlines() if ln.strip()]
    x0, _ = (float(v) for v in lines[0].split()[:2])
    x1, _ = (float(v) for v in lines[1].split()[:2])
    freq_khz = float(lines[2].split()[0])
    vt, vn = (float(v) for v in lines[3].split()[:2])
    return {"x_min_mm": x0 * 1e3, "x_max_mm": x1 * 1e3,
            "frequency_hz": freq_khz * 1e3,
            "tangential_velocity_mps": vt, "normal_velocity_mps": vn}


def _prot0_value(txt: str, pattern: str):
    m = re.search(pattern, txt, flags=re.IGNORECASE)
    return float(m.group(1).replace("D", "E")) if m else None


def parse_prot0(path: Path) -> dict:
    txt = _read(path)
    return {
        "molecular_mass_kg": _prot0_value(txt, r"Molecular masses\s*=\s*([0-9.+\-EeD]+)") or DEFAULT_N2_MASS_KG,
        "gas_constant_J_kgK": _prot0_value(txt, r"Gas constant\s*=\s*([0-9.+\-EeD]+)"),
        "gamma": _prot0_value(txt, r"ratio Cp/Cv\s*=\s*([0-9.+\-EeD]+)"),
        "mass_density_kgm3": _prot0_value(txt, r"Density\s*=\s*([0-9.+\-EeD]+)\s*kg/m\^3"),
        "pressure_Pa": _prot0_value(txt, r"Pressure\s*=\s*([0-9.+\-EeD]+)\s*Pa"),
    }


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, default=DEFAULT_SOURCE,
                    help="source DSMC case root (contains npy_files/ and *.dan)")
    ap.add_argument("--out", type=Path, default=TUTORIAL_DIR / "data",
                    help="output data directory")
    ap.add_argument("--step-start", type=int, default=750000)
    ap.add_argument("--step-end", type=int, default=800000)
    ap.add_argument("--subsample", type=int, default=2,
                    help="keep every Nth snapshot in time (2 gives about 501 snapshots, "
                         "enough for a clean DMD of the 250 kHz tone)")
    ap.add_argument("--x-min-mm", type=float, default=None, help="optional streamwise crop")
    ap.add_argument("--x-max-mm", type=float, default=None)
    ap.add_argument("--depth-m", type=float, default=1.0)
    args = ap.parse_args()

    src = args.source
    npy_root = src / "npy_files"
    snap = parse_snap_data(src / "snap.data")
    data = parse_data_dan(src / "data.dan")
    code = parse_code_dan(src / "code.dan")
    surface = parse_surface_move(src / "surfaceMove.data")
    prot0 = parse_prot0(src / "Prot0.dat")

    tau = data["tau_seconds"]
    macp = code.get("MACP", 50)
    first_step = snap["step_marker_start"] + macp           # step of the first stored snapshot
    interval = macp

    n_counts = np.load(npy_root / "N_snapshots.npy", mmap_mode="r")
    p_pa = np.load(npy_root / "P_snapshots.npy", mmap_mode="r")
    nt, ny, nx = (int(v) for v in n_counts.shape)
    last_step = first_step + (nt - 1) * interval

    # pick the time window (the converged tail)
    if not (first_step <= args.step_start <= args.step_end <= last_step):
        raise SystemExit(f"requested {args.step_start}-{args.step_end} outside "
                         f"available {first_step}-{last_step}")
    i0 = int(np.ceil((args.step_start - first_step) / interval))
    i1 = int(np.floor((args.step_end - first_step) / interval)) + 1
    sel = np.arange(i0, i1, args.subsample)
    shipped_first_step = first_step + i0 * interval
    shipped_interval = interval * args.subsample
    dt_seconds = shipped_interval * tau

    # crop in space along the streamwise direction
    x_edges = np.linspace(snap["x0_m"] * 1e3, snap["x1_m"] * 1e3, nx + 1)
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    xmask = np.ones(nx, dtype=bool)
    if args.x_min_mm is not None:
        xmask &= x_centers >= args.x_min_mm
    if args.x_max_mm is not None:
        xmask &= x_centers <= args.x_max_mm
    xi = np.flatnonzero(xmask)
    x_lo_mm, x_hi_mm = float(x_centers[xi[0]]), float(x_centers[xi[-1]])
    # rebuild a clean, uniform [min, max] grid description for the cropped block
    grid_x_min = snap["x0_m"] * 1e3 + xi[0] * (snap["x1_m"] - snap["x0_m"]) * 1e3 / nx
    grid_x_max = snap["x0_m"] * 1e3 + (xi[-1] + 1) * (snap["x1_m"] - snap["x0_m"]) * 1e3 / nx
    nx_out = xi.size

    # convert counts into number density
    dx_m = (snap["x1_m"] - snap["x0_m"]) / nx
    dy_m = (snap["y1_m"] - snap["y0_m"]) / ny
    cell_volume_m3 = dx_m * dy_m * args.depth_m
    fnum = data["fnum"]
    n2mass = prot0["molecular_mass_kg"]
    counts_to_n = fnum / cell_volume_m3

    args.out.mkdir(parents=True, exist_ok=True)
    print(f"[source] {src.name}")
    print(f"[time]   stored {first_step}-{last_step} step {interval}; "
          f"selecting {args.step_start}-{args.step_end} every {args.subsample} "
          f"-> {sel.size} snapshots, dt={dt_seconds:.3e}s, fs={1/dt_seconds:.3e}Hz")
    print(f"[space]  {ny} x {nx_out} cells (x in [{grid_x_min:.3f},{grid_x_max:.3f}] mm)")
    print(f"[convert] n = counts * {counts_to_n:.6e}  (FNUM={fnum:.4e}, V={cell_volume_m3:.4e} m^3)")

    # build the trimmed arrays
    n_out = (np.asarray(n_counts[sel][:, :, xi], dtype=np.float64) * counts_to_n).astype(np.float32)
    p_out = np.asarray(p_pa[sel][:, :, xi], dtype=np.float32)
    np.save(args.out / "number_density_m3.npy", n_out)
    np.save(args.out / "pressure_Pa.npy", p_out)

    # freestream sanity check (top of the domain, upstream of the actuator)
    y_edges = np.linspace(snap["y0_m"], snap["y1_m"], ny + 1) * 1e3
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    x_centers_out = 0.5 * (np.linspace(grid_x_min, grid_x_max, nx_out + 1)[:-1]
                           + np.linspace(grid_x_min, grid_x_max, nx_out + 1)[1:])
    fs_x_mm = min(30.0, x_centers_out[-1]) if x_centers_out[0] <= 30.0 else x_centers_out[0]
    fx = int(np.argmin(np.abs(x_centers_out - fs_x_mm)))
    fy = ny - 3
    n_fs = float(np.mean(n_out[:, fy, fx]))
    ref = data.get("number_density_m3", 1.2e24)
    print(f"[check]  freestream-ish n ~ {n_fs:.3e} m^-3  (inflow n_inf = {ref:.3e}, "
          f"ratio {n_fs/ref:.3f})")

    mach = None
    if prot0.get("gas_constant_J_kgK") and prot0.get("gamma") and data.get("T_K"):
        a = (prot0["gamma"] * prot0["gas_constant_J_kgK"] * data["T_K"]) ** 0.5
        mach = data["U_mps"] / a

    meta = {
        "title": "Hypersonic flat-plate DSMC: portable DMD tutorial dataset",
        "description": ("Trimmed, time-subsampled, dimensionalized excerpt of a "
                        "Mach ~6 N2 flat-plate boundary layer forced by a localized "
                        "wall actuator. Two fields: number density and pressure."),
        "provenance": {
            "source_case": src.name,
            "original_first_step": first_step, "original_interval_steps": interval,
            "original_tau_seconds": tau, "original_n_snapshots": nt,
            "selected_step_start": int(args.step_start), "selected_step_end": int(args.step_end),
            "time_subsample_factor": int(args.subsample),
            "x_crop_mm": [args.x_min_mm, args.x_max_mm],
        },
        "grid": {"x_min_mm": grid_x_min, "x_max_mm": grid_x_max,
                 "y_min_mm": snap["y0_m"] * 1e3, "y_max_mm": snap["y1_m"] * 1e3,
                 "nx": nx_out, "ny": ny, "note": "uniform cell centres in mm"},
        "time": {"n_snapshots": int(sel.size), "step_first": int(shipped_first_step),
                 "step_interval": int(shipped_interval), "tau_seconds": tau,
                 "dt_seconds": dt_seconds, "sampling_frequency_hz": 1.0 / dt_seconds,
                 "nyquist_hz": 0.5 / dt_seconds},
        "fields": {
            "number_density": {"file": "number_density_m3.npy", "name": "Number density",
                               "symbol": "n", "units": "m^-3", "dtype": "float32"},
            "pressure": {"file": "pressure_Pa.npy", "name": "Pressure",
                         "symbol": "p", "units": "Pa", "dtype": "float32"},
        },
        "forcing": {"type": "localized wall actuator (periodic blowing/suction)", **surface},
        "freestream": {"gas": "N2", "mach": mach, "U_mps": data.get("U_mps"),
                       "T_K": data.get("T_K"), "number_density_m3": ref,
                       "mass_density_kgm3": prot0.get("mass_density_kgm3"),
                       "pressure_Pa": prot0.get("pressure_Pa"),
                       "molecular_mass_kg": n2mass,
                       "gas_constant_J_kgK": prot0.get("gas_constant_J_kgK"),
                       "gamma": prot0.get("gamma")},
        "conversion": {"fnum": fnum, "cell_volume_m3": cell_volume_m3,
                       "counts_to_number_density_factor": counts_to_n,
                       "molecular_mass_kg": n2mass,
                       "note": "n = counts * FNUM / V_cell;  rho = n * molecular_mass"},
        "suggested_probes": {
            "boundary": {"x_mm": 60.8, "y_mm": 0.5,
                         "note": "just downstream of the actuator, inside the boundary layer"},
            "freestream": {"x_mm": float(fs_x_mm), "y_mm": float(y_centers[ny - 3]),
                           "note": "top of the domain, outside the boundary layer"},
        },
    }
    (args.out / "dataset_metadata.json").write_text(
        json.dumps(meta, indent=2, default=float), encoding="utf-8")

    mb = lambda p: p.stat().st_size / 1e6
    print(f"[write]  number_density_m3.npy  {mb(args.out/'number_density_m3.npy'):.1f} MB")
    print(f"[write]  pressure_Pa.npy        {mb(args.out/'pressure_Pa.npy'):.1f} MB")
    print(f"[write]  dataset_metadata.json")
    print("[done]")


if __name__ == "__main__":
    main()
