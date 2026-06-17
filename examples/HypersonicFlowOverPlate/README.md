# Example: hypersonic flow over a flat plate (second-mode reference case)

This is a complete worked example for the toolkit. It holds the data for one
case and all of the analysis results produced from it, so you can browse the
outputs here directly or regenerate them yourself.

## The case

A Direct Simulation Monte Carlo (DSMC) run of a Mach 6 nitrogen flow over a flat
plate. This is the unforced reference case that the second-mode analysis in the
paper was built on. There is no wall actuator and no external forcing, so the
boundary layer develops a naturally growing second-mode (Mack mode) instability,
the dominant transition mechanism at hypersonic speeds. The point of the example
is to recover that instability straight from the data.

The data in `data/` is the last 1000 snapshots of the run, the converged and
statistically steady tail, calibrated to physical units against the validated
freestream of this case.

| Property | Value |
|----------|-------|
| Fields | number density (m⁻³) and pressure (Pa) |
| Array shape | (1000, 63, 1501) = (time, y, x) |
| Domain | x from 40 to 100 mm, y from 0 to 2.5 mm |
| Snapshots | 1000 (the last 1000 of the run) |
| Timestep | dt = 1e-7 s, so fs = 10 MHz and Nyquist = 5 MHz |
| Freestream | Mach 6 N2, U = 860.6 m/s, T = 50 K, n = 1.2e24 m⁻³, p = 828 Pa |
| Forcing | none (natural second-mode instability) |

Number density was calibrated from raw particle counts and pressure from the raw
sampled field, each anchored to the validated freestream (n = 1.2e24 m⁻³ and
p = 828 Pa). The freestream check in `00_inspect` reads 1.198e24 m⁻³, within
0.2 percent.

## How the results are organized

```
results/
├── 00_inspect/         data report, time-mean and fluctuation maps
├── 01_animation/       movies of each field and its fluctuation
├── 02_rms/             RMS fluctuation maps and a summary table
├── 03_psd/             boundary vs. freestream spectra and tables
├── 04_dmd/
│   ├── pressure/       SVD, spectrum, mode maps, mode .npz, summary CSV, metadata
│   └── number_density/ the same set for number density
└── 05_reconstruction/
    └── pressure/       reconstruction from a few modes, convergence, movies
```

## What the analysis finds

Because the instability is natural rather than driven, it is broadband. The DMD
of pressure picks up a cluster of second-mode frequencies from roughly 210 to
325 kHz, with the most energetic mode near 280 kHz. The paper's representative
second mode near 247 kHz sits inside this band. Every mode in the cluster has a
phase speed close to 800 m/s, about 0.9 of the freestream velocity, and a
streamwise wavelength near 3 mm, which are the hallmarks of the second mode.

![Dominant pressure mode near 280 kHz](results/04_dmd/pressure/modes/mode_01.png)

The mode is concentrated near the wall and grows as it travels downstream, the
classic picture of a spatially amplifying second-mode wavepacket. Number density
shows the same band and the same phase speeds.

The point spectra contrast a probe inside the boundary layer with one in the
free stream. The boundary probe carries the second-mode energy across its band,
while the freestream probe is a flat broadband floor, which here is DSMC
statistical scatter.

![Pressure PSD, boundary vs. freestream](results/03_psd/figures/psd_pressure.png)

Reconstructing the pressure field from only the first few conjugate mode-pairs
recovers the travelling wave. The raw fluctuation is very noisy because in an
unforced DSMC run most of the instantaneous fluctuation energy is incoherent
statistical noise, so the reconstruction is best read as the coherent instability
extracted from that noise. The error against the full coherent DMD model drops
quickly with the first few pairs.

![Reconstruction from a few modes](results/05_reconstruction/pressure/figures/snapshot_compare.png)

## Regenerating these results

The scripts default to this example, so from the repository root you can run:

```bash
python run_all.py
```

or run any single step, for example:

```bash
python scripts/04_dmd_analysis.py --field pressure
python scripts/05_reconstruct_from_modes.py --field pressure --pairs 3
```

Everything writes back into `results/`. The dataset itself was produced from the
reference run by `tools/prepare_reference_data.py`.
