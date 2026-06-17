# Portable DMD Tutorial for DSMC Flow-Field Analysis

A small, self-contained toolkit and worked example for analysing 2-D flow-field
snapshot data. It shows you how to structure the data, visualize it, take RMS
fields and point spectra, run a Dynamic Mode Decomposition (DMD), and rebuild
the flow from a handful of DMD modes.

It is a deliberately general and purely dimensional distillation of a
hypersonic boundary-layer pipeline. There is no boundary-layer-edge tracking, no
linear-stability-theory overlay, no `R = sqrt(Re·x)` or `F` non-dimensional
axes, and no freestream-percent normalization. Everything stays in physical
units (Hz and kHz, m/s, mm, m⁻³, Pa). You hand it a stack of snapshots plus a
one-file description of the grid and the timestep, and everything else follows.

## The example data, in one paragraph

The shipped dataset is an excerpt of a Direct Simulation Monte Carlo (DSMC)
simulation of a Mach 6 nitrogen flow over a flat plate. A small patch of the
wall, from `x = 59.5` to `60.5 mm`, oscillates with periodic blowing and suction
at 250 kHz. It acts as an actuator that excites a travelling instability wave,
the second (Mack) mode, in the boundary layer. We ship two fields over the
converged, statistically steady tail of the run: number density in m⁻³ and
pressure in Pa, both on a 60 by 2402 grid. Because the forcing frequency is
known, this is an ideal teaching case. The DMD should recover a dominant mode
right at 250 kHz, and rebuilding the flow from just that one conjugate pair
already reproduces the travelling wave.

## Folder layout

```
DMD_Tutorial/
├── README.md
├── requirements.txt
├── run_all.py                     # runs scripts 00 to 05 end to end
├── dmdkit/                        # the portable library (the reusable part)
│   ├── dataset.py                 # load snapshots and metadata, windows, dimensionalize
│   ├── fields.py                  # time-mean, RMS fields, fluctuation cubes
│   ├── psd.py                     # point power-spectral density, RMS from PSD
│   ├── dmd.py                     # detrend, SVD, DMD, physical modes, reconstruction
│   ├── viz.py                     # plotting and MP4/GIF animations (journal-minimum fonts)
│   └── io_utils.py                # CSV, .npz, JSON, README writers
├── scripts/                       # runnable, heavily commented tutorials
│   ├── 00_inspect_data.py         # how the data is structured and dimensionalized
│   ├── 01_animate_field.py        # movies of the field and its fluctuation
│   ├── 02_rms_fields.py           # where the flow fluctuates (RMS maps)
│   ├── 03_point_psd.py            # spectra: boundary tone vs. freestream floor
│   ├── 04_dmd_analysis.py         # full DMD with structured outputs
│   └── 05_reconstruct_from_modes.py  # rebuild the flow from a few modes
├── tools/
│   └── prepare_tutorial_data.py   # how data/ was produced from a full DSMC case
└── data/                          # the shipped dataset (versioned with Git LFS)
    ├── number_density_m3.npy      # (n_time, n_y, n_x) float32, m⁻³
    ├── pressure_Pa.npy            # (n_time, n_y, n_x) float32, Pa
    └── dataset_metadata.json      # grid (mm), timestep (s), units, forcing, freestream
```

## Install and run

The data arrays live in the repo through Git LFS, so install Git LFS once and a
normal clone pulls them down with everything else:

```bash
git lfs install
git clone https://github.com/BARDAK1995/dmd-flowfield-tutorial.git
cd dmd-flowfield-tutorial
pip install -r requirements.txt
python run_all.py            # everything (a few minutes; writes animations)
python run_all.py --quick    # skip the slow animations
```

You can also run a single step:

```bash
python scripts/00_inspect_data.py
python scripts/04_dmd_analysis.py --field pressure --x-min-mm 40 --x-max-mm 100 --y-max-mm 2.5
python scripts/05_reconstruct_from_modes.py --field pressure --pairs 3 \
       --x-min-mm 40 --x-max-mm 100 --y-max-mm 2.5
```

All outputs land in `outputs/<step>/`.

**Animation format.** Movies are written as MP4 when `ffmpeg` is available (about
3 to 5 MB each), and fall back to GIF otherwise. The fluctuation movies and the
reconstruction comparison are the visual payoff.

**The data.** The dataset is about 578 MB total, two `float32` cubes of 501
snapshots each on a 60 by 2402 grid. It is versioned with Git LFS rather than
committed as plain blobs, which is why you need Git LFS installed before
cloning. A copy also lives on Box at
https://uofi.app.box.com/folder/391220154153 if you ever want to grab the files
directly. The 501-snapshot count is deliberate. It is the amount that gives a
clean DMD of the 250 kHz tone, because fewer snapshots smear the mode across
nearby frequencies. If you have the source DSMC case you can regenerate `data/`
yourself with `tools/prepare_tutorial_data.py`, described at the end.

## How the data is structured (the one thing to understand)

A case is one array per field plus one JSON file.

* The **field array** has shape `(n_time, n_y, n_x)` and is a `.npy` that loads
  memory-mapped. Axis 0 is time (the snapshots), and axes 1 and 2 are the
  spatial grid `(y, x)`.
* **`dataset_metadata.json`** carries the physical context:
  * `grid`: `x_min_mm, x_max_mm, y_min_mm, y_max_mm, nx, ny` for uniform cell centres.
  * `time`: `dt_seconds` (the snapshot spacing), `step_first`, `step_interval`, `n_snapshots`.
  * `fields`: for each field, its `file`, `name`, `symbol`, and `units`.
  * `forcing`, `freestream`, `conversion`, and `suggested_probes` for context and sanity checks.

From `dt` alone, every frequency becomes dimensional:

```
fs      = 1 / dt                     sampling frequency
Nyquist = fs / 2                     highest resolvable frequency
df      = fs / N                     PSD frequency resolution
```

That is exactly why the DMD and PSD results come out in kHz. A snapshot timestep
of `dt = 1e-7 s` gives `fs = 10 MHz`, so a 250 kHz wave is sampled about 40 times
per period and is resolved with plenty of room to spare.

To run on your own data, write a `dataset_metadata.json` with the same keys and
drop your `.npy` cubes next to it. No code changes are needed.

## What each script teaches

| Script | Question it answers | Key outputs |
|---|---|---|
| `00_inspect_data` | What is in the data? How does dt become kHz? Is number density truly dimensional? | console report, mean and fluctuation maps, `data_report.json` |
| `01_animate_field` | What does the unsteady flow look like? | `<field>_raw` and `<field>_fluct` movies |
| `02_rms_fields` | Where does it fluctuate? | `rms_<field>.png`, `rms_summary.csv` |
| `03_point_psd` | Tone vs. noise, boundary vs. freestream | `psd_<field>.png`, `psd_summary.csv`, `psd_data.csv` |
| `04_dmd_analysis` | The modes: frequency, wavelength, phase speed | SVD, spectrum, and mode figures, `dmd_mode_summary.csv`, `mode_XX.npz`, metadata, README |
| `05_reconstruct_from_modes` | Can a few modes rebuild the flow? | error vs. pairs, snapshot comparison, side-by-side and full-field movies |

### RMS convention
The RMS field removes a least-squares line in time at each point before taking
the RMS, so a slow residual transient is not mistaken for unsteadiness.

### PSD convention
One-sided, density-scaled (`[signal]²/Hz`), with a Hann window and a linear
detrend, at `fs = 1/dt`. The spectrum integrates back to the variance, so
`RMS = sqrt(∫ PSD df)`. The boundary probe shows a sharp peak at the forcing
frequency, while the freestream probe shows a flat broadband floor, which here
is DSMC statistical scatter.

### DMD convention
The chain is detrend, then SVD, then forward-backward exact DMD at rank 30 with
optimal amplitudes. For each eigenvalue `λ`:

```
frequency   f = Im(log λ) / (2π·dt)      [Hz], reported in kHz
growth rate g = Re(log λ) / dt           [1/s]
wavelength  λ_x  from a spatial FFT of the mode's real part   [mm]
phase speed c = f · λ_x                   [m/s]
```

Positive-frequency modes are ranked by amplitude and exported with their
real-part and amplitude maps, a CSV summary, and `.npz` files.

## Reconstruction from a few modes (the headline demo)

Real-valued data produces DMD eigenvalues in complex-conjugate pairs, a `+f`
mode and its `−f` twin. To rebuild a real field you have to keep both members of
each pair. So when `05_reconstruct_from_modes.py` says "reconstruct from the
first N modes", it means the first N positive-frequency modes together with
their conjugates, which is 2N complex modes that recombine into a real field:

```
field(t) ≈ mean + Σ over selected pairs of  φ_k · b_k · λ_k^t
```

The toolkit finds the conjugate partner for you (`DmdResult.reconstruct`). The
script sweeps from 1 up to `max_pairs`, reports the relative reconstruction
error, and writes a side-by-side movie of the original fluctuation next to the
reconstruction. For this forced case one or two conjugate pairs already
reproduce the 250 kHz travelling wave, and the error drops sharply over the
first few pairs.

**Selecting any modes you like.** The reconstruction is not limited to the first
N. You can rebuild from any subset of positive-frequency mode ranks:

```bash
python scripts/05_reconstruct_from_modes.py --field pressure --pairs 3      # first 3 pairs
python scripts/05_reconstruct_from_modes.py --field pressure --ranks 1 3 5  # exactly modes 1, 3, 5
```

You can also call the library directly. Ranks are 1-based and ordered by
descending amplitude, and the conjugate twin of each is added automatically so
the field stays real:

```python
from dmdkit import dataset, dmd
f   = dataset.load_field("data", "pressure").crop(40, 100, None, 2.5)
res = dmd.run_dmd(f.data, f.dt_seconds, f.x_mm, f.y_mm, units=f.units, symbol=f.symbol)

cube = res.reconstruct([1, 3, 5], add_mean=True)      # full field from modes 1, 3, 5
fluc = res.reconstruct([1],       add_mean=False)     # just the dominant wave's fluctuation
for r in (1, 3, 5):                                   # see which modes you are combining
    print(res.mode_info(r)["frequency_khz"], "kHz")
```

## Provenance and regenerating the data

`tools/prepare_tutorial_data.py` records exactly how `data/` was produced from a
full DSMC case. It reads the solver input and output (`data.dan`, `code.dan`,
`snap.data`, `surfaceMove.data`, `Prot0.dat`), selects the converged step
window, subsamples in time down to a portable size, and converts raw particle
counts to number density with `n = counts · FNUM / V_cell`. Re-run it to pick a
different window or keep more snapshots:

```bash
python tools/prepare_tutorial_data.py --step-start 750000 --step-end 800000 --subsample 2
```

## Citation and acknowledgments

This toolkit is a thin, dimensional wrapper around well-established methods and
the PyDMD library, so please cite the originals. Machine-readable entries are in
[`CITATION.cff`](CITATION.cff) and [`references.bib`](references.bib).

DMD was introduced for fluid flows by Schmid (2010). The exact DMD variant used
here comes from Tu et al. (2014), and the noise-robust forward-backward option
(`forward_backward=True`) follows Dawson et al. (2016).

- P. J. Schmid, "Dynamic mode decomposition of numerical and experimental data," *J. Fluid Mech.* **656**, 5 to 28 (2010). doi:10.1017/S0022112010001217
- J. H. Tu, C. W. Rowley, D. M. Luchtenburg, S. L. Brunton, J. N. Kutz, "On dynamic mode decomposition: theory and applications," *J. Comput. Dyn.* **1**(2), 391 to 421 (2014). doi:10.3934/jcd.2014.1.391
- S. T. M. Dawson, M. S. Hemati, M. O. Williams, C. W. Rowley, "Characterizing and correcting for the effect of sensor noise in the dynamic mode decomposition," *Exp. Fluids* **57**, 42 (2016). doi:10.1007/s00348-016-2127-7

PyDMD is the DMD engine this toolkit calls:

- N. Demo, M. Tezzele, G. Rozza, "PyDMD: Python Dynamic Mode Decomposition," *J. Open Source Softw.* **3**(22), 530 (2018). doi:10.21105/joss.00530
- S. M. Ichinaga, F. Andreuzzi, N. Demo, M. Tezzele, K. Lapo, G. Rozza, S. L. Brunton, J. N. Kutz, "PyDMD: A Python package for robust dynamic mode decomposition," *J. Mach. Learn. Res.* (2024). arXiv:2402.07463

PyDMD is BSD-licensed and is a runtime dependency listed in `requirements.txt`.
This repository does not redistribute it.

## License

MIT, see [`LICENSE`](LICENSE).
