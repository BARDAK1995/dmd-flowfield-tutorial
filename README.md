# Portable DMD Tutorial — DSMC flow-field analysis

A small, self-contained toolkit and worked example for analysing 2-D
flow-field snapshot data: **structuring it, visualizing it, taking RMS fields
and point spectra, running a Dynamic Mode Decomposition (DMD), and
reconstructing the flow from a few DMD modes.**

It is a deliberately *general and purely dimensional* distillation of a
hypersonic-boundary-layer analysis pipeline. There is **no** boundary-layer-edge
tracking, no linear-stability-theory overlay, no `R = sqrt(Re·x)` / `F`
non-dimensional axes, no freestream-percent normalization — just physical units
(Hz/kHz, m/s, mm, m⁻³, Pa). Hand it a stack of snapshots plus a one-file
description of the grid and timestep, and everything follows.

---

## The example data, in one paragraph

The shipped dataset is an excerpt of a Direct Simulation Monte Carlo (DSMC)
simulation of a **Mach ≈ 6 nitrogen flow over a flat plate**. A small patch of
the wall (`x ∈ [59.5, 60.5] mm`) oscillates (periodic blowing/suction) at
**250 kHz**, acting as an actuator that excites a travelling instability wave
(a second / "Mack" mode) in the boundary layer. We ship two fields over the
converged, statistically steady tail of the run — **number density** (m⁻³) and
**pressure** (Pa) — on a `60 × Nx` grid. Because the forcing frequency is known,
this is an ideal teaching case: the DMD should recover a dominant mode at
**250 kHz**, and reconstructing from that one conjugate pair already reproduces
the travelling wave.

---

## Folder layout

```
DMD_Tutorial/
├── README.md
├── requirements.txt
├── run_all.py                     # run scripts 00–05 end to end
├── dmdkit/                        # the portable library (the reusable part)
│   ├── dataset.py                 # load snapshots + metadata, windows, dimensionalize
│   ├── fields.py                  # time-mean, RMS fields, fluctuation cubes
│   ├── psd.py                     # point power-spectral density, RMS-from-PSD
│   ├── dmd.py                     # detrend → SVD → DMD → physical modes → reconstruction
│   ├── viz.py                     # plotting + GIF animations (journal-minimum fonts)
│   └── io_utils.py                # CSV / .npz / JSON / README writers
├── scripts/                       # runnable, heavily-commented tutorials
│   ├── 00_inspect_data.py         # how the data is structured & dimensionalized
│   ├── 01_animate_field.py        # movies of the field and its fluctuation
│   ├── 02_rms_fields.py           # where the flow fluctuates (RMS maps)
│   ├── 03_point_psd.py            # spectra: boundary tone vs. freestream floor
│   ├── 04_dmd_analysis.py         # full DMD with structured outputs
│   └── 05_reconstruct_from_modes.py  # rebuild the flow from a few modes
├── tools/
│   └── prepare_tutorial_data.py   # how data/ was produced from a full DSMC case
└── data/                          # the shipped dataset
    ├── number_density_m3.npy      # (n_time, n_y, n_x) float32, m^-3
    ├── pressure_Pa.npy            # (n_time, n_y, n_x) float32, Pa
    └── dataset_metadata.json      # grid (mm), timestep (s), units, forcing, freestream
```

---

## Install & run

```bash
pip install -r requirements.txt
python run_all.py            # everything (a few minutes; writes GIFs)
python run_all.py --quick    # skip the slow animations
```

Or run a single step, e.g.:

```bash
python scripts/00_inspect_data.py
python scripts/04_dmd_analysis.py --field pressure --x-min-mm 40 --x-max-mm 100 --y-max-mm 2.5
python scripts/05_reconstruct_from_modes.py --field pressure --pairs 3 \
       --x-min-mm 40 --x-max-mm 100 --y-max-mm 2.5
```

All outputs land in `outputs/<step>/`.

**Animation format.** Movies are written as **MP4** when `ffmpeg` is available
(small, ~3–5 MB) and fall back to GIF otherwise. The fluctuation movies and the
reconstruction comparison are the visual payoff.

**Data size & download.** The dataset is ~580 MB (two `float32` cubes, 501
snapshots × 60 × 2402) — too large for a git repository, so it is hosted on Box,
**not** committed here:

> **Box folder:** https://uofi.app.box.com/folder/391220154153  (`DMD_Tutorial_Data`)

Download `number_density_m3.npy`, `pressure_Pa.npy`, and `dataset_metadata.json`
into this repo's `data/` directory, then run the scripts. (`dataset_metadata.json`
is also committed in the repo, so you really only need the two `.npy` files.)
The 501-snapshot count is deliberate — it is what gives a *clean* DMD of the
250 kHz tone (fewer snapshots smear the mode across nearby frequencies). With
access to the source DSMC case you can instead regenerate `data/` with
`tools/prepare_tutorial_data.py` (see below).

---

## How the data is structured (the one thing to understand)

A *case* is **one array per field** plus **one JSON file**:

* **field array** — shape `(n_time, n_y, n_x)`, a `.npy` (memory-map friendly).
  Axis 0 is time (snapshots), axes 1–2 are the spatial grid `(y, x)`.
* **`dataset_metadata.json`** — the physical context:
  * `grid` : `x_min_mm, x_max_mm, y_min_mm, y_max_mm, nx, ny` (uniform cell centres).
  * `time` : `dt_seconds` (snapshot spacing), `step_first`, `step_interval`, `n_snapshots`.
  * `fields` : for each field, its `file`, `name`, `symbol`, `units`.
  * `forcing`, `freestream`, `conversion`, `suggested_probes` (context + sanity checks).

From `dt` alone, every frequency is dimensional:

```
fs      = 1 / dt                     sampling frequency
Nyquist = fs / 2                     highest resolvable frequency
df      = fs / N                     PSD frequency resolution
```

That is exactly why DMD/PSD results come out in **kHz**: a snapshot timestep of
`dt = 2.5e-7 s` ⇒ `fs = 4 MHz`, so a 250 kHz wave is sampled ~16×/period and is
comfortably resolved.

To run on **your own data**, write a `dataset_metadata.json` with the same keys
and drop your `.npy` cubes next to it — no code changes needed.

---

## What each script teaches

| Script | Question it answers | Key outputs |
|---|---|---|
| `00_inspect_data` | What's in the data? How does dt → kHz? Is number density truly dimensional? | console report, mean & fluctuation maps, `data_report.json` |
| `01_animate_field` | What does the unsteady flow look like? | `<field>_raw.gif`, `<field>_fluct.gif` |
| `02_rms_fields` | *Where* does it fluctuate? | `rms_<field>.png`, `rms_summary.csv` |
| `03_point_psd` | Tone vs. noise — boundary vs. freestream | `psd_<field>.png`, `psd_summary.csv`, `psd_data.csv` |
| `04_dmd_analysis` | The modes: frequency, wavelength, phase speed | SVD/spectrum/mode figures, `dmd_mode_summary.csv`, `mode_XX.npz`, metadata, README |
| `05_reconstruct_from_modes` | Can a few modes rebuild the flow? | error-vs-pairs, snapshot compare, side-by-side & full-field GIFs |

### RMS convention
The RMS field removes a **least-squares line in time** at each point before
taking the RMS — so a slow residual transient isn't mistaken for unsteadiness.

### PSD convention
One-sided, **density**-scaled (`[signal]²/Hz`), Hann window, linear detrend,
`fs = 1/dt`. The spectrum integrates back to the variance:
`RMS = sqrt(∫ PSD df)`. The boundary probe shows a sharp peak at the forcing
frequency; the freestream probe shows a flat broadband floor (here, DSMC
statistical scatter).

### DMD convention
`detrend → SVD → forward-backward exact DMD (rank 30, optimal amplitudes)`.
For each eigenvalue `λ`:

```
frequency   f = Im(log λ) / (2π·dt)      [Hz]   → reported in kHz
growth rate g = Re(log λ) / dt           [1/s]
wavelength  λ_x  from a spatial FFT of the mode's real part   [mm]
phase speed c = f · λ_x                   [m/s]
```

Positive-frequency modes are ranked by amplitude and exported with their
real-part and amplitude maps, a CSV summary, and `.npz` files.

---

## Reconstruction from a few modes (the headline demo)

Real-valued data produces DMD eigenvalues in **complex-conjugate pairs** (a `+f`
mode and its `−f` twin). To rebuild a *real* field you must keep both members of
each pair. So in `05_reconstruct_from_modes.py`, "reconstruct from the first *N*
modes" means the first *N* positive-frequency modes **and** their conjugates —
i.e. `2N` complex modes that recombine into a real field:

```
field(t) ≈ mean + Σ_{selected pairs}  φ_k · b_k · λ_k^t
```

The toolkit finds the conjugate partner automatically (`DmdResult.reconstruct`).
The script sweeps `1…max_pairs`, reports the relative reconstruction error, and
writes a side-by-side animation (original fluctuation vs. reconstruction). For
this forced case, **1–2 conjugate pairs already reproduce the 250 kHz travelling
wave** and the error drops sharply with the first few pairs.

**Selecting any modes you like.** The reconstruction is not limited to "the first
N" — you can rebuild from *any* subset of positive-frequency mode ranks:

```bash
python scripts/05_reconstruct_from_modes.py --field pressure --pairs 3      # first 3 pairs
python scripts/05_reconstruct_from_modes.py --field pressure --ranks 1 3 5  # exactly modes 1,3,5
```

or call the library directly (ranks are 1-based, by descending amplitude; the
conjugate twin of each is added automatically so the field stays real):

```python
from dmdkit import dataset, dmd
f   = dataset.load_field("data", "pressure").crop(40, 100, None, 2.5)
res = dmd.run_dmd(f.data, f.dt_seconds, f.x_mm, f.y_mm, units=f.units, symbol=f.symbol)

cube = res.reconstruct([1, 3, 5], add_mean=True)      # full field from modes 1,3,5
fluc = res.reconstruct([1],       add_mean=False)      # just the dominant wave's fluctuation
# inspect which modes you're combining:
for r in (1, 3, 5):
    print(res.mode_info(r)["frequency_khz"], "kHz")
```

---

## Provenance / regenerating the data

`tools/prepare_tutorial_data.py` documents exactly how `data/` was produced from
a full DSMC case: it reads the solver input/output (`data.dan`, `code.dan`,
`snap.data`, `surfaceMove.data`, `Prot0.dat`), selects the converged step
window, time-subsamples to a portable size, and converts raw particle **counts**
to **number density** via `n = counts · FNUM / V_cell`. Re-run it to re-trim a
different window or keep more snapshots:

```bash
python tools/prepare_tutorial_data.py --step-start 750000 --step-end 800000 --subsample 5
```

---

## Citation & acknowledgments

This toolkit is a thin, dimensional wrapper around well-established methods and
the **PyDMD** library — please cite the originals. Machine-readable entries are
in [`CITATION.cff`](CITATION.cff) and [`references.bib`](references.bib).

**Dynamic Mode Decomposition.** DMD was introduced for fluid flows by
**Schmid (2010)**; the *exact* DMD variant used here is from **Tu et al. (2014)**,
and the noise-robust **forward–backward** option (`forward_backward=True`) follows
**Dawson et al. (2016)**.

- P. J. Schmid, "Dynamic mode decomposition of numerical and experimental data," *J. Fluid Mech.* **656**, 5–28 (2010). doi:10.1017/S0022112010001217
- J. H. Tu, C. W. Rowley, D. M. Luchtenburg, S. L. Brunton, J. N. Kutz, "On dynamic mode decomposition: theory and applications," *J. Comput. Dyn.* **1**(2), 391–421 (2014). doi:10.3934/jcd.2014.1.391
- S. T. M. Dawson, M. S. Hemati, M. O. Williams, C. W. Rowley, "Characterizing and correcting for the effect of sensor noise in the dynamic mode decomposition," *Exp. Fluids* **57**, 42 (2016). doi:10.1007/s00348-016-2127-7

**PyDMD** (the DMD engine this toolkit calls):

- N. Demo, M. Tezzele, G. Rozza, "PyDMD: Python Dynamic Mode Decomposition," *J. Open Source Softw.* **3**(22), 530 (2018). doi:10.21105/joss.00530
- S. M. Ichinaga, F. Andreuzzi, N. Demo, M. Tezzele, K. Lapo, G. Rozza, S. L. Brunton, J. N. Kutz, "PyDMD: A Python package for robust dynamic mode decomposition," *J. Mach. Learn. Res.* (2024). arXiv:2402.07463

PyDMD is BSD-licensed and is a runtime dependency (`requirements.txt`); this
repository does not redistribute it.

## License

MIT — see [`LICENSE`](LICENSE).
