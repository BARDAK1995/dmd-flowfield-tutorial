# Agent context for this repository

Claude Code reads this file automatically when it works in a clone of this repo,
so any agent on any machine starts oriented. Read the top-level `README.md` and
`examples/HypersonicFlowOverPlate/README.md` for the full story; this file is the
quick brief and the working conventions.

## What this repo is

A small, portable, purely dimensional toolkit for analysing 2-D flow-field
snapshot data: structuring it, visualizing it, taking RMS fields and point
spectra, running a Dynamic Mode Decomposition (DMD), and rebuilding the flow from
a few DMD modes. It was distilled from a hypersonic boundary-layer pipeline, with
all case-specific scaffolding removed.

## Layout

- `dmdkit/` is the reusable library: `dataset`, `fields`, `psd`, `dmd`, `viz`, `io_utils`.
- `scripts/` are the numbered tutorials, 00 through 05.
- `tools/prepare_tutorial_data.py` shows how the data was produced from a full DSMC case.
- `examples/HypersonicFlowOverPlate/` is a complete worked example: `data/` plus a
  full set of committed `results/`. The scripts default to this example.

## How a case is defined

A case is one `.npy` array per field with shape `(n_time, n_y, n_x)`, plus a
`dataset_metadata.json` that gives the grid in mm, the timestep in seconds, the
units, and some context. Everything physical follows from those two things. To
run on a new case, write a `dataset_metadata.json` with the same keys and drop the
arrays next to it, then point `--data` and `--out` at it. No code changes needed.

## Conventions to keep

- Data is NPY snapshot cubes. Ignore any legacy parquet probe files.
- `N_snapshots`-style arrays hold particle counts; number density is `counts * FNUM / V_cell`.
- Frequencies always come from the real timestep: `f = Im(log lambda) / (2 pi dt)`.
  Never hard-code a frequency scaling factor tied to one dataset's dt.
- The DMD is forward-backward exact at rank 30 with optimal amplitudes. Reconstruction
  keeps complex-conjugate pairs together so the rebuilt field stays real.
- Stay purely dimensional (Hz and kHz, m/s, mm, m⁻³, Pa). Do not add boundary-layer-edge
  tracking, stability-theory overlays, R or F non-dimensional axes, or freestream-percent scaling.
- The data arrays and the `.mp4` animations are versioned with Git LFS. Install Git LFS before cloning.
- Animations are written as MP4 when ffmpeg is available, with a GIF fallback.

## Writing style

Write all prose, docstrings, comments, and commit messages in clear, concise,
conversational language that reads the way people talk. Do not use the em-dash or a
double hyphen as sentence punctuation; the only `--` allowed is a real command-line
flag. Use commas, periods, parentheses, and plain words instead.

## Running it

```bash
git lfs install
pip install -r requirements.txt
python run_all.py          # regenerates the whole example into results/
```
