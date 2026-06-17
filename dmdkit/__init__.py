"""
dmdkit -- a small, portable toolkit for analysing 2-D DSMC/CFD snapshot data.

This package is a deliberately *general and dimensional* distillation of a
hypersonic-boundary-layer analysis pipeline.  It knows nothing about boundary
layers, linear stability theory, or any particular experiment.  It only knows:

    * how a stack of 2-D field snapshots is laid out in time,
    * how to turn a snapshot index into a physical time / frequency,
    * how to take means, RMS fields and point power-spectra,
    * how to run a Dynamic Mode Decomposition (DMD) and report the modes in
      *physical* units (Hz / kHz, m/s, mm), and
    * how to reconstruct the flow field from a chosen subset of DMD modes.

Everything is driven by a single ``dataset_metadata.json`` sidecar so the same
code runs on any case as long as that file describes the grid and timestep.

Modules
-------
dataset   : load snapshot arrays + metadata, select time windows, dimensionalise
fields    : time-mean, RMS fields, fluctuation cubes, animation frames
psd        : point power-spectral density (freestream vs. boundary), RMS from PSD
dmd        : detrend -> SVD -> DMD -> physical frequencies -> modes -> reconstruction
viz        : publication-quality plotting helpers (fonts >= journal minimums)
io_utils  : structured output writers (CSV summaries, mode .npz, JSON, README)
"""

from . import dataset, fields, psd, dmd, viz, io_utils  # noqa: F401

__all__ = ["dataset", "fields", "psd", "dmd", "viz", "io_utils"]
__version__ = "1.0.0"
