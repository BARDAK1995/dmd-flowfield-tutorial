"""
dmdkit is a small, portable toolkit for analysing 2-D DSMC/CFD snapshot data.

This package is a deliberately general and dimensional distillation of a
hypersonic-boundary-layer analysis pipeline.  It knows nothing about boundary
layers, linear stability theory, or any particular experiment.  All it knows is:

    * how a stack of 2-D field snapshots is laid out in time,
    * how to turn a snapshot index into a physical time or frequency,
    * how to take means, RMS fields and point power-spectra,
    * how to run a Dynamic Mode Decomposition (DMD) and report the modes in
      physical units (Hz or kHz, m/s, mm), and
    * how to reconstruct the flow field from a chosen subset of DMD modes.

Everything is driven by a small ``fieldinputs.dat`` next to the arrays (dt, the x
and y ranges, and the unit), so the same code runs on any case. A richer
``dataset_metadata.json`` sidecar is also supported for extra context.

Modules
-------
dataset   : load snapshot arrays and metadata, select time windows, dimensionalise
fields    : time-mean, RMS fields, fluctuation cubes, animation frames
psd        : point power-spectral density (freestream vs. boundary), RMS from PSD
dmd        : detrend, then SVD, then DMD, then physical frequencies, modes, reconstruction
viz        : publication-quality plotting helpers (fonts at or above journal minimums)
io_utils  : structured output writers (CSV summaries, mode .npz, JSON, README)
"""

from . import dataset, fields, psd, dmd, viz, io_utils  # noqa: F401

__all__ = ["dataset", "fields", "psd", "dmd", "viz", "io_utils"]
__version__ = "1.0.0"
