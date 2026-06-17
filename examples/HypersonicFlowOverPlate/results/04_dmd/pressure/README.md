# DMD analysis: Pressure

- Field: `pressure` (Pressure, Pa)
- Window: steps 750000..800000, 501 snapshots, dt = 1.000e-07 s (fs = 10000 kHz).
- Spatial crop: x in [40.0, 100.0] mm, y in [0.03, 2.48] mm.
- Detrend: mean; DMD: forward-backward exact, rank 30, optimal amplitudes, eigenvalue-sorted by magnitude.
- Snapshot-matrix condition number: 8.030e+13.

## Method
Physical frequency f = Im(log lambda)/(2*pi*dt); growth = Re(log lambda)/dt. Positive-frequency modes are ranked by |amplitude|. Each mode's streamwise wavelength comes from a spatial FFT of its real part, and phase speed c = f*lambda.

## Leading modes
- Mode 01: 251.3 kHz, norm. amp 1.000, lambda 3.16 mm, c 793 m/s
- Mode 02: 312.0 kHz, norm. amp 0.285, lambda 2.50 mm, c 780 m/s
- Mode 03: 355.0 kHz, norm. amp 0.279, lambda 2.22 mm, c 789 m/s
- Mode 04: 338.7 kHz, norm. amp 0.218, lambda 2.31 mm, c 782 m/s
- Mode 05: 80.0 kHz, norm. amp 0.161, lambda 7.50 mm, c 600 m/s
- Mode 06: 469.9 kHz, norm. amp 0.160, lambda 1.67 mm, c 783 m/s
- Mode 07: 202.2 kHz, norm. amp 0.156, lambda 4.00 mm, c 809 m/s
- Mode 08: 264.5 kHz, norm. amp 0.152, lambda 2.73 mm, c 721 m/s

Dominant mode: 251.3 kHz (phase speed 793 m/s).
