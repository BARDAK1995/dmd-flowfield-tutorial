# DMD analysis: Number density

- Field: `number_density` (Number density, m^-3)
- Window: steps 0..999, 1000 snapshots, dt = 1.000e-07 s (fs = 10000 kHz).
- Spatial crop: x in [40.0, 100.0] mm, y in [0.02, 2.48] mm.
- Detrend: mean; DMD: forward-backward exact, rank 30, optimal amplitudes, eigenvalue-sorted by magnitude.
- Snapshot-matrix condition number: 6.762e+13.

## Method
Physical frequency f = Im(log lambda)/(2*pi*dt); growth = Re(log lambda)/dt. Positive-frequency modes are ranked by |amplitude|. Each mode's streamwise wavelength comes from a spatial FFT of its real part, and phase speed c = f*lambda.

## Leading modes
- Mode 01: 280.5 kHz, norm. amp 1.000, lambda 2.86 mm, c 801 m/s
- Mode 02: 415.6 kHz, norm. amp 0.670, lambda 1.88 mm, c 779 m/s
- Mode 03: 246.5 kHz, norm. amp 0.669, lambda 3.16 mm, c 778 m/s
- Mode 04: 325.1 kHz, norm. amp 0.637, lambda 2.40 mm, c 780 m/s
- Mode 05: 449.6 kHz, norm. amp 0.612, lambda 1.76 mm, c 793 m/s
- Mode 06: 210.0 kHz, norm. amp 0.468, lambda 3.75 mm, c 787 m/s
- Mode 07: 338.1 kHz, norm. amp 0.256, lambda 2.14 mm, c 725 m/s
- Mode 08: 264.9 kHz, norm. amp 0.179, lambda 3.16 mm, c 836 m/s

Dominant mode: 280.5 kHz (phase speed 801 m/s).
