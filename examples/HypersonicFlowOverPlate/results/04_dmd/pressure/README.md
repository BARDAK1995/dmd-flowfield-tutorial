# DMD analysis: Pressure

- Field: `pressure` (Pressure, Pa)
- Window: steps 0..999, 1000 snapshots, dt = 1.000e-07 s (fs = 10000 kHz).
- Spatial crop: x in [40.0, 100.0] mm, y in [0.02, 2.48] mm.
- Detrend: mean; DMD: forward-backward exact, rank 30, optimal amplitudes, eigenvalue-sorted by magnitude.
- Snapshot-matrix condition number: 4.276e+13.

## Method
Physical frequency f = Im(log lambda)/(2*pi*dt); growth = Re(log lambda)/dt. Positive-frequency modes are ranked by |amplitude|. Each mode's streamwise wavelength comes from a spatial FFT of its real part, and phase speed c = f*lambda.

## Leading modes
- Mode 01: 279.8 kHz, norm. amp 1.000, lambda 2.86 mm, c 799 m/s
- Mode 02: 323.1 kHz, norm. amp 0.527, lambda 2.40 mm, c 776 m/s
- Mode 03: 244.5 kHz, norm. amp 0.468, lambda 3.16 mm, c 772 m/s
- Mode 04: 213.5 kHz, norm. amp 0.321, lambda 3.75 mm, c 801 m/s
- Mode 05: 293.1 kHz, norm. amp 0.234, lambda 2.61 mm, c 765 m/s
- Mode 06: 264.2 kHz, norm. amp 0.196, lambda 3.16 mm, c 834 m/s
- Mode 07: 223.9 kHz, norm. amp 0.170, lambda 3.33 mm, c 746 m/s
- Mode 08: 375.0 kHz, norm. amp 0.084, lambda 1.88 mm, c 703 m/s

Dominant mode: 279.8 kHz (phase speed 799 m/s).
