# DMD analysis: Number density

- Field: `number_density` (Number density, m^-3)
- Window: steps 750000..800000, 501 snapshots, dt = 1.000e-07 s (fs = 10000 kHz).
- Spatial crop: x in [40.0, 100.0] mm, y in [0.03, 2.48] mm.
- Detrend: mean; DMD: forward-backward exact, rank 30, optimal amplitudes, eigenvalue-sorted by magnitude.
- Snapshot-matrix condition number: 4.850e+13.

## Method
Physical frequency f = Im(log lambda)/(2*pi*dt); growth = Re(log lambda)/dt. Positive-frequency modes are ranked by |amplitude|. Each mode's streamwise wavelength comes from a spatial FFT of its real part, and phase speed c = f*lambda.

## Leading modes
- Mode 01: 349.5 kHz, norm. amp 1.000, lambda 2.22 mm, c 777 m/s
- Mode 02: 252.7 kHz, norm. amp 0.894, lambda 3.16 mm, c 798 m/s
- Mode 03: 376.2 kHz, norm. amp 0.692, lambda 2.07 mm, c 778 m/s
- Mode 04: 411.9 kHz, norm. amp 0.301, lambda 1.94 mm, c 797 m/s
- Mode 05: 315.2 kHz, norm. amp 0.183, lambda 2.61 mm, c 822 m/s
- Mode 06: 236.8 kHz, norm. amp 0.177, lambda 4.00 mm, c 947 m/s
- Mode 07: 486.5 kHz, norm. amp 0.164, lambda 1.71 mm, c 834 m/s
- Mode 08: 482.0 kHz, norm. amp 0.130, lambda 1.62 mm, c 782 m/s

Dominant mode: 349.5 kHz (phase speed 777 m/s).
