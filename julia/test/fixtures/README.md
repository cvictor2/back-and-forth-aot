# MATLAB reference fixtures

`matlab_short_runs.mat` was generated with MATLAB R2025b from the MATLAB implementation shipped in this repository on 21 November 2026. It contains 31 short experiments, not sampled images or interpolated plot coordinates.

- Twelve heat, transport and Burgers cases cover standard, cutoff, diffusive and Voigt backward models on 128 points.
- Four Lorenz cases cover continuous observations and third-window observations with half-open intervals with fractions 0, 1/2 and 1.
- Eight NSE cases cover standard, diffusive and Voigt models, filtered variants and full observations. They use a normalized Fourier restriction of the supplied restart to 64 by 64 points.
- Seven KdV cases cover every paper panel on 128 points.

Each record stores the complete numerical choices as JSON, initial reference, terminal reference, final assimilated state, six error histories and completed-cycle endpoint errors. The test uses explicit absolute and relative tolerances, including an absolute tolerance near roundoff. Primary errors are not floored for display.

To regenerate with MATLAB, add `julia/test` to the MATLAB path and call `generate_matlab_fixtures`. Its output is this small fixture file. Intermediate reference runs go to a temporary directory. Normal Julia tests only read the supplied fixture and require no MATLAB installation.

NSE states use the convention velocity=(-psi_y,psi_x), vorticity=Delta psi. The fixture states are streamfunction arrays, and the diagnostics are velocity norms.
