# Numerical validation

The tested environment is MATLAB R2025b and Julia 1.10.10 on macOS arm64. The Julia dependency manifest records the library versions. README.md gives the test commands, and MATLAB_CODE_GUIDE.md describes the numerical implementation.

## MATLAB numerical checks

`tests/verify_numerics.m` checks the translated invariant Lorenz solution, exact heat forward/backward factors, retention of coefficients beyond the backward viscosity cutoff, full and filtered Voigt feedback, Parseval normalization and the observed/unobserved error decomposition.

NSE checks include an analytic velocity field and nonzero two-mode advection. With velocity=(-psi_y,psi_x) and vorticity=Delta psi, the profile psi=cos(x)+0.2*cos(2y) has unforced inviscid nonlinear derivative -0.24*sin(x)*sin(2y). The Basdevant expression is checked against this analytic result, an independent vorticity-advection evaluation and the simultaneous velocity-sign symmetry of its quadratic products.

All seven KdV variants execute on short resolved windows and satisfy the squared L2 error decomposition. These tests use small problems and are separate from grid/time-step convergence studies.

## Run identity and persistent output

`tests/verify_resume.m` checks that an identical calculation can reuse its saved output. Changes to the input's initial data, forcing or recorded source identity require a new calculation. A cached result must contain input hashes, and a missing input file raises an error. Binary SHA-256 hashing is checked independently.

Julia tests create suite and experiment directories in a child process, then verify and reload the saved output after that process exits. This tests persistence separately from an in-process MAT-file round trip. A command-line heat suite also verifies numerical files and all 27 plot exports for its three profiles.

## Cross-language checks

The Julia numerical suite passes 217 assertions. It compares 31 short experiments against MATLAB state arrays and complete error histories, checks independent numerical identities, and verifies MAT/JSON/CSV output. Details and tolerances are in julia/VALIDATION.md.

Full-resolution cross-language validation covers four experiments listed in julia/VALIDATION.md. The remaining configurations are covered by short trajectory comparisons and smoke suites. All 65 smoke configurations and all 52 saved-curve assets have been exercised.

## Manuscript experiments

Full-resolution MATLAB runs cover the Lorenz parameter families and observation-window cases, twelve transport/Burgers cases, eight NSE configurations and seven KdV panels. Standard viscous KdV BFN is the expected failed case. Incomplete cycles do not supply reconstruction endpoints.

The manuscript reference figures and endpoint tables contain computed numerical data. Supplemental raster panels and display-coordinate archives are labeled separately in PROVENANCE.md. Visual agreement with a raster image is distinct from a numerical state-array comparison.

FIGURE_MANIFEST.csv accounts for 25 numbered figures and 52 plot assets. Each reference asset has FIG, PNG and PDF exports. MATLAB suites check required export coverage, and Julia suites check PNG, PDF and SVG coverage. Coverage verifies the completeness of outputs, not numerical accuracy.
