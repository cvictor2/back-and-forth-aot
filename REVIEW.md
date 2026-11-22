# Reproducing and checking the experiments

FIGURE_INDEX.html displays all 52 manuscript plot assets. FIGURE_MANIFEST.csv maps them to manuscript figures and experiment runners. METHODS.md specifies the equations and numerical settings. MATLAB_CODE_GUIDE.md explains the implementation and the time-stepping paths.

## MATLAB checks

From the repository folder

```matlab
addpath('tests')
verify_numerics
verify_resume
run_nse('smoke')
```

The analytical checks cover signed physical steps, explicit feedback, coefficient retention, norm normalization and an independent nonzero NSE advection calculation. Resume tests verify unchanged-input reuse and invalidation after changes to initial data, forcing or source identity. Smoke suites check smaller/shorter versions of the experiments and their output formats.

## Julia checks

```sh
julia --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=julia julia/test/plot_checks.jl
```

The numerical tests include MATLAB reference comparisons, independent identities and a subprocess check that saved results survive Julia exit. The optional plotting checks cover all saved-curve assets and the three smoke suites. See julia/VALIDATION.md for the tested cases and limits.

## Full experiments and reference data

The three runners use full manuscript settings in paper mode. Each new run saves its numerical settings, diagnostics and completed-cycle states in a separate output folder. Full-resolution PDE reference windows require substantial memory. Standard viscous KdV BFN is expected to fail during backward integration and must remain labeled as failed.

Reference figures and computed endpoint CSVs are supplied for comparison. Supplemental archive coordinates and KdV raster panels are labeled separately in PROVENANCE.md. Changes to a numerical method should be checked against the affected full experiment before replacing reference data.
