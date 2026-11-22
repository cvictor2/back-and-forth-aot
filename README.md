# Back-and-forth nudging experiments

MATLAB and Julia implementations for *A Critical Note on the Back-and-Forth Nudging Algorithm for Data Assimilation*, by Aseel Farhat, Edriss S. Titi and Collin Victor. The code is distributed under the MIT license.

The three experiment families cover Lorenz, one-dimensional transport/Burgers/KdV models, and the two-dimensional incompressible Navier–Stokes equations. [METHODS.md](METHODS.md) specifies their equations and parameters. [FIGURE_INDEX.html](FIGURE_INDEX.html) displays all 52 revised-paper plot assets, and [FIGURE_MANIFEST.csv](FIGURE_MANIFEST.csv) maps them to manuscript figures and runners.

## MATLAB

Open this folder in MATLAB R2021a or newer with the JVM enabled. No optional toolboxes or paths into another project are required. The `+paperbfn` folder contains the shared implementation.

```matlab
addpath('tests')
verify_numerics
verify_resume
run_lorenz('smoke')
run_1d('smoke')
run_nse('smoke')
```

Smoke mode uses shorter windows and smaller PDE grids to check installation and exercise each variant. It does not reproduce the full paper trajectories. For the manuscript settings

```matlab
run_lorenz('paper')
run_1d('paper')
run_nse('paper')
```

Select individual 1D families as needed.

```matlab
run_1d('paper','transport')
run_1d('paper','burgers')
run_1d('paper','kdv')
run_1d('paper','heat')       % optional examples, outside the paper figure set
```

[MATLAB_CODE_GUIDE.md](MATLAB_CODE_GUIDE.md) follows an NSE calculation from configuration to its actual time step. Detailed comments beside the operations explain state variables, Fourier conventions, masks, backward regularization and error calculations. Lorenz assimilation has its own loop in `lorenz_leg.m`, and KdV has a self-contained driver in `kdv_run.m`.

## Julia

Julia provides the same experiment families without requiring MATLAB. Install the dependencies from the repository root and run the short checks first.

```sh
julia --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=julia julia/run_lorenz.jl smoke
julia --project=julia julia/run_1d.jl smoke
julia --project=julia julia/run_nse.jl smoke
```

Replace `smoke` with `paper` for full manuscript settings. The 1D runner accepts a second argument such as `transport`, `burgers` or `kdv`. See [julia/README.md](julia/README.md) for the configuration and plotting API.

## Outputs and diagnostics

Fresh runs save configurations, source identities, numerical histories, completed-cycle endpoint states and plots in isolated directories under `results/`. Julia uses `results/julia/`. These directories persist after the numerical program exits and are ignored by Git. MATLAB exports PNG, vector PDF and editable FIG files. Julia exports PNG, PDF and editable SVG files, with ordinary MAT arrays and CSV diagnostics.

One cycle consists of a forward pass followed by a backward pass. Each pass occupies half a cycle on the horizontal axis. History errors compare the current estimate with the reference at the same physical time. Completed-cycle errors compare against the reference at physical time zero.

Primary diagnostics use absolute integral L2 and H1-seminorm errors, velocity norms for NSE, and Euclidean norms for Lorenz. KdV figures display total, observed and unobserved L2 errors. KdV result files also retain a separately labeled legacy diagnostic for source comparison. Primary saved errors are not clipped to a plotting floor.

An unstable run is recorded as failed. Incomplete cycles do not produce recovered endpoints. Standard BFN for viscous KdV is expected to fail during its first backward pass at the paper settings.

Full-resolution reference storage for one PDE window requires approximately 0.66–0.86 GB plus working arrays. Experiments run sequentially. Lorenz paper sweeps include 100-cycle families and take substantially longer than smoke checks.

## Custom experiments and resume

For one NSE experiment in MATLAB

```matlab
c = paperbfn.config('nse');
c.backward.scheme = 'voigt';
c.backward.alpha = 1e-3;
c.backward.filtered = true;
[result, folder] = paperbfn.run(c);
```

For damped KdV with the diffusive backward correction

```matlab
options = struct('model','damped','variant','diffusive BFN', ...
                 'damping',0.5,'stabilization',1e-5);
result = paperbfn.kdv_run(options);
```

The optional MATLAB resume shortcut reuses only completed results with matching settings, source hashes and input-file hashes. Changing the NSE initial state or forcing at the same filename invalidates reuse. Older results without input hashes are recomputed. Input and source identity are recorded in the saved result.

Both languages use the NSE convention velocity=(-psi_y,psi_x), vorticity=Delta psi, with saved curl-forcing coefficients. Nonlinear masks act on product inputs. Backward viscosity truncation changes only the viscous multiplier and retains all state coefficients.

## Saved and historical plots

Julia `saved` mode redraws all final revised-paper numerical curves without integration.

```sh
julia --project=julia julia/run_lorenz.jl saved
julia --project=julia julia/run_1d.jl saved
julia --project=julia julia/run_nse.jl saved
```

MATLAB `archive` mode instead redraws historical display coordinates and copies original KdV raster panels. These predate the final corrected Lorenz window plots. `run_lorenz('historical')` performs a separate numerical audit of the older observation-window convention. Neither historical plots nor saved-curve redraws substitute for a fresh solver test.

## Validation and data

[VALIDATION.md](VALIDATION.md) and [julia/VALIDATION.md](julia/VALIDATION.md) describe the tests performed and the full-resolution comparisons. [PROVENANCE.md](PROVENANCE.md) distinguishes current reference data from historical archives and credits the numerical sources and forcing data. [REVIEW.md](REVIEW.md) provides a reproducibility checklist. Citation information is in CITATION.bib.
