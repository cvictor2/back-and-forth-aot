# Julia experiments

These runners implement the same experiments as the MATLAB code using an open source numerical and plotting stack. MATLAB is not required. Julia 1.10.10 was used for validation. Project.toml and Manifest.toml record the dependencies.

Install [Julia](https://julialang.org/downloads/) first if needed. From the repository root, install the dependencies once.

```sh
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
```

The first load compiles the numerical and plotting packages and takes longer than subsequent runs. The three entry scripts accept the same experiment families as MATLAB.

```sh
julia --project=julia julia/run_lorenz.jl smoke
julia --project=julia julia/run_1d.jl smoke
julia --project=julia julia/run_nse.jl smoke
```

Replace `smoke` with `paper` for the manuscript parameters and full cycle budgets. Lorenz parameter sweeps are substantial computations. Each complete cycle has a forward and a backward pass. Smoke runs deliberately use shorter windows and smaller PDE grids and do not reproduce the paper figures.

```sh
julia --project=julia julia/run_lorenz.jl paper
julia --project=julia julia/run_1d.jl paper transport
julia --project=julia julia/run_1d.jl paper burgers
julia --project=julia julia/run_1d.jl paper kdv
julia --project=julia julia/run_nse.jl paper
```

The default 1D selection is `all`, which includes transport, Burgers, and all seven KdV panels. `heat` is optional. `historical` Lorenz mode selects an alternate schedule with closed right endpoints and backward masks evaluated at t-dt. Paper runs use half-open observation windows and current-time backward feedback.

To redraw the complete manuscript figure set without rerunning simulations, use `saved` mode for each runner. This reads the numerical curves in `data/reference_curves.mat`, exported from the final manuscript MATLAB FIG files. It includes every Lorenz, one-dimensional and NSE manuscript asset. These are numerical coordinates, not estimates from PNG images. Saved mode does not validate the solver or produce new state trajectories.

```sh
julia --project=julia julia/run_lorenz.jl saved
julia --project=julia julia/run_1d.jl saved
julia --project=julia julia/run_nse.jl saved
```

Every plot is exported as PNG, vector PDF, and editable SVG. Numerical results are stored as ordinary arrays in MAT files, with CSV histories and cycle endpoints and JSON settings. Julia does not create MATLAB FIG files. Use the SVG exports or redraw the arrays to edit Julia figures. MATLAB reference FIG files are available in `reference_figures/`.

All output goes under `results/julia/`, with a separate directory for each fresh suite and each experiment. `runs.json` lists the outcomes. `coverage.json` lists the root paper assets. Fresh experiment and suite output directories persist after Julia exits.

An unstable experiment has status `failed` and retains its last state and recorded history. It does not contribute invented cycle endpoints. Standard BFN for viscous KdV is expected to fail during backward integration at paper resolution.

## Change a run

```julia
using BFNPaper
c = preset(:nse)
c.scheme = :voigt
c.alpha = 1e-3
c.filtered = true
r = run_experiment(c)
```

For damped KdV with the paper's diffusive backward model

```julia
c = preset(:kdv_damped)
c.scheme = :diffusive
c.gamma = 0.5
c.diffusion_correction = 1e-5
c.cancel_damping = false
r = run_experiment(c)
```

`Config` exposes the time step, physical window, cycles, viscosity, damping, feedback gains, observation cutoff, Voigt parameter, output stride, and other choices explicitly. Presets are in `src/config.jl`. Paper suites are in `src/suites.jl`. Voigt regularizes both the physical terms and feedback. Spectral masks act on nonlinear evaluations, never by projecting the evolving state. Backward viscosity truncation changes only the viscous multiplier.

For NSE the velocity convention is `(-psi_y, psi_x)` with vorticity `Delta*psi`. The saved `f_hat` is `curl(f)`. MATLAB uses the same convention.

The full PDE reference occupies approximately 0.66–0.86 GB plus working arrays. The memory limit is checked before allocating the trajectory. Runs are sequential. NSE reads the shared restart and forcing in `data/`; only smoke mode explicitly restricts that restart to a smaller grid with its FFT normalization adjusted.

## Data and plotting API

`diagnostics.csv` contains cycle coordinate, physical time, cycle number, direction, step, and total/observed/unobserved L2 and H1-seminorm errors. `cycle_errors.csv` reports only completed cycles, against physical time zero. NSE norms are velocity norms. Lorenz norms are Euclidean. KdV plots display L2 errors only. The primary diagnostics have no display floor.

To redraw a fresh result

```julia
replot_result("results/julia/your_run/result.mat")
```

To rerun the numerical checks

```sh
julia --project=julia -e 'using Pkg; Pkg.test()'
```

See `VALIDATION.md` for the checks actually performed and their limits. The MATLAB reference fixtures are small and included, so the tests do not require a MATLAB installation. Numerical agreement is tested with tolerances. Floating point FFTs and long chaotic trajectories can differ across platforms, so bitwise identity is not assumed.

The authored MATLAB and Julia code is distributed under the repository MIT license. Dependencies retain their respective licenses.
