# Sources and data

The MATLAB and Julia implementations accompany *A Critical Note on the Back-and-Forth Nudging Algorithm for Data Assimilation*, by Aseel Farhat, Edriss S. Titi and Collin Victor. METHODS.md specifies the equations, parameters, discretizations and observation operators.

## Numerical sources

The implementations draw on the Lorenz, one-dimensional, streamfunction and KdV routines identified in SOURCE_MANIFEST.json. The distributed source files are identified in RELEASE_SOURCE_MANIFEST.json. Both languages implement the same numerical formulations and experiment settings.

Lorenz uses translated coordinates, physical RK4 stages and a separate explicit feedback term. The parameter families and observation schedules are listed in METHODS.md. Transport and Burgers use integrating-factor Euler with explicit feedback. Their viscous backward cutoff changes the diffusion multiplier without projecting the state.

NSE evolves streamfunction with velocity=(-psi_y,psi_x), vorticity=Delta psi and f_hat=curl(f). Its nonlinear term uses the Basdevant quadratic formulation. Voigt regularizes the full backward right-hand side, including feedback.

KdV uses integrating-factor RK4 for the physical equation and an integrating-factor-weighted explicit feedback term. The damped and viscous models have separate parameter presets. The viscous diffusion panel cancels viscosity with a correction nu; NSE dBFN instead reverses its sign with a correction 2 nu.

## Supplied datasets

| Location | Contents and use |
| --- | --- |
| `data/NSE_reference_and_forcing.mat` | Initial streamfunction, curl-forcing coefficients and spectral operators for the NSE experiments |
| `reference_figures/` | All 52 manuscript plot assets, each in FIG, PNG and PDF format |
| `data/reference_curves.mat` | Numerical curve coordinates exported from the reference FIG files for Julia saved mode |
| `data/reproduced_endpoints/` | Computed cycle-end errors, including the KdV tables |
| `data/manuscript_plot_data.json` | Supplemental display-coordinate archive for 45 Lorenz, transport/Burgers and NSE assets |
| `data/kdv_original_panels/` | Supplemental KdV raster panels for visual comparison |
| `data/NSE_snapshot_diagnostics.json` | NSE snapshot diagnostics with run identifiers and metadata descriptions |
| `julia/test/fixtures/` | MATLAB reference trajectories and comparison metrics for Julia tests |

The manuscript KdV tables and reference curves are computed from numerical trajectories. The supplemental KdV PNG archive contains raster images without their trajectory arrays. Those images support visual comparison, while the supplied computed data support numerical comparisons.

The NSE forcing is the deterministic forcing used in Gesho, Olson and Titi (2015). Its coefficients follow the Olson–Titi forcing described by Olson and Titi (2008), supplied through `OlsonTitiForcing.m` by Adam Larios. The companion MAT file contains the exact saved forcing array and initialization.

## Output identity

Each MATLAB run records the configuration, current source hashes and input-file hashes. Resume requires agreement of all three. Julia records its configuration and source identity. Output histories and endpoint tables identify the diagnostic convention and the time/cycle coordinate.

Paper mode performs numerical integration. MATLAB archive mode redraws the supplemental plot data, while Julia saved mode redraws the manuscript reference curves. These modes are labeled separately because redrawing data is distinct from testing a numerical solver.

The authored code is distributed under the MIT license. Dependencies retain their respective licenses. Citation information is provided in CITATION.bib.
