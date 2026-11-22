# Manuscript figure coverage

The manuscript contains 25 numbered figures comprising 52 plot assets. Every reference asset is supplied in FIG, PNG and PDF format under `reference_figures/`.

| Family | Plot assets | MATLAB runner |
| --- | ---: | --- |
| Lorenz | 6 | run_lorenz |
| Transport and Burgers | 36 | run_1d |
| Damped and viscous KdV | 7 | run_1d |
| Two-dimensional Navier–Stokes | 3 | run_nse |
| Total | 52 | |

The three NSE files contain the initial energy spectrum and two stacked observed/unobserved H1-error comparisons. Table 4 summarizes the endpoint comparisons, including full observations. Individual NSE run folders also contain numerical diagnostics and separate error plots.

Heat is available as an optional numerical example and has no manuscript figure. Alpha sensitivity is discussed in the manuscript without a separate plot.

FIGURE_MANIFEST.csv identifies each asset's manuscript figure, runner and reference-export checksums. MATLAB paper and smoke suites check their expected FIG, PDF and PNG files and save `figure_coverage.csv`.

## Julia exports

The three Julia runners provide the same 52 manuscript asset names, with PNG, PDF and SVG exports. There are 6 Lorenz, 43 one-dimensional and 3 NSE assets. The `julia_runner` column in FIGURE_MANIFEST.csv identifies each entry point. Saved-curve exports and fresh numerical runs are checked separately. See julia/VALIDATION.md.
