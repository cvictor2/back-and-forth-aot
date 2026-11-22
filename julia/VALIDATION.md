# Julia validation

The tested environment is Julia 1.10.10 on macOS arm64 and MATLAB R2025b. The supplied Manifest.toml pins the dependencies. `Pkg.test()` passes 217 numerical and result-file assertions, including independent NSE advection and output persistence after a subprocess exits. MATLAB is not required to repeat these tests.

## Numerical agreement

Thirty-one short MATLAB experiments were regenerated from the MATLAB code distributed in this repository. Julia was compared against the initial reference, terminal reference, final assimilated state, all six error histories, and completed-cycle errors for each case. This includes twelve heat/transport/Burgers configurations, four Lorenz observation configurations, all eight NSE configurations, and all seven KdV panels.

Across these cases, the largest absolute final-state coefficient difference was 8.11e-15. The largest absolute difference in the six error histories or completed-cycle error tables was 3.56e-15. Individual metrics are in `test/fixtures/parity_metrics.json`. Assertions use explicit absolute and relative tolerances, so near-zero errors do not require misleading relative accuracy.

Independent checks verify the L2 and H1 observed/unobserved decompositions, Parseval normalization, a known NSE velocity field, Nyquist derivative conventions, a nonzero analytic NSE advection field, cutoff retention of out-of-band state coefficients, Voigt factors and filtered observed modes, half-open observation windows including the zero-observation case, parameter/memory guards, and failed-run handling. A result-file round trip verifies MAT, JSON and CSV output.

Four selected cases were also run at the full manuscript grid, time step, physical window and five-cycle budget, and compared with the saved full MATLAB state arrays.

| Full-resolution experiment | Relative final-state difference |
| --- | ---: |
| NSE standard BFN, 1024 × 1024 | 9.97e-19 |
| Viscous Burgers, mixed profile, 4096 points | 5.45e-15 |
| Damped KdV standard BFN, 512 points | 3.34e-15 |
| Damped KdV Voigt BFN, 512 points | 1.11e-14 |

The metric is the Euclidean/Frobenius norm of the difference between the final spectral arrays divided by `max(norm(MATLAB array), 1)`. NSE arrays represent streamfunction. All four completed their five cycles. Summary metrics are in `test/fixtures/full_resolution_metrics.json`. The large full-resolution state files remain outside the distribution package.

## Runner and figure checks

All 65 short-suite experiments completed successfully, with configurations saved for each run. The three suites exported every required paper asset name, with PNG, PDF and SVG files.

| Suite | Experiments | Paper assets |
| --- | ---: | ---: |
| Lorenz | 38 | 6 |
| Transport, Burgers and KdV | 19 | 43 |
| NSE | 8 | 3 |
| Total | 65 | 52 |

All 52 final manuscript numerical curve records were also redrawn through Julia, producing 156 image/vector exports. This is checked separately from running the solver. The records come from the final MATLAB FIG files, not PNG digitization. A saved fresh result was successfully loaded and replotted without integration.

The Julia root export names match FIGURE_MANIFEST.csv. The renderer, colors, text and vector format differ from MATLAB. Exported figures remain editable through SVG and numerical data. No Julia-generated MATLAB FIG file is claimed.

## Repeating the checks

From the repository root

```sh
julia --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=julia julia/test/plot_checks.jl
```

The optional plotting check redraws all saved assets and runs all three smoke suites in temporary folders. It takes longer than the numerical unit checks. The included MATLAB fixture generator can regenerate the small reference arrays when MATLAB is available.

## Limits

Full-resolution cross-language validation covers the four cases listed above. Other configurations are covered by short MATLAB trajectory comparisons and smoke suites.

Long chaotic trajectories, nearly exact reconstructions and unstable backward equations amplify floating point differences. Pixel-for-pixel or bitwise agreement is not an acceptance criterion. In particular standard BFN for viscous KdV can fail at a slightly different step across FFT implementations. Failed runs retain their status and last valid state, and never supply endpoints for incomplete cycles. Full-resolution MATLAB experiment coverage is described in the root VALIDATION.md.

## Output persistence and NSE convention

Saved experiment and suite directories persist after process exit. The parent test reloads a result after its producing subprocess exits, verifying both directory levels and the saved numerical files.

Both languages use velocity=(-psi_y,psi_x), vorticity=Delta psi and f_hat=curl(f). The two-mode advection test checks the sign against a nonzero analytic result and an independent vorticity-advection expression. All 31 MATLAB reference-trajectory comparisons use this convention.
