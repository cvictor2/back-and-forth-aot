# Reading and reviewing the MATLAB implementation

The numerical code is in `+paperbfn/`. MATLAB treats that folder as the `paperbfn` namespace. A call such as `paperbfn.step(...)` executes `+paperbfn/step.m`. The three top-level runners select experiments. The equations and integrators are in the helpers, with detailed comments beside the actual operations.

## Follow one NSE experiment

Read the files in this order.

| File | What to inspect |
| --- | --- |
| `run_nse.m` | Public modes and a map of the solver files |
| `+paperbfn/config.m` | Default NSE parameters and definitions of configuration fields |
| `+paperbfn/suite_nse.m` | The eight experiment choices and the explicit smoke-mode overrides |
| `+paperbfn/model.m` | Fourier ordering, derivative multipliers, masks, viscosity, restart and forcing |
| `+paperbfn/reference.m` | Unnudged physical reference at every discrete time in one window |
| `+paperbfn/run.m` | Cycle loops, signed steps, reference indices and endpoint recording |
| `+paperbfn/step.m` | Physical integrating-factor step, backward modification and explicit nudging |
| `+paperbfn/nonlinear.m` | Velocities and the Basdevant streamfunction nonlinearity |
| `+paperbfn/errors.m` | Velocity L2 and H1-seminorm errors and their observation split |

The configuration structure is called `c`. The derived model arrays are in `m`. The state evolved for NSE is the Fourier streamfunction. `v` in `run.m` means the assimilated streamfunction estimate. The names `u` and `v` inside `nonlinear.m` instead mean the two local physical velocity fields. They are local names in different functions.

The velocity convention is `u = -psi_y`, `v = psi_x`, with vorticity `Delta*psi`. The Fourier multipliers are `m.dx = i*kx` and `m.dy = i*ky`, with zero Nyquist first derivatives. The streamfunction state is transformed to physical velocity before forming the quadratic products. The inverse Laplacian converts the resulting expression to the streamfunction right-hand side. The saved forcing is `curl(f)=partial_x f_2-partial_y f_1`, and the equation is `psi_t + Delta^(-1)(velocity dot grad(vorticity)) = nu*Delta*psi + Delta^(-1)*curl(f)`. This convention is used in both languages.

### The actual NSE time step

The physical decomposition used by the solver is

```text
psi_t = L .* psi + N(psi)
L = -nu * (kx^2 + ky^2)
```

Here `N` includes both the nonlinear term and the saved forcing in streamfunction units. For the default forward NSE step, `step.m` evaluates

```matlab
E = exp(h*L);
physical = E.*(z + h*paperbfn.nonlinear(z,m,c));
feedback = h*gain*obs.*(uLeft-z);
next = physical + feedback;
```

`z` and `uLeft` belong to the same old physical time. The reference at the new time is used to measure error after the step, not to form feedback. This is integrating-factor Euler for the physical PDE with a separate explicit observation correction. It is not an NSE RK4 implementation.

For a backward pass, `h=-dt` and `gain=-muBackward`. Their product is positive, so feedback still attracts the estimate toward the reference. The physical backward variant determines `L` and the additional multiplier `V`.

| Backward variant | Multiplier and interpretation |
| --- | --- |
| Standard BFN | `L=-nu*k^2`; negative `h` amplifies high modes |
| Default diffusive BFN | `L=+nu*k^2`; negative `h` gives diffusion |
| Voigt BFN | `V=1/(1+alpha^2*k^2)` multiplies the physical linear term, nonlinear term and feedback |
| Filtered variants | Apply the regularization only on the complement of the observed modes |

Voigt reduces the high-wavenumber backward growth rate. It does not reverse the viscous sign in the way diffusive BFN does. The comments in `step.m` explain these distinctions beside the corresponding expressions.

### Three different masks/cutoffs

`m.obs` selects where data enter feedback. NSE uses the strict circular band `0 < radius < observationCutoff`, or all nonzero modes in the full-observation cases.

`m.keep` selects the inputs used to evaluate nonlinear products. NSE uses a rectangular per-coordinate cutoff. The local multiplication in `nonlinear.m` does not project the caller's evolving state.

The backward viscosity cutoff used in the 1D viscous examples alters only the viscous multiplier. It does not remove state coefficients. These three roles should not be interchanged.

### Physical time and plotted time

Reference column `j+1` belongs to physical time `j*dt`. A forward pass runs from physical zero to `window`. The backward pass runs from `window` to zero. The estimate is carried between passes.

The plotted completed-cycle coordinate increases throughout both passes. Each pass takes half a cycle, so a completed backward endpoint has an integer cycle coordinate and physical time zero. History errors compare against the reference at the current physical time. Endpoint table errors compare against the initial reference.

## Inspect a single update in MATLAB

Open `+paperbfn/step.m` in the MATLAB editor and place a breakpoint on the `physical=...` line in its `ifeuler` branch. From the repository folder, run

```matlab
run_nse('smoke')
```

The first calls construct the reference with `truthOnly=true` and `gain=0`. Later calls perform assimilation. Inspect `z`, `h`, `L`, `V`, `obs`, `uLeft`, `physical` and `feedback`. Step into `paperbfn.nonlinear` to inspect the physical velocities and Fourier products.

To inspect one chosen model without the complete eight-run suite, configure a small NSE experiment explicitly. This uses the same restricted restart as smoke mode.

```matlab
c = paperbfn.config('nse');

% Explicitly restrict the supplied target/forcing for a short review run.
a = load(c.initial.restartFile,'psi_hat_initial','f_hat');
n = 64;
N = size(a.psi_hat_initial,1);
take = [1:n/2, N-n/2+1:N];
psi_hat_initial = a.psi_hat_initial(take,take)*(n/N)^2;
f_hat = a.f_hat(take,take)*(n/N)^2;
restartFile = [tempname '.mat'];
save(restartFile,'psi_hat_initial','f_hat');

c.N = n;
c.nonlinearCutoff = n/3;
c.observe.cutoff = 12;
c.initial.restartFile = restartFile;
c.window = .004;
c.dt = .001;
c.cycles = 2;
c.backward.scheme = 'voigt';   % 'standard', 'diffusive', or 'voigt'
c.backward.alpha = 1e-3;
c.backward.filtered = false;
[result, folder] = paperbfn.run(c);
```

These are short review settings. Full paper settings remain in `config.m` and `suite_nse.m`.

## Lorenz and KdV

For Lorenz assimilation, read `lorenz_leg.m`. It contains the ordinary physical RK4 stages, old-time explicit feedback and half-open observation windows. Its `rhs` uses the shifted variable `z=Z-(rho+sigma)`, explaining the constant term in the third equation. `suite_lorenz.m` specifies the parameter-sweep baselines and beta labels.

For KdV, read `kdv_run.m` from top to bottom. It is a separate driver with configuration, Fourier setup, reference generation, pass loops and endpoint output in one file. The local `physicalStep` helper contains integrating-factor RK4, and `rhs` contains the physical-space square and Fourier derivative. The feedback line is separate and uses an additional integrating-factor weighting. RK4 describes the physical stages, not a fourth-order accuracy claim for the complete method with its explicit feedback correction.

## Verification

From the repository folder, run

```matlab
addpath('tests')
verify_numerics
verify_resume
```

The comments beside each assertion explain its analytical purpose. The checks cover signed heat steps, retention of out-of-band coefficients, full/filtered Voigt feedback, Parseval normalization, an analytic NSE velocity field, and all seven short KdV variants.

Archive mode redraws historical data and cannot verify the solver. Resume compares the full configuration, the current MATLAB source hashes, and SHA-256 hashes of the input files. Replacing the NSE restart or forcing at the same filename invalidates reuse. A cached result must include an input manifest to be eligible for reuse. New runs record both source and input manifests. Reference cache identities include the settings, input arrays and source text.
