# Equations, discretization and experiment map

## Common conventions

Fourier arrays use MATLAB's unnormalized FFT. Feedback uses observations at the current time step and the old assimilated state. It is applied separately from the physical integrator. Backward passes use a negative physical step and a negative feedback gain, so their product remains positive. No interpolation of observations to RK stages is performed. This mixed scheme is not claimed to be fourth order in the full nudged equations.

Dealiasing masks act on inputs to the nonlinear evaluation, as in the supplied drivers. They do not reset the evolving state after each step. This distinction matters for reproducing the specified viscosity-only regularization.

## Lorenz

The evolved coordinates are the translated Lorenz coordinates used in the paper

    x' = sigma (y-x)
    y' = -sigma x - y - x z
    z' = -b z + x y - b (rho+sigma)

Physical terms use RK4. The default baseline is sigma=10, rho=28, b=8/3. The time step is 1e-5 and the physical window is one. Generic target and guess are (20,30,40) and (30,40,50). Only the first two components are observed in invariant and parameter tests.

- Invariant tests use target (0,0,0), guesses (0,0,phi) with phi=1, 1e-5 and 1e-10, gain 1000 and five cycles.
- Parameter sweeps use gain 100 and 100 cycles. Sigma takes 10,20,30,40,50,100. Rho takes 10,20,28,30,40,50,100. b takes 1,2,8/3,3,4,5,10. The sigma sweep uses rho=28 and b=8/3. The rho sweep uses sigma=20 and b=8/3. The b sweep uses sigma=20 and rho=28.
- Windowed tests use ten cycles and half-open physical-time component windows [(i-1)/3, ((i-1)+fraction)/3). Gains are 10,100,1000,10000,100000 at fraction one. Fractions are 0,1/8,1/4,1/2,3/4,7/8,15/16,31/32,63/64,1 at gain 10000.

Observation masks are evaluated at the current physical time, including backward time traversal. A fixed 1e-5 grid does not align exactly with thirds. Endpoints are not silently moved and no time-step change is made to force alignment.

## Transport, Burgers and heat

Paper runs use L=2, N=4096, gain 1000 and five cycles. Integer indices 0<|m|<16 are observed. Inputs to the nonlinear term retain nonzero indices through |m|=682. Transport references are first spun up for 101 steps, or 0.0101 time units. Burgers uses 10001 steps, or 0.010001 time units. These counts specify the reference initialization before the assimilation window. The assimilating guess is zero.

| Model | dt | Window | Viscosity | Transport coefficient |
|---|---|---|---|---|
| Inviscid transport | 1e-4 | 1 | 0 | 1 |
| Viscous transport | 1e-4 | 1 | 0.05 | 1 |
| Inviscid Burgers | 1e-6 | 0.01 | 0 | 0 |
| Viscous Burgers | 1e-6 | 0.01 | 0.05 | 0 |

Before spin-up, the three profiles are A cos(30 pi x), cos(pi x), and cos(pi x)+0.05 cos(30 pi x). A=1 for viscous transport and 0.05 otherwise.

The linear physical multiplier is -nu kappa^2 + i a kappa. Thus transport follows the specified direction, with +a u_x. Burgers uses G(z)=-i kappa FFT(real(IFFT(Dz))^2)/2. The update is

    reference_next = exp(h A) (reference + h G(reference))
    estimate_next = exp(h A_s) (estimate + h G(estimate))
                    + h signed_gain P_obs (reference - estimate)

In a viscous backward pass A_s sets only the viscous multiplier to zero for |m|>=50. All state coefficients remain. No integrating-factor weight multiplies feedback in these models. Heat is available separately through the same linear solver, with a=0 and viscosity 0.05. Its optional numerical tests are not manuscript figure claims.

## NSE

The state is the streamfunction FFT. The convention is velocity=(-psi_y,psi_x) and vorticity=Delta psi. The initial state and forcing are loaded from the supplied `data/NSE_reference_and_forcing.mat`. The saved f_hat is curl(f)=partial_x f_2-partial_y f_1, not a velocity forcing array. The evolved equation is psi_t + Delta^(-1)(velocity dot grad(vorticity)) = nu*Delta*psi + Delta^(-1)*curl(f). The nonlinear term is evaluated with the Basdevant quadratic formula.

The domain has side lengths 2pi, N=1024, nu=1e-4, dt=1e-3, window=0.05, gain 1000 and five cycles. Initialization is the observed projection. Partial observations retain 0<|k|<50 with the Euclidean radial cutoff. Full observations include all nonzero grid coefficients. The rectangular nonlinear input mask uses the same index boundaries as Parameters.m. Zero and Nyquist first-derivative conventions follow that source.

Physical stepping is integrating-factor Euler. Standard feedback is added without integrating-factor weighting. The backward variants are

- Standard BFN retains the physical viscosity.
- dBFN adds -2 nu Delta, reversing the viscous sign in physical time.
- VBFN multiplies the entire physical RHS and feedback by 1/(1+alpha^2 |k|^2), with alpha=1e-3.
- Filtered variants apply the modification only in unobserved modes.

Full-observation comparisons run the three unfiltered methods. The output includes velocity L2 and H1-seminorm histories and endpoint tables, spectra and state arrays. No additional alpha figure is generated. Alpha can be edited in the saved explicit config.

For streamfunction coefficient errors, the squared velocity magnitude is (|Dx|^2+|Dy|^2)|error|^2. Its sum is multiplied by L^2/N^4 in two dimensions. The H1 seminorm includes the additional physical |k|^2 factor.

## KdV

`kdv_run.m` is adapted from the validated coefficient runner based on KdV.m, not KdV_refactor.m. All seven panels use L=2, N=512, dt=1e-5, window=1, gain 1000 and five cycles. The exact initial state is cos(pi x)+0.1 cos(20 pi x). Forcing is exp(cos(pi x)) with its zero Fourier coefficient removed. Observations retain 0<|m|<15. The nonlinear input mask retains 0<|m|<=floor(N/3), excluding Nyquist.

The physical integrator is IF-RK4. Its separate feedback is multiplied by the signed-step integrating factor as in KdV.m, unlike transport/Burgers and NSE. Voigt smooths physical terms and feedback. The default dispersion coefficient is one.

| Panel | Forward model | Backward modification |
|---|---|---|
| bfn_normal | scalar damping 0.5 | standard BFN |
| bfn_damped | scalar damping 0.5 | reverse scalar damping using +2 gamma v |
| bfn_diffusive | scalar damping 0.5 | retain scalar damping and add -1e-5 v_xx in physical time |
| bfn_VOIGT | scalar damping 0.5 | alpha=1e-7 |
| bfn_model_diffusive | viscosity 1e-3 | standard BFN, expected instability |
| bfn_viscous_model_diffusive | viscosity 1e-3 | add -nu v_xx, canceling the physical viscosity |
| bfn_VOIGT_model_diffusive | viscosity 1e-3 | alpha=1e-2 |

KdV initialization is exact to test stability and preservation, not recovery from an unknown initial state. Failed cycles are left uncompleted. Error histories use a 0-to-5 completed-cycle coordinate for the five-cycle experiments. Supplemental raster panels retain their own axis labels.

## Alternate Lorenz observation schedule

Paper mode uses current-time half-open observation masks in both directions. The optional `run_lorenz('historical')` schedule uses closed right endpoints and a backward mask evaluated at t-dt. These are different discrete observation operators. The zero-fraction setting in the closed-endpoint schedule can apply feedback at isolated endpoints, whereas the half-open zero-fraction interval is empty. Report the chosen observation schedule when comparing results.
