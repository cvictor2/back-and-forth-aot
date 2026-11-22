function c=config(model,scenario)
%CONFIG Paper defaults. All cutoffs are integer Fourier indices.
% Return a configuration structure c. This routine defines DEFAULTS only.
% The suite files override individual values to select each paper experiment.
% KdV settings are separate and appear at the beginning of kdv_run.m.
%
% Review an NSE experiment by reading config('nse'), suite_nse.m, model.m,
% reference.m, run.m, step.m, nonlinear.m and errors.m in that order.
% All times below are physical PDE/ODE times. A plotted BFN cycle is different
% and always contains one forward pass plus one backward pass.
if nargin<2, scenario='mixed'; end
root=fileparts(fileparts(mfilename('fullpath')));
c.name=[model '_' scenario]; c.model=model; c.scenario=scenario;
% Shared 1D baseline. N is points per coordinate (not total points in 2D),
% L is the periodic side length, dt the step, window the physical pass length,
% and cycles the number of forward/backward pairs. Model overrides follow.
c.N=4096; c.L=2; c.dt=1e-4; c.window=1; c.cycles=5;
% Physical coefficients. nu is viscosity, gamma is scalar damping,
% transportSpeed multiplies +u_x, and dispersion is used by the KdV driver.
c.nu=0; c.gamma=0; c.transportSpeed=1; c.dispersion=1;
% ifeuler uses an exact diagonal linear factor and an explicit Euler remainder.
% feedbackFactor=false places the old-time nudging outside that factor.
% The physical integrator and explicit nudging are separate choices.
c.integrator='ifeuler'; c.feedbackFactor=false;
% Cutoffs are mode indices, not physical wavenumbers 2*pi*index/L.
% The nonlinear cutoff masks product inputs. It does not truncate the state.
% NSE has a separate per-coordinate rectangular input mask in model.m.
c.nonlinearCutoff=682; c.nonlinearRule='1d';
% Nonnegative feedback strengths. step.m applies the backward sign internally.
% Their explicit correction alone is stable only if mu*dt <= 2.
c.muForward=1000; c.muBackward=1000;
% PDE observation mask is strictly 0<radius<cutoff. full=true observes every
% nonzero mode. For Lorenz, components picks coordinates and fractions gives
% the observation fraction within each component's third of the time window.
% The rule/includeMean fields record the fixed convention used by model.m.
c.observe=struct('cutoff',16,'rule','strict','includeMean',false, ...
 'components',[1 2],'fractions',[1 1 1],'schedule','continuous','full',false);
% 1D profiles are cos(mode 1), highAmplitude*cos(highMode), or their sum.
% reference/guess are explicit Lorenz states. guessType selects zero,
% observed projection, exact reference, or an explicit vector as appropriate.
% spinup evolves the physical reference before defining assimilation time zero.
% NSE reads restartFile instead of constructing these trigonometric profiles.
c.initial=struct('kind',scenario,'highMode',30,'highAmplitude',0.05, ...
 'reference',[20;30;40],'guess',[30;40;50],'guessType','zero', ...
 'spinup',0.0101,'restartFile','','restartVariable','psi_hat_initial');
% The shared 1D examples have zero forcing. NSE uses the supplied f_hat array,
% not a newly generated forcing from G. Some fields retain source metadata.
% KdV builds its own exp(cos) forcing and explicitly removes its mean.
c.forcing=struct('kind','none','amplitude',1,'G',1e6);
% Lorenz defaults use sigma=10,rho=28,beta=8/3. suite_lorenz.m documents the
% sigma=20 baseline for the rho and beta sweeps. backwardObservationShift
% and windowRule select the alternate observation schedule in historical mode.
c.lorenz=struct('sigma',10,'rho',28,'beta',8/3,'originalArithmetic',true,'windowRule','halfOpen','backwardObservationShift',0);
% standard leaves the physical operator unchanged on backward integration.
% cutoff removes viscosity above viscosityCutoff without deleting state modes.
% diffusive changes the viscous/damping factors on backward passes only.
% voigt applies 1/(1+alpha^2*k^2) to the full backward RHS, including feedback.
% filtered=true restricts regularization to the unobserved complement.
% Default viscosityFactor=-1 reverses the viscous sign in the signed-time RHS.
% The alpha in the multiplier is squared, as in the manuscript convention.
c.backward=struct('scheme','standard','filtered',false,'alpha',1e-3, ...
 'viscosityFactor',-1,'dampingFactor',-1,'viscosityCutoff',50);
% Output stride samples diagnostics only. Every solver step is still taken.
% Cycle states store completed t=0 reconstructions. Reference caching is optional
% and its byte limit applies to one window, not the total cycle count.
c.output=struct('root',fullfile(root,'results'),'plot',true, ...
 'saveCycleStates',true,'diagnosticStride',10);
c.reference=struct('useCache',false,'maxBytes',2e9,'root',fullfile(root,'cache'));
% Stop on nonfinite/very large errors. Empty tolerance requests the full cycle
% budget. A nonempty tolerance is checked at completed backward endpoints.
c.stop=struct('maxNorm',1e100,'tolerance',[]);
switch model
 case 'lorenz'
  c.N=3; c.dt=1e-5; c.window=1; c.cycles=100;
  c.integrator='rk4'; c.initial.spinup=0; c.initial.guessType='explicit';
  c.muForward=100; c.muBackward=100; c.output.diagnosticStride=500;
 case 'heat'
  c.transportSpeed=0; c.nu=0.05; c.backward.scheme='cutoff';
 case 'transport'
 % Burgers uses a short physical window and a small time step.
 % The reference spinup comprises 10001 steps, giving 0.010001 time units.
 case 'burgers'
  c.transportSpeed=0; c.dt=1e-6; c.window=0.01; c.initial.spinup=0.010001;
 % Paper NSE settings. The domain is [0,2*pi)^2 with 1024 points per side.
 % The observed initial projection is the initial estimate, not the full truth.
 % The strict observation radius is 50, and nonlinear inputs use an
 % N/3 per-coordinate mask. Forcing and initial streamfunction are supplied data.
 case 'nse'
  c.N=1024; c.L=2*pi; c.dt=1e-3; c.window=0.05; c.cycles=5;
  c.nu=1e-4; c.transportSpeed=0; c.initial.spinup=0;
  c.initial.guessType='projection'; c.observe.cutoff=50;
  c.nonlinearRule='nse'; c.nonlinearCutoff=c.N/3;
  c.initial.restartFile=fullfile(root,'data','NSE_reference_and_forcing.mat');
  c.forcing.kind='saved'; c.output.diagnosticStride=1;
 otherwise
  error('paperbfn:Model','Unknown model %s.',model);
end
end
