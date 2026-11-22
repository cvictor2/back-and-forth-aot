function c=validate(c)
%VALIDATE Reject silently changed grids and unstable isolated feedback choices.
% Guard numerical settings before reference allocation or integration.
% This routine checks configuration consistency, not convergence of the PDE
% method. Passing these checks does not guarantee a stable backward problem.
% For example full-resolution viscous KdV BFN is intentionally unstable.
for name={'dt','window','L'}
 validateattributes(c.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
for name={'cycles','N'}
 validateattributes(c.(name{1}),{'numeric'},{'scalar','integer','positive'});
end
for name={'nu','gamma','muForward','muBackward'}
 validateattributes(c.(name{1}),{'numeric'},{'scalar','real','finite','nonnegative'});
end
% Use an integer number of dt steps to land exactly at the requested window.
% A mismatch is rejected rather than changing the final observation time.
% Spinup must likewise correspond to an integer number of steps.
c.stepsPerPass=round(c.window/c.dt);
assert(c.stepsPerPass>=1 && abs(c.stepsPerPass*c.dt-c.window)<1e-10*c.window, ...
 'paperbfn:Grid','window must be an integer multiple of dt.');
assert(abs(c.initial.spinup/c.dt-round(c.initial.spinup/c.dt))<1e-8);
assert(ismember(c.backward.scheme,{'standard','cutoff','diffusive','voigt'}));
assert(ismember(c.integrator,{'rk4','ifeuler','ifrk4'}));
assert(ismember(c.initial.guessType,{'zero','projection','exact','explicit'}));
assert(ismember(c.observe.schedule,{'continuous','thirds'}));
assert(numel(c.observe.fractions)==3 && all(c.observe.fractions>=0 & c.observe.fractions<=1));
validateattributes(c.output.diagnosticStride,{'numeric'},{'scalar','integer','positive'});
% For an isolated observed error e, explicit feedback gives e_new=(1-mu*dt)e.
% Its necessary scalar stability condition is 0<=mu*dt<=2. This bound alone
% does not control nonlinear, advective or backward-diffusion instabilities.
assert(max(c.muForward,c.muBackward)*c.dt<=2, ...
 'paperbfn:Feedback','Explicit feedback requires mu*dt<=2 even in isolation.');
% Lorenz states have three entries. PDE grids must be even for the stored
% FFT ordering, with cutoffs and initial high modes below the Nyquist limit.
if strcmp(c.model,'lorenz')
 assert(c.N==3 && strcmp(c.backward.scheme,'standard'));
 assert(numel(c.initial.reference)==3 && numel(c.initial.guess)==3);
else
 assert(mod(c.N,2)==0 && c.N>=16);
 assert(c.observe.cutoff<c.N/2 && c.observe.cutoff>=0);
 assert(c.nonlinearCutoff<c.N/2);
 assert(c.initial.highMode<c.N/2);
end
% A Voigt run must use a finite positive alpha. alpha is the length parameter
% inside 1+alpha^2*k^2, not the already-squared coefficient.
if strcmp(c.backward.scheme,'voigt')
 assert(c.backward.alpha>0 && isfinite(c.backward.alpha));
end
end
