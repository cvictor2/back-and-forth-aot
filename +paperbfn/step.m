function next=step(z,h,t,uLeft,~,m,c,backward,truthOnly)
%STEP Signed physical update with separate OLD-TIME explicit feedback.
% One signed physical step followed by an explicit observation correction.
% This is the shared PDE stepper for transport, heat, Burgers and NSE.
% Lorenz uses this routine for its reference and lorenz_leg.m for assimilation.
% KdV has its own physicalStep function in kdv_run.m.
%
% Inputs and units
%   z          current assimilated state, in Fourier coefficients for PDEs
%   h          +dt on the forward pass, -dt on the backward pass
%   t          physical time of z, before taking this step
%   uLeft      reference at that same old physical time, not the next time
%   ~          unused next-time reference argument retained in the interface
%   m          model geometry, linear operator, forcing and observation mask
%   c          experiment settings from config.m and the selected suite
%   backward   whether this is the backward assimilation pass
%   truthOnly  generate a reference with the physical model and no feedback
% For NSE, z and uLeft are N-by-N arrays of streamfunction coefficients.
%
% The decomposition is z_t = L*z + N(z), with diagonal Fourier operator L.
% The linear part is exponentiated. N(z) includes nonlinear terms and forcing.
% Nudging uses old-time data separately from the physical RK stages.
% RK4 describes the physical stages, not fourth-order accuracy of the complete
% mixed method with its explicit observation corrector.
% See run.m for the pass loops and nonlinear.m for the actual equations.
% h<0 and gain<0 backward, so the feedback remains attractive.
if nargin<9,truthOnly=false;end
% Start from the forward physical model on every call. V=1 means no smoothing.
% For NSE, m.L=-nu*(kx^2+ky^2). Nothing here initially reverses that operator.
L=m.L;V=ones(m.shape);
% The backward gain is negative because h is negative. Their product h*gain
% is positive on BOTH passes, attracting the estimate toward observations.
% Reference generation sets gain=0 and bypasses backward regularization.
if truthOnly,gain=0;elseif backward,gain=-c.muBackward;else,gain=c.muForward;end
if backward && ~truthOnly
 % Only modify backward assimilation. "filtered" confines regularization to
 % unobserved coefficients, leaving the observed coefficients unregularized.
 % It does not change which coefficients are observed or delete any state.
 active=ones(m.shape);if c.backward.filtered,active=double(~m.obs);end
 switch c.backward.scheme
  % Viscosity-only truncation. m.viscousL=-nu*k^2, so subtracting it removes
  % viscosity at radius >= viscosityCutoff. Advection and state coefficients
  % remain. This is distinct from the input mask in nonlinear.m.
  case 'cutoff'
   L=L-m.viscousL.*(m.radius>=c.backward.viscosityCutoff);
  % For the default viscosityFactor=-1, L changes from -nu*k^2 to +nu*k^2.
  % Since h=-dt backward, exp(h*L) then damps rather than amplifies these modes.
  % The analogous dampingFactor acts on a physical scalar damping term, if any.
  % With filtering, these alterations apply only outside the observation mask.
  case 'diffusive'
   L=L+active.*((c.backward.viscosityFactor-1)*m.viscousL+ ...
    (c.backward.dampingFactor-1)*m.dampedL);
  % The backward Voigt operator is (I-alpha^2*Delta)^(-1), whose multiplier
  % is V=1/(1+alpha^2*k^2). Apply V to L here, to N below, and to feedback.
  % Smoothing all three terms implements the confirmed Voigt formulation.
  % V reduces the high-wavenumber backward growth rate without changing its sign.
  % Voigt is not the same operation as reversing viscosity to add diffusion.
  case 'voigt'
   V=1./(1+c.backward.alpha^2*m.k2.*active);L=V.*L;
 end
end
% E advances the diagonal linear equation over one signed step. H is its
% half-step factor for integrating-factor RK4. Array products are elementwise.
% For backward standard viscous BFN, E grows at high k. This is the source of
% the backward instability. Regularized variants alter the multiplier above.
E=exp(h*L);H=exp(h*L/2);
switch c.integrator
 % The paper's transport/Burgers/NSE driver uses integrating-factor Euler.
 % Evaluate nonlinear terms and forcing at OLD z, multiply them by h*V,
 % then advance this physical increment with E. This is first order overall.
 % For standard NSE the line below reads
 %   physical = exp(-h*nu*k^2) .* (z + h*N(z)).
 case 'ifeuler'
  physical=E.*(z+h*V.*paperbfn.nonlinear(z,m,c));
 % Four integrating-factor RK stages for the physical equation. a, b, d, e
 % are stage derivatives, not observations. With L=0 and V=1 these reduce
 % to ordinary RK4. No interpolation of observation data is used here.
 case {'rk4','ifrk4'}
  a=V.*paperbfn.nonlinear(z,m,c);
  b=V.*paperbfn.nonlinear(H.*(z+h*a/2),m,c);
  d=V.*paperbfn.nonlinear(H.*z+h*b/2,m,c);
  e=V.*paperbfn.nonlinear(E.*z+h*H.*d,m,c);
  physical=E.*z+h/6*(E.*a+2*H.*(b+d)+e);
  % Use a fixed Lorenz summation order for reproducible arithmetic.
  % Floating point reassociation can affect long parameter sweeps.
  if strcmp(c.model,'lorenz') && isfield(c.lorenz,'originalArithmetic') && c.lorenz.originalArithmetic
   physical=z+h/6*(a+2*b+2*d+e);
  end
end
% Feedback uses the current state's physical time. Paper Lorenz runs use
% shift=0 and half-open windows. The historical option selects an alternate
% discrete convention with a backward observation-time offset of -dt.
obsTime=t;
if strcmp(c.model,'lorenz') && h<0 && strcmp(c.observe.schedule,'thirds') && isfield(c.lorenz,'backwardObservationShift')
 obsTime=(round(t/c.dt)+c.lorenz.backwardObservationShift)*c.dt;
end
% obs is 0/1 on coefficients or components. The residual uLeft-z is formed
% before the step, so this is an explicit correction, not stagewise nudging.
% For NSE it is a streamfunction residual on the same observed Fourier modes.
obs=paperbfn.observed_mask(m,c,obsTime);
feedback=h*gain*V.*obs.*(uLeft-z);
% feedbackFactor places the explicit correction inside the integrating
% factor by multiplying it by E. It is false for the shared
% NSE/transport/Burgers defaults and true in the separate KdV implementation.
if c.feedbackFactor,feedback=E.*feedback;end
% Return the complete updated state. No post-step projection is performed.
% In particular, unobserved and out-of-band coefficients are retained.
next=physical+feedback;
end
