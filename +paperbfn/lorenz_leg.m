function [v,rows,lastIdx,failed]=lorenz_leg(v,h,truth,indices,m,c,cycle,leg)
%LORENZ_LEG Same RK4 and explicit feedback, with scalar work inside the loop.
% Time loop for ONE Lorenz forward or backward pass. Unlike PDEs there is
% no Fourier grid. v and each truth column are three real shifted coordinates.
% h sets the direction. indices lists the reference's old-time sample indices.
% rows uses the same eleven diagnostic columns as run.m.
%
% The physical stages implement ordinary RK4. Feedback is evaluated once
% using the old v and the old reference, after assembling the physical update.
% The separate implementation avoids repeated general model dispatch in the
% large Lorenz sweeps. It is the assimilation path actually called by run.m.
% Avoid rebuilding spectral arrays and diagnostic structs at every ODE step.
s=c.lorenz.sigma;r=c.lorenz.rho;b=c.lorenz.beta;
stride=c.output.diagnosticStride;steps=numel(indices);limit=c.stop.maxNorm;
% Backward h and gain both change sign. h*gain stays attractive on both passes.
gain=c.muForward;if h<0,gain=-c.muBackward;end
% Precompute a component-by-time feedback mask for this pass. Continuous
% observations repeat the fixed component mask. The windowed experiment instead
% observes each selected component during its own allocated third of the window.
mask=repmat(m.obs,1,steps);
if strcmp(c.observe.schedule,'thirds')
 % Subtract 1 because reference column 1 is physical t=0. Paper runs use
 % shift=0, matching the current state. The alternate historical schedule
 % uses shift=-1 to evaluate its backward mask at physical time t-dt.
 timeIndices=indices-1;
 if h<0 && isfield(c.lorenz,'backwardObservationShift'),timeIndices=timeIndices+c.lorenz.backwardObservationShift;end
 t=timeIndices*c.dt;
 mask(:)=false;
 for component=c.observe.components
  % A fraction f observes this component over [left,left+f*window/3).
  % The half-open right endpoint prevents overlap between successive thirds.
  % With f=0 this interval is empty, including at t=0. The legacyClosed option
  % instead includes right endpoints and can observe isolated endpoints.
  left=(component-1)*c.window/3;
  if isfield(c.lorenz,'windowRule') && strcmp(c.lorenz.windowRule,'legacyClosed')
   if component==1,leftTest=t>=left;else,leftTest=t>left;end
   mask(component,:)=leftTest & t<=left+c.observe.fractions(component)*c.window/3;
  else
   mask(component,:)=(t>=left & t<left+c.observe.fractions(component)*c.window/3);
  end
 end
end
rows=nan(ceil(steps/stride)+1,11);count=0;failed=false;
% Advance one signed physical step. lastIdx points to the NEW reference time.
% u is still the reference at the OLD time for explicit feedback.
for j=1:steps
 idx=indices(j);lastIdx=idx+sign(h);u=truth(:,idx);
 % RK4 stages at old time, the two half-step states, and the full-step state.
 % There is no feedback inside these stages, so observations need no RK-stage
 % interpolation. a,d,e,f correspond to the usual k1,k2,k3,k4 derivatives.
 a=rhs(v,s,r,b);d=rhs(v+h*a/2,s,r,b);
 e=rhs(v+h*d/2,s,r,b);f=rhs(v+h*e,s,r,b);
 physical=v+h/6*(a+2*(d+e)+f);
 % Choose the ungrouped RK4 sum when originalArithmetic is enabled.
 % The two sums agree algebraically, while roundoff depends on addition order.
 if isfield(c.lorenz,'originalArithmetic') && c.lorenz.originalArithmetic,physical=v+h/6*(a+2*d+2*e+f);end
 % Explicit corrector with old-time residual u-v. The new physical state is
 % not substituted into the residual, which would define a different method.
 v=physical+h*gain*mask(:,j).*(u-v);
 % Measure error at the NEW physical time. For diagnostic splitting use m.obs,
 % the full selected component set, rather than the time-dependent feedback mask.
 err=v-truth(:,lastIdx);power=abs(err).^2;
 total=sqrt(sum(power));bad=~isfinite(total) || total>limit;
 % Save sampled errors plus endpoints/failures. Each pass advances the plot's
 % cycle coordinate by 1/2, even though backward physical time decreases.
 if mod(j,stride)==0 || j==steps || bad
  count=count+1;obs=m.obs;
  rows(count,:)=[cycle-1+(leg-1)/2+j/(2*steps),(lastIdx-1)*c.dt, ...
   cycle,sign(h),j,total,sqrt(sum(power(obs))),sqrt(sum(power(~obs))),0,0,0];
 end
 if bad,failed=true;break;end
end
rows=rows(1:count,:);
end
% Shifted Lorenz equations. Let z=Z-(rho+sigma) in the standard Lorenz system.
% Then y_t=-sigma*x-y-x*z and z_t=x*y-beta*z-beta*(rho+sigma).
% This coordinate convention explains the constant term and the missing rho*x.
function q=rhs(v,s,r,b)
q=[s*(v(2)-v(1));-s*v(1)-v(2)-v(1)*v(3);-b*v(3)+v(1)*v(2)-b*(r+s)];
end
