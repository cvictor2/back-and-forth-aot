function [truth,info,m] = reference(c,m)
%REFERENCE Store exactly ONE assimilation window, independent of cycle budget.
% Generate the target trajectory with the FORWARD physical model and no nudging.
% This is a discretely generated reference, not an analytic exact solution.
% Its time integrator is the same physical step used in forward assimilation.
%
% truth has m.count rows and stepsPerPass+1 columns. Column j+1 is the complete
% state at t=j*dt. A 2D NSE state is flattened columnwise for storage.
% info records the reference identity. m.u0 is updated to the post-spinup state.
% The allocation is for one window only, regardless of the number of BFN cycles.
% Describe everything that determines this reference. Cache keys include both
% settings and actual initial/forcing arrays, so a changed restart is detected.
key=struct('version',3,'model',c.model,'N',c.N,'L',c.L,'dt',c.dt, ...
    'window',c.window,'nu',c.nu,'gamma',c.gamma,'speed',c.transportSpeed, ...
    'dispersion',c.dispersion,'integrator',c.integrator, ...
    'nonlinearCutoff',c.nonlinearCutoff,'initial',c.initial,'forcing',c.forcing,'lorenz',c.lorenz);
% Hash actual initial/forcing arrays as well: restart contents may change at a path.
% Hash source text too. Even comment edits change this conservative identity
% and require a matching reference-cache identity.
key.solverHash=digest(uint8([fileread(which('paperbfn.step')) fileread(which('paperbfn.nonlinear')) fileread(which('paperbfn.model')) fileread(which('paperbfn.reference'))]));
key.initialHash=digest(m.u0); key.forcingHash=digest(m.forcing);
textKey=jsonencode(key); id=digest(uint8(textKey));
info=struct('key',textKey,'id',id,'steps',c.stepsPerPass,'dt',c.dt,'window',c.window);
file=fullfile(c.reference.root,[id '.mat']);
% A complex double requires 16 bytes. Check the reference budget before making
% the large array. Full NSE stores N^2 coefficients at every reference time.
bytes=16*m.count*(c.stepsPerPass+1);
assert(bytes<=c.reference.maxBytes,'bfn:Memory', ...
    'Truth history needs %.2f GB; reduce N/window, or explicitly raise reference.maxBytes.',bytes/1e9);
% Caching is optional and disabled by default. A cache must match both its
% settings key and array shape. It never replaces a differently configured run.
if c.reference.useCache && isfile(file)
    a=load(file,'truth','info');
    assert(strcmp(a.info.key,textKey),'bfn:Cache','Reference key mismatch.');
    assert(isequal(size(a.truth),[m.count,c.stepsPerPass+1]),'bfn:Cache','Reference shape mismatch.');
    truth=a.truth; info=a.info; m.u0=reshape(truth(:,1),m.shape); return
end
% Spin up only the physical reference. truthOnly=true disables all feedback
% and backward regularization. The configured PDE initial profile can therefore
% be evolved briefly before the assimilation window starts at physical t=0.
u=m.u0;
for j=1:round(c.initial.spinup/c.dt)
    u=paperbfn.step(u,c.dt,(j-1)*c.dt,u,u,m,c,false,true);
    assert(all(isfinite(u(:))),'bfn:Truth','Nonfinite spin-up state.');
end
% Define the reference initial state AFTER spinup. This is the target against
% which backward endpoint errors will be measured. Preserve every coefficient.
m.u0=u; truth=complex(zeros(m.count,c.stepsPerPass+1)); truth(:,1)=u(:);
% March the physical model forward and save every time level. Old-time
% observations in the BFN loops can then be accessed exactly by index, without
% interpolation. j=stepsPerPass ends at the physical observation window.
for j=1:c.stepsPerPass
    u=paperbfn.step(u,c.dt,(j-1)*c.dt,u,u,m,c,false,true);
    assert(all(isfinite(u(:))) && norm(u(:))<1e100,'bfn:Truth','Reference diverged.');
    truth(:,j+1)=u(:);
end
if c.reference.useCache
    if ~isfolder(c.reference.root), mkdir(c.reference.root); end
    save(file,'truth','info','-v7.3');
end
end
% Local SHA-256 helper. Hash real and imaginary parts as bytes, not rounded
% printed numbers, so the reference identity tracks the actual numerical data.
function s=digest(x)
if ~isa(x,'uint8')
    x=[typecast(real(x(:)),'uint8');typecast(imag(x(:)),'uint8')];
end
md=java.security.MessageDigest.getInstance('SHA-256'); md.update(x(:));
s=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
