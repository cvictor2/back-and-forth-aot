function verify_numerics
%VERIFY_NUMERICS Analytical properties independent of end-to-end plotting.
% Small numerical checks that run without recreating full manuscript figures.
% Run from the repository root with addpath('tests'); verify_numerics.
% Assertions test analytical identities and coefficient-level behavior, rather
% than whether a plot merely looks similar. A failing assertion stops execution.
% These checks are useful after arithmetic edits. They are not a proof of
% convergence or a full-resolution reproduction of every paper experiment.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
% On x=y=0, shifted Lorenz z satisfies a scalar affine linear ODE with a
% known exponential solution. Compare one physical RK4 step against that.
% truthOnly=true ensures no feedback enters this physical accuracy check.
c=paperbfn.config('lorenz');c=paperbfn.validate(c);m=paperbfn.model(c);
z=[0;0;1];u=[0;0;0];h=c.dt;
z1=paperbfn.step(z,h,0,u,u,m,c,false,true);
assert(norm(z1-[0;0;(1+c.lorenz.rho+c.lorenz.sigma)*exp(-c.lorenz.beta*h)-c.lorenz.rho-c.lorenz.sigma])<1e-13);
% Integrating-factor heat cycle leaves every unobserved coefficient unchanged.
% The high mode is unobserved. At this resolution its mode lies below the
% viscosity cutoff. Forward and backward exact linear factors therefore cancel.
% This checks signed-time conventions independently of a nonlinear model.
c=paperbfn.config('heat','high');c.N=128;c.nonlinearCutoff=21;c=paperbfn.validate(c);m=paperbfn.model(c);
z=m.u0;u=zeros(size(z));z1=paperbfn.step(z,h,0,u,u,m,c,false,false);
z2=paperbfn.step(z1,-h,h,u,u,m,c,true,false);
assert(norm(z2-z)<1e-10);
% Viscosity-only cutoff must retain out-of-band coefficients exactly.
% Put a nonzero Fourier coefficient explicitly above the viscosity cutoff.
% It must survive unchanged when that mode's diffusion multiplier is removed.
% An accidental projection of the state would make this assertion fail.
c.backward.scheme='cutoff';c.backward.viscosityCutoff=20;m=paperbfn.model(c);
z=zeros(m.shape);z(31)=1;
z1=paperbfn.step(z,-h,0,u,u,m,c,true,false);assert(z1(31)==1);
% Full Voigt smooths feedback; filtered Voigt restores observed-mode feedback.
% With zero viscosity and a single observed residual, the update isolates
% feedback itself. Full Voigt must divide it by 1+alpha^2*k^2. Filtered Voigt
% must leave that same observed residual unsmoothed.
c.backward.scheme='voigt';c.backward.alpha=.1;c.nu=0;m=paperbfn.model(c);
z=zeros(m.shape);u=z;u(2)=1;
a=paperbfn.step(z,-h,0,u,u,m,c,true,false);
assert(abs(a(2)-h*c.muBackward/(1+.1^2*pi^2))<1e-14);
c.backward.filtered=true;
b=paperbfn.step(z,-h,0,u,u,m,c,true,false);assert(abs(b(2)-h*c.muBackward)<1e-14);
% Parseval and orthogonal observed/unobserved split, without diagnostic floors.
z=fft(cos((0:127)'*2*pi/128)+.05*cos(30*(0:127)'*2*pi/128));
d=paperbfn.errors(z,zeros(128,1),m);
assert(abs(d(1)^2-1.0025)<1e-12);assert(abs(d(1)^2-d(2)^2-d(3)^2)<1e-12);
% NSE streamfunction derivative orientation and integral velocity norm.
% Use analytic psi=cos(x)*cos(y) without loading a restart. Its velocities
% are u=cos(x)*sin(y), v=-sin(x)*cos(y), and integral velocity L2=sqrt(2)*pi.
% These checks detect swapped axes, derivative signs and FFT normalization.
c=paperbfn.config('nse');c.N=64;c.nonlinearCutoff=64/3;
m=paperbfn.model(c,true);[x,y]=ndgrid((0:63)*2*pi/64);psi=fft2(cos(x).*cos(y));
u=-real(ifft2(m.dy.*psi));v=real(ifft2(m.dx.*psi));
assert(norm(u-cos(x).*sin(y),'fro')<1e-12);assert(norm(v+sin(x).*cos(y),'fro')<1e-12);
d=paperbfn.errors(psi,zeros(64),m);assert(abs(d(1)-sqrt(2)*pi)<1e-12);
% Independent NSE advection check with two Laplacian eigenvalues. A single
% Taylor-Green eigenfunction has zero advection and cannot detect its sign.
% For psi=cos(x)+.2*cos(2*y) and velocity=(-psi_y,psi_x), the inviscid,
% unforced nonlinear psi_t is -.24*sin(x)*sin(2*y).
m.invLap=zeros(m.shape);nz=m.k2>0;m.invLap(nz)=-1./m.k2(nz);
m.forcing=zeros(m.shape);
psi=fft2(cos(x)+.2*cos(2*y));q=paperbfn.nonlinear(psi,m,c);
expected=-.24*sin(x).*sin(2*y);
assert(norm(real(ifft2(q))-expected,'fro')<1e-12);
% The quadratic velocity products are invariant under simultaneous reversal
% of both components. Check this algebraic symmetry independently.
negU=real(ifft2(m.dy.*(psi.*m.keep)));negV=real(ifft2(-m.dx.*(psi.*m.keep)));
symmetricQ=m.invLap.*((m.kx.^2-m.ky.^2).*fft2(negU.*negV)+m.kx.*m.ky.*fft2(negV.^2-negU.^2));
assert(isequal(q,symmetricQ));
% Compare the Basdevant expression with a separately evaluated vorticity
% advection. omega=Delta*psi and psi_t=-Delta^{-1}(velocity dot grad omega).
u=-real(ifft2(m.dy.*psi));v=real(ifft2(m.dx.*psi));omega=-m.k2.*psi;
advection=u.*real(ifft2(m.dx.*omega))+v.*real(ifft2(m.dy.*omega));
independent=-m.invLap.*fft2(advection);
assert(norm(q-independent,'fro')/max(norm(independent,'fro'),1)<1e-12);
% Seven KdV variants on short resolved windows, including error decomposition.
variants={'BFN','damped BFN','diffusive BFN','Voigt BFN','BFN','diffusive BFN','Voigt BFN'};
% These very short KdV windows remain finite even for the variant that fails
% at the full paper window. Check that all seven update paths execute and
% that observed/unobserved squared L2 errors sum to the total at each sample.
for i=1:7
 model='damped';if i>=5,model='viscous';end
 o=struct('model',model,'variant',variants{i},'N',128,'T',.0001, ...
 'dt',1e-5,'cycles',2,'sampleEvery',1,'makePlot',false,'saveOutput',false);
 r=paperbfn.kdv_run(o);assert(~r.failed && height(r.endpoints)==2);
 assert(max(abs(r.history.L2.^2-r.history.observedL2.^2-r.history.unobservedL2.^2))<1e-12);
end
fprintf('Analytical solver checks and all seven KdV variants passed.\n');
end
