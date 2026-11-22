function m=model(c,geometryOnly)
%MODEL Fourier geometry. NSE evolves streamfunction, as in BFN_runner.m.
% Build the representation used by step.m, nonlinear.m and errors.m.
% c holds user settings. m holds arrays derived from those settings.
%
% Important fields
%   L            diagonal forward physical operator in Fourier space
%   viscousL     viscosity contribution -nu*k^2, stored separately for BFN
%   dampedL      scalar damping contribution -gamma
%   obs          static observation mask, not a mask of the evolving state
%   keep         mask used only on inputs to nonlinear products
%   u0           reference initial state before any configured spinup
%   forcing      forcing in the evolved variable's Fourier representation
%   shape/count  array dimensions and total number of coefficients
%   volume       physical domain length in 1D, area in 2D
%
% For NSE the state is streamfunction psi, with velocity (-psi_y,psi_x)
% and vorticity Delta*psi. It is not a velocity or vorticity array.
% geometryOnly=true builds masks/multipliers for norms and plots without
% loading the restart or constructing initial data. Early returns matter.
if nargin<2, geometryOnly=false; end
m.name=c.model;
% Lorenz has three real components and no diagonal linear Fourier operator.
% Its entire shifted right-hand side is evaluated in nonlinear.m / lorenz_leg.m.
if strcmp(c.model,'lorenz')
 m.shape=[3 1];m.count=3;m.keep=true(3,1);m.k2=zeros(3,1);
 m.obs=false(3,1);m.obs(c.observe.components)=true;
 m.u0=c.initial.reference(:);m.forcing=zeros(3,1);
 m.L=zeros(3,1);m.viscousL=m.L;m.dampedL=m.L;return
end
% MATLAB FFT order stores nonnegative frequencies first, then negative ones.
% index is an integer mode number. k=2*pi*index/L is a physical wavenumber.
% Cutoff settings use integer indices/radii, while derivatives use physical k.
n=c.N; index=[0:n/2-1 -n/2:-1]; k=2*pi/c.L*index;
% NSE lives on a periodic L-by-L square. ndgrid makes kx vary along rows
% and ky along columns. These orientations must agree with the saved restart.
% The source assigns a positive Nyquist entry in the Laplacian grid.
if strcmp(c.model,'nse')
 k(n/2+1)=pi*n/c.L;
 [m.kx,m.ky]=ndgrid(k,k); [ix,iy]=ndgrid(index,index);
 m.dx=1i*m.kx;m.dy=1i*m.ky;
 % Set the unpaired Nyquist first derivative to zero to preserve the real-field
 % convention of the source. The Laplacian still uses its nonzero k^2.
 % Do not replace the Laplacian by abs(dx)^2+abs(dy)^2 at these entries.
 m.dx(n/2+1,:)=0;m.dy(:,n/2+1)=0;
 m.k2=m.kx.^2+m.ky.^2;m.radius=sqrt(ix.^2+iy.^2);
 % Same index mask as Parameters.m, including its endpoint convention.
 % Rectangular nonlinear-input mask, approximately a per-coordinate N/3 cutoff.
 % MATLAB stores the excluded interval between the low positive and negative
 % frequencies. The explicit endpoints below define the retained rectangle.
 % The mean is excluded from nonlinear inputs. Nothing zeros the stored state.
 % This mask is NOT the circular observed-mode region defined below.
 keep=true(n);lo=ceil(n/3);hi=n-floor(n/3)+1;
 keep(lo:hi,:)=false;keep(:,lo:hi)=false;keep(1)=false;
 m.keep=keep;m.shape=[n n];m.count=n*n;m.volume=c.L^2;
else
 % For 1D models, store Fourier states as N-by-1 columns. Retain nonzero
 % modes up to nonlinearCutoff only when evaluating nonlinear products.
 m.k=k(:);m.k2=m.k.^2;m.radius=abs(index(:));
 m.keep=m.radius<=c.nonlinearCutoff & m.radius>0;
 m.shape=[n 1];m.count=n;m.volume=c.L;
end
% Observed modes obey the STRICT inequality 0 < radius < cutoff.
% For NSE radius=sqrt(index_x^2+index_y^2), giving a circular observed set.
% "full" includes every nonzero Fourier mode. The mean is still excluded.
m.obs=m.radius<c.observe.cutoff & m.radius>0;
if c.observe.full,m.obs=true(m.shape);m.obs(m.radius==0)=false;end
if geometryOnly,return;end
% Construct the forward physical operator. Keeping the viscous and damping
% parts separate lets step.m alter only the requested backward contribution.
% NSE uses L=-nu*(kx^2+ky^2), and Burgers/heat use L=-nu*k^2.
m.viscousL=-c.nu*m.k2;m.dampedL=-c.gamma*ones(m.shape);
m.L=m.viscousL+m.dampedL;m.forcing=zeros(m.shape);
% Transport evolves +a*u_x, giving +i*a*k in the linear multiplier.
% With positive a, profiles translate in the negative x direction.
if strcmp(c.model,'transport'),m.L=m.L+1i*c.transportSpeed*m.k;end
if strcmp(c.model,'nse')
 % Load the paper's NSE target and forcing rather than generate a new random
 % flow. The restart must have the requested grid resolution. suite_nse.m
 % explicitly constructs a smaller, normalized restart for smoke mode.
 a=load(c.initial.restartFile,'psi_hat_initial','f_hat');
 assert(isequal(size(a.psi_hat_initial),m.shape),'Restart resolution must match N.');
 m.u0=a.psi_hat_initial;
 % Invert Delta only away from its zero mode. Setting invLap(0)=0 defines the
 % undefined mean of streamfunction and avoids division by zero.
 % The saved f_hat is curl(f)=d_x f_2-d_y f_1. invLap*f_hat is the
 % forcing that belongs in the streamfunction evolution equation.
 m.invLap=zeros(m.shape);nz=m.k2>0;m.invLap(nz)=-1./m.k2(nz);
 invLap=zeros(m.shape);nz=m.k2>0;invLap(nz)=-1./m.k2(nz);
 m.forcing=invLap.*a.f_hat;
else
 % Periodic 1D grid excludes the duplicated right endpoint. The low component
 % has integer mode 1 and amplitude 1. The high component has the configured
 % mode/amplitude. high, low and mixed select the three observation examples.
 % FFT produces MATLAB's unnormalized coefficients. Norms account for that.
 x=(0:n-1)'*c.L/n;
 low=cos(2*pi*x/c.L);high=c.initial.highAmplitude*cos(c.initial.highMode*2*pi*x/c.L);
 switch c.initial.kind
  case 'low',u=low;
  case 'high',u=high;
  case 'mixed',u=low+high;
  otherwise,error('Unknown profile.');
 end
 m.u0=fft(u);
end
end
