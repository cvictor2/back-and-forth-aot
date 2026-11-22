function q=nonlinear(z,m,c)
%NONLINEAR Mask inputs only. Do not project the evolving state.
% Right-hand side left after the diagonal linear operator has been separated.
% For PDEs, input z and output q are Fourier coefficients of the same state.
% The name "nonlinear" also includes the additive forcing in m.forcing.
% Viscosity, transport and scalar damping live in m.L, not in this function.
% NSE evolves streamfunction psi, with velocity (u,v)=(-psi_y,psi_x)
% and vorticity Delta*psi.
% Names are local to each function. In run.m, v denotes the estimated state.
% Here the NSE local variables u and v denote the two physical velocities.
% These are the shifted Lorenz variables used in the paper. If (X,Y,Z) are
% standard Lorenz coordinates, use x=X, y=Y, z=Z-(rho+sigma). The shift
% explains both -sigma*x in the second equation and -b*(rho+sigma) in the third.
if strcmp(m.name,'lorenz')
 s=c.lorenz.sigma;r=c.lorenz.rho;b=c.lorenz.beta;
 q=[s*(z(2)-z(1));-s*z(1)-z(2)-z(1)*z(3); ...
  -b*z(3)+z(1)*z(2)-b*(r+s)];return
end
% Mask a LOCAL input copy before forming the nonlinear products.
% The caller's evolving state is not projected.
% m.keep is a 1D wavenumber mask or an NSE rectangular mask.
% This input filtering is separate from both observation and viscosity cutoffs.
z=z.*m.keep;
switch m.name
 % Heat and transport have no nonlinear terms. Their physical evolution is
 % already represented by the linear Fourier multiplier in model.m.
 case {'heat','transport'},q=m.forcing;
 % Use u*u_x = (1/2)*d_x(u^2). Transform the masked state to real space,
 % square there, transform back, and differentiate with i*k.
 % The minus sign gives u_t = -u*u_x + nu*u_xx when combined with m.L.
 case 'burgers'
  q=-0.5i*m.k.*fft(real(ifft(z)).^2)+m.forcing;
 case 'nse'
  % Basdevant streamfunction RHS for velocity (-psi_y,psi_x).
  % NSE velocity from the filtered streamfunction. dx=i*kx and dy=i*ky, so
  % u=-psi_y and v=psi_x. ifft2 supplies MATLAB's 1/N^2 normalization.
  % real removes roundoff-sized imaginary parts before physical-space products.
  u=-real(ifft2(m.dy.*z));v=real(ifft2(m.dx.*z));
  % Basdevant form evaluates quadratic velocity products without introducing a
  % separate pressure solve. Form u*v and v^2-u^2 in physical space, FFT them,
  % then combine their derivative multipliers and apply the inverse Laplacian.
  % invLap=-1/(kx^2+ky^2) off the zero mode, and 0 at the zero mode.
  % The saved f_hat is curl(f)=d_x f_2-d_y f_1. model.m applies invLap.
  % Thus psi_t=-invLap*(velocity dot grad(vorticity))+nu*Delta*psi
  %            +invLap*curl(f), with vorticity=Delta*psi.
  % The returned q is psi_t's nonlinear-plus-forcing part, NOT velocity_t.
  invLap=m.invLap;
  q=invLap.*((m.kx.^2-m.ky.^2).*fft2(u.*v)+ ...
    m.kx.*m.ky.*fft2(v.^2-u.^2))+m.forcing;
end
end
