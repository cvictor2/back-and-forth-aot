function d = errors(v,u,m)
%ERRORS Absolute integral L2 and H1-seminorm errors; no display floors.
% Input v is the estimated state and u is the reference at the SAME time.
% Return six absolute errors in this order
%   [L2 total, L2 observed, L2 unobserved,
%    H1 seminorm total, H1 seminorm observed, H1 seminorm unobserved].
% Lorenz uses Euclidean norms and returns zeros in the H1 slots.
% NSE input states are streamfunction coefficients, but reported norms are
% VELOCITY norms. In 1D the norms are those of the evolved scalar field.
%
% The H1 quantity is the gradient seminorm, not sqrt(L2^2+gradient^2).
% No relative normalization, coefficient clipping or plotting floor is used.
% Split the error with the fixed observation mask. In windowed Lorenz runs,
% this labels the selected component set, not its instantaneous active window.
% Observed and complementary squared errors add to the total squared error.
e=v-u; obs=m.obs;
if strcmp(m.name,'lorenz')
    q=abs(e).^2; hq=zeros(size(q)); scale=1;
else
    % MATLAB's forward FFT is unnormalized. Parseval gives integral L2 squared
    % as volume/count^2 times the sum of coefficient powers. Taking the square
    % root yields this scale. For NSE it is L/N^2 on the L-by-L periodic square.
    scale=sqrt(m.volume)/m.count;
    % Velocity Fourier coefficients are (-dy*psi,dx*psi). Sum their powers before
    % forming norms. This respects the zero Nyquist first-derivative convention.
    % Do not take the norm of psi itself and call it a velocity error.
    if strcmp(m.name,'nse')
        % Streamfunction coefficients give velocity power |k|^2 |psi|^2.
        q=zeros(size(e)); k=m.k2>0; q=(abs(m.dx).^2+abs(m.dy).^2).*abs(e).^2;
    else
        q=abs(e).^2;
    end
    % One further spatial derivative contributes k^2 to the squared norm.
    % For NSE this measures velocity gradients rather than streamfunction gradients.
    hq=m.k2.*q;
end
% Flatten all spatial dimensions so summation covers every coefficient.
% The observation split is orthogonal in Fourier space, so no cross terms enter.
q=q(:);hq=hq(:);obs=obs(:);
d=scale*[sqrt(sum(q)),sqrt(sum(q(obs))),sqrt(sum(q(~obs))), ...
    sqrt(sum(hq)),sqrt(sum(hq(obs))),sqrt(sum(hq(~obs)))];
end
