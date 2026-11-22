function result = kdv_run(options)
% Coefficient tests using the Fourier and IF-RK4 formulas in BFN/KdV.m.
% Self-contained KdV experiment driver, separate from the shared PDE step.m.
% Read this file from configuration, through Fourier setup and reference,
% to the forward/backward loops. physicalStep and rhs at the bottom contain
% the numerical method. All arrays here use row-vector Fourier states.
%
% Forward physical equations represented by this driver
%   damped   u_t + u*u_x + dispersion^2*u_xxx + damping*u = forcing
%   viscous  u_t + u*u_x + dispersion^2*u_xxx = viscosity*u_xx + forcing
% Both use integrating-factor RK4 for the physical terms and an explicit
% old-time nudging correction. Observations are not evaluated at RK stages.
% RK4 here describes the physical substep rather than a fourth-order claim
% for the complete method after the explicit nudging correction is included.
%
% The backward variant is chosen independently of the forward physical model.
% Defaults specify the paper experiments. In particular the damped
% diffusive panel RETAINS damping=.5 and adds stabilization=1e-5. Its damping
% is not canceled unless cancelDamping is explicitly requested.
% result = kdv_run(struct('model','damped','damping',0.5));
% Physical model and backward regularization are independent choices.
if nargin < 1, options = struct(); end
% All user-adjustable settings are collected here. T and dt are physical
% time, cycles counts forward/backward pairs, N is spatial points, Lx is
% periodic length. alpha is squared inside the Voigt multiplier. The initial
% estimate defaults to the exact reference to test backward stability itself.
cfg = struct('model','damped','variant','BFN','filtered',false, ...
    'N',512,'Lx',2,'dt',1e-5,'T',1,'cycles',5, ...
    'cancelDamping',false,'damping',0.5,'viscosity',1e-3,'stabilization',[], ...
    'alpha',[],'dispersion',1,'mu',1000,'observationCutoff',15, ...
    'lowAmplitude',1,'highAmplitude',0.1,'highMode',20, ...
    'forcingAmplitude',1,'initialization','exact', ...
    'diagnostic','integral','sampleEvery',100,'maxAmplitude',1e100, ...
    'makePlot',true,'saveOutput',true,'outputRoot', ...
    fullfile(fileparts(fileparts(mfilename('fullpath'))),'results','kdv'));
assert(isstruct(options) && isscalar(options),'Options must be a scalar struct.');
% Apply only recognized overrides. Reject misspelled parameters rather than
% silently running a default experiment that differs from the intended one.
f = fieldnames(options);
for j=1:numel(f)
    assert(isfield(cfg,f{j}), 'Unknown option %s.',f{j});
    cfg.(f{j}) = options.(f{j});
end
cfg.model=validatestring(cfg.model,{'damped','viscous'});
% The paper uses different Voigt parameters for damped and viscous KdV.
% Both defaults are explicit so selecting a physical model sets its alpha.
if isempty(cfg.alpha)
 if strcmp(cfg.model,'damped'),cfg.alpha=1e-7;else,cfg.alpha=1e-2;end
end
assert(isscalar(cfg.cancelDamping) && ismember(cfg.cancelDamping,[false true]));
cfg.variant=validatestring(cfg.variant,{'BFN','damped BFN','diffusive BFN','Voigt BFN'});
cfg.initialization=validatestring(cfg.initialization,{'exact','zero','observed'});
cfg.diagnostic=validatestring(cfg.diagnostic,{'legacy','integral'});
validateattributes(cfg.N,{'numeric'},{'scalar','integer','>=',8});
assert(mod(cfg.N,2)==0,'N must be even.');
for name={'Lx','dt','T','maxAmplitude'}
    validateattributes(cfg.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
for name={'damping','viscosity','alpha','mu'}
    validateattributes(cfg.(name{1}),{'numeric'},{'scalar','real','finite','nonnegative'});
end
for name={'cycles','sampleEvery'}
    validateattributes(cfg.(name{1}),{'numeric'},{'scalar','integer','positive'});
end
assert(cfg.observationCutoff>1 && cfg.observationCutoff<cfg.N/2, ...
    'Observation cutoff must lie between 1 and N/2.');
assert(cfg.highMode<cfg.N/2 && mod(cfg.highMode,1)==0,'highMode must be a resolved integer.');
% The reference and both passes must land on the same discrete time grid.
% There is no interpolation of reference values between these grid points.
steps=round(cfg.T/cfg.dt);
assert(steps>=1 && abs(steps*cfg.dt-cfg.T)<=1e-10*cfg.T, ...
    'T must be an integer multiple of dt.');
% For the viscous model, a correction equal to viscosity cancels the physical
% anti-diffusion on the backward pass. This differs from the shared NSE dBFN
% choice, which reverses the viscous sign with a correction of twice viscosity.
% For the damped model, the additional diffusion is 1e-5.
if isempty(cfg.stabilization)
    if strcmp(cfg.model,'viscous'), cfg.stabilization=cfg.viscosity;
    else, cfg.stabilization=1e-5; end
end
validateattributes(cfg.stabilization,{'numeric'},{'scalar','finite','real','nonnegative'});
assert(~strcmp(cfg.variant,'damped BFN') || strcmp(cfg.model,'damped'), ...
    'Use damped BFN with the damped model. Use diffusive BFN for viscosity.');
% Periodic grid [0,Lx) excludes the duplicated endpoint. Unlike the shared
% 1D model, this source sets the Nyquist derivative multiplier to zero.
% m is integer FFT mode order, and k=2*pi*m/Lx is physical wavenumber.
N=cfg.N; dx=cfg.Lx/N; x=(0:N-1)*dx;
m=[0:N/2-1 0 -N/2+1:-1]; k=(2*pi/cfg.Lx)*m;
% mask filters inputs to the nonlinear product only. obs selects the strict
% nonzero observation band. Neither is applied as a post-step state projection.
% Explicitly exclude Nyquist despite its zero entry in m for this source.
mask=(abs(m)<=floor(N/3)) & (m~=0); mask(N/2+1)=false;
obs=(abs(m)<cfg.observationCutoff) & (m~=0); obs(N/2+1)=false;
% Initial reference is the prescribed low cosine plus a high cosine. Fourier
% coefficients are unnormalized MATLAB fft values. The forcing is exp(cos)
% with its zero Fourier mode removed.
u0=cfg.lowAmplitude*cos(2*pi*x/cfg.Lx)+ ...
    cfg.highAmplitude*cos(cfg.highMode*2*pi*x/cfg.Lx);
u0hat=fft(u0);
fhat=fft(cfg.forcingAmplitude*exp(cos(2*pi*x/cfg.Lx))); fhat(1)=0;
% Separate the REAL diagonal physical multiplier from imaginary dispersion.
% Scalar damping contributes -damping. Viscosity contributes -viscosity*k^2.
% These physical coefficients do not change between forward reference steps.
if strcmp(cfg.model,'damped'), realForward=-cfg.damping*ones(1,N);
else, realForward=-cfg.viscosity*k.^2; end
% Standard BFN integrates the same physical RHS with negative dt. Therefore
% negative realForward becomes amplification on a backward pass. The following
% variants alter the backward RHS before that negative step is applied.
realBackward=realForward;
switch cfg.variant
    % Reversing the damping sign makes backward negative-dt evolution dissipative.
    % This is equivalent to adding twice the forward scalar damping to the RHS.
    case 'damped BFN'
        realBackward=+cfg.damping*ones(1,N);
    % Add +stabilization*k^2 to the signed-time RHS. With h=-dt this acts as
    % diffusion, even though the multiplier appears with a positive sign here.
    % The optional cancelDamping affects only the damped physical model.
    case 'diffusive BFN'
        if strcmp(cfg.model,'damped')
            % Retain the forward scalar damping and add the diffusion correction.
            realBackward=realForward+cfg.stabilization*k.^2;
            if cfg.cancelDamping,realBackward=cfg.stabilization*k.^2;end
        else
            realBackward=realForward+cfg.stabilization*k.^2;
        end
end
% Voigt smooths the whole backward RHS through V, including the nonlinear
% stages, linear operator and feedback. Filtering restores V=1 on obs.
% The mean and Nyquist multiplier remain 1 under the source convention.
V=ones(1,N);
if strcmp(cfg.variant,'Voigt BFN')
    V=1./(1+cfg.alpha^2*k.^2);
    if cfg.filtered, V(obs)=1; end
    V(1)=1; V(N/2+1)=1;
end
% For filtered diffusion/damping corrections, retain the forward physical
% real multiplier on observed modes. Unobserved modes keep the correction.
if cfg.filtered, realBackward(obs)=realForward(obs); end
% The Fourier multiplier of -dispersion^2*u_xxx is +i*dispersion^2*k^3.
% Combine it with the real physical/regularized part. Backward Voigt multiplies
% both real and imaginary components by V, not just the nonlinear term.
Lforward=realForward+1i*cfg.dispersion^2*k.^3;
Lbackward=V.*(realBackward+1i*cfg.dispersion^2*k.^3);
% Store the reference exactly as KdV.m does. No state projection is applied.
% Reference storage is time-by-mode in THIS driver, unlike the column-per-time
% layout in reference.m. Each row is one full Fourier state. Store one physical
% window, shared by all cycles, with no nudging or state projection.
truth=complex(zeros(steps+1,N)); truth(1,:)=u0hat;
% Precompute integrating factors for the forward reference. physicalStep uses
% these same physical stages during forward assimilation.
H=exp(cfg.dt/2*Lforward); E=exp(cfg.dt*Lforward);
for j=1:steps
    truth(j+1,:)=physicalStep(truth(j,:),cfg.dt,H,E,ones(1,N));
    assert(all(isfinite(truth(j+1,:))), 'Reference became nonfinite at step %d.',j);
end
% "exact" starts the estimate from u0hat. "observed" projects only the initial
% guess. Even exact initialization need not survive unstable backward evolution.
switch cfg.initialization
    case 'exact', v=u0hat;
    case 'zero', v=zeros(1,N);
    case 'observed', v=u0hat.*obs;
end
% Save sampled physical-time errors and completed-cycle initial-state errors
% separately. Endpoint rows are initially NaN, so failed/incomplete cycles
% cannot be mistaken for exact zero error. The stored state rows follow suit.
capacity=1+2*cfg.cycles*(ceil(steps/cfg.sampleEvery)+1);
history=nan(capacity,12); nh=0;
endpoints=nan(cfg.cycles,8); endpointStates=complex(nan(cfg.cycles,N));
failed=false; failureMessage=''; failureCycleCoordinate=NaN;
lastPhysicalTime=0; lastCycleCoordinate=0;
record(0,0,0,u0hat,v);
% Odd pass numbers move 0->T; even pass numbers move T->0. The estimate v is
% carried from pass to pass. A pair of passes advances plot coordinate by 1.
for leg=1:2*cfg.cycles
    forward=mod(leg,2)==1;
    if forward
        h=cfg.dt; signedMu=cfg.mu; L=Lforward; Vk=ones(1,N); idx=1;
    else
        h=-cfg.dt; signedMu=-cfg.mu; L=Lbackward; Vk=V; idx=steps+1;
    end
    % Use the selected signed step and backward multiplier for this pass. Negative
    % h does not change the observation attraction, because signedMu also flips.
    H=exp(h/2*L); E=exp(h*L);
    for j=1:steps
        % Keep the OLD state for the explicit residual and for failure recovery.
        % physicalStep advances only the physical/regularized equation. Its stages
        % never evaluate observations. Feedback is then added as a separate corrector.
        old=v;
        physical=physicalStep(old,h,H,E,Vk);
        % Old-time observations and explicit feedback, as in KdV.m.
        % This scheme weights the explicit feedback by E as well as V.
        % This E weighting is part of this KdV scheme and differs from shared NSE
        % feedbackFactor=false. Residual truth(idx,:)-old uses the old physical time.
        v=physical+h*signedMu*Vk.*E.*obs.*(truth(idx,:)-old);
        % Now advance the reference index to the NEW physical time for diagnostics.
        % The coordinate is in cycles, so each complete pass occupies half a unit.
        if forward, idx=idx+1; else, idx=idx-1; end
        physicalTime=(idx-1)*cfg.dt;
        coordinate=((leg-1)+j/steps)/2;
        % Stop a blow-up and restore the last valid Fourier state rather than store
        % a fabricated endpoint. Failure time and the last valid time are kept apart.
        % At paper resolution the viscous standard-BFN panel is expected to fail.
        if ~all(isfinite(v)) || max(abs(ifft(v)))>cfg.maxAmplitude
            failed=true; failureCycleCoordinate=coordinate;
            failureMessage=sprintf('State exceeded the finite/amplitude limit in pass %d, step %d.',leg,j);
            v=old; break
        end
        lastPhysicalTime=physicalTime; lastCycleCoordinate=coordinate;
        if mod(j,cfg.sampleEvery)==0 || j==steps
            record(coordinate,physicalTime,leg,truth(idx,:),v);
        end
    end
    if failed, break; end
    % A successful backward pass returns to physical t=0. Compare against u0hat,
    % not against the terminal reference at T. These rows populate endpoint tables.
    if ~forward
        c=leg/2; d=diagnostics(v-u0hat);
        endpoints(c,:)=[c d(1:6) max(abs(ifft(v)))]; endpointStates(c,:)=v;
    end
end
% History stores six integral norms and three source-style scaled/floored
% norms for comparison only. New paper L2 plots use the integral diagnostics.
% The H1 columns are available numerically, but KdV figures show L2 only.
names={'cycleCoordinate','physicalTime','pass','L2','observedL2','unobservedL2', ...
    'H1seminorm','observedH1','unobservedH1','legacyL2','legacyObservedL2','legacyUnobservedL2'};
result=struct('config',cfg,'history',array2table(history(1:nh,:),'VariableNames',names), ...
    'endpoints',array2table(endpoints,'VariableNames', ...
    {'cycle','L2','observedL2','unobservedL2','H1seminorm','observedH1','unobservedH1','maxAmplitude'}), ...
    'endpointStates',endpointStates,'initialReferenceHat',u0hat, ...
    'terminalReferenceHat',truth(end,:),'lastEstimateHat',v, ...
    'x',x,'k',k,'observedMask',obs,'nonlinearMask',mask, ...
    'failed',failed,'failureMessage',failureMessage, ...
    'failureCycleCoordinate',failureCycleCoordinate, ...
    'lastPhysicalTime',lastPhysicalTime,'lastCycleCoordinate',lastCycleCoordinate, ...
    'outputFolder','');
result.sourceManifest=paperbfn.source_manifest();
result.matlabVersion=version;
fig=[];
% Plot copies of diagnostics. Zeros/nonpositive values cannot appear on a log
% axis, so hide them in the displayed copy without changing saved errors.
if cfg.makePlot
    fig=figure('Color','w','Visible','off','Position',[100 100 1100 650]);
    if strcmp(cfg.diagnostic,'legacy')
        % Same constant scaling and coefficient floor as the source plot.
        y=history(1:nh,4:6)*sqrt(cfg.Lx); y(:,1)=history(1:nh,10);
        % Obtain observed/unobserved floored diagnostics stored separately below.
        y(:,2:3)=history(1:nh,11:12);
    else, y=history(1:nh,4:6); end
    y(y<=0)=NaN;
    curves=semilogy(history(1:nh,1),y,'LineWidth',1.5);
    curves(1).LineWidth=1.8;curves(3).LineStyle=':';grid on;
    xlabel('Completed-cycle coordinate'); ylabel('L^2 error');xlim([0 cfg.cycles]);
    legend('Total','Observed','Unobserved','Location','best');
    title(sprintf('%s KdV, %s',cfg.model,cfg.variant));
    set(gca,'FontSize',14);
    if failed, subtitle(failureMessage,'Interpreter','none'); end
end
% Write raw states, settings, histories and endpoints beside editable FIG,
% PNG and vector PDF exports. Numerical data are retained even for failed runs.
if cfg.saveOutput
    if ~isfolder(cfg.outputRoot), mkdir(cfg.outputRoot); end
    folder=tempname(cfg.outputRoot); mkdir(folder); result.outputFolder=folder;
    save(fullfile(folder,'result.mat'),'result','-v7.3');
    writetable(result.history,fullfile(folder,'history.csv'));
    writetable(result.endpoints,fullfile(folder,'endpoints.csv'));
    fid=fopen(fullfile(folder,'settings.json'),'w');
    assert(fid>=0,'Cannot write settings.');
    fprintf(fid,'%s',jsonencode(cfg,PrettyPrint=true)); fclose(fid);
    if ~isempty(fig)
        savefig(fig,fullfile(folder,'error.fig'));
        exportgraphics(fig,fullfile(folder,'error.png'),'Resolution',250);
        exportgraphics(fig,fullfile(folder,'error.pdf'),'ContentType','vector');
        close(fig);
    end
    fprintf('Results saved in %s\n',folder);
end
if failed, warning('KdV:Stopped','%s',failureMessage); end

    % Integrating-factor RK4 for y_t=L*y+Vk*N(y). H=exp(h*L/2), E=exp(h*L).
    % Stages a,b,c,d evaluate the physical nonlinear-plus-forcing RHS at the
    % transformed intermediate states. Vk=1 forward and is Voigt smoothing when
    % selected backward. This helper contains NO nudging term.
    function z=physicalStep(y,h,H,E,Vk)
        a=Vk.*rhs(y);
        b=Vk.*rhs(H.*(y+h*a/2));
        c=Vk.*rhs(H.*y+h*b/2);
        d=Vk.*rhs(E.*y+h*H.*c);
        z=E.*y+h/6*(E.*a+2*H.*(b+c)+d);
    end
    % Use u*u_x=(1/2)*d_x(u^2). Mask only the input y, transform to real physical
    % space, square, FFT back, then differentiate. The coefficient is -i*k/2.
    % Add the prescribed mean-free forcing. Linear dispersion/dissipation live in L.
    function z=rhs(y)
        z=-0.5i*k.*fft(real(ifft(y.*mask)).^2)+fhat;
    end
    % Parseval normalization gives sqrt(Lx)/N for the integral L2 norm. H1
    % seminorm additionally multiplies error coefficients by k. Split obs and
    % its full complement without clipping coefficients to a display threshold.
    function d=diagnostics(e)
        f=sqrt(cfg.Lx)/N;
        d=[f*norm(e),f*norm(e(obs)),f*norm(e(~obs)), ...
            f*norm(k.*e),f*norm(k(obs).*e(obs)),f*norm(k(~obs).*e(~obs))];
    end
    % Record errors against the matching physical-time reference. Store the
    % coefficient-scaled Lx/N diagnostic and its 1e-14 floor in separate columns
    % so they cannot
    % be mistaken for the primary integral diagnostics used in the new tables.
    function record(coordinate,time,pass,u,vv)
        nh=nh+1; e=vv-u; d=diagnostics(e); floored=max(abs(e),1e-14);
        history(nh,:)=[coordinate,time,pass,d,dx*norm(floored), ...
            dx*norm(floored(obs)),dx*norm(floored(~obs))];
    end
end
