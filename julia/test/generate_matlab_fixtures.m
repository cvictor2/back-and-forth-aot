function generate_matlab_fixtures
% Regenerate the small Julia parity fixtures using the repository MATLAB code.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);
dest=fullfile(root,'julia','test','fixtures');folder=tempname;mkdir(folder);
cases={};
for model={'transport','burgers','heat'}
 for scheme={'standard','cutoff','voigt','diffusive'}
  c=paperbfn.config(model{1});c.N=128;c.nonlinearCutoff=42;c.dt=1e-5;c.window=.0002;c.initial.spinup=.0001;c.cycles=2;c.output.diagnosticStride=2;
  c.nu=.05;c.backward.scheme=scheme{1};c.output.plot=false;c.output.root=folder;
  c.name=[model{1} '_' scheme{1}];r=paperbfn.run(c);j=convert_config(c);
  cases{end+1}=struct('name',c.name,'config_json',jsonencode(j),'initial',r.initialReference,'terminal',terminal(c),'final',r.finalState,'errors',r.diagnostics(:,6:11),'endpoints',r.cycleErrors);
 end
end
for fraction=[-1 0 .5 1]
 c=paperbfn.config('lorenz');c.window=.006;c.cycles=2;c.output.diagnosticStride=10;c.output.plot=false;c.output.root=folder;c.name=sprintf('lorenz_%g',fraction);
 if fraction>=0,c.observe.components=1:3;c.observe.schedule='thirds';c.observe.fractions=repmat(fraction,1,3);c.muForward=1e4;c.muBackward=1e4;end
 r=paperbfn.run(c);j=convert_config(c);
 cases{end+1}=struct('name',c.name,'config_json',jsonencode(j),'initial',r.initialReference,'terminal',terminal(c),'final',r.finalState,'errors',r.diagnostics(:,6:11),'endpoints',r.cycleErrors);
end
s=load(fullfile(root,'data','NSE_reference_and_forcing.mat'),'psi_hat_initial','f_hat');n=64;N=size(s.psi_hat_initial,1);take=[1:n/2 N-n/2+1:N];factor=(n/N)^2;
psi_hat_initial=s.psi_hat_initial(take,take)*factor;f_hat=s.f_hat(take,take)*factor;restart=fullfile(folder,'restart.mat');save(restart,'psi_hat_initial','f_hat');
schemes={'standard','diffusive','voigt','diffusive','voigt','standard','diffusive','voigt'};
for i=1:8
 c=paperbfn.config('nse');c.N=n;c.nonlinearCutoff=n/3;c.observe.cutoff=12;c.initial.restartFile=restart;c.window=.004;c.cycles=2;c.output.plot=false;c.output.root=folder;
 c.name=sprintf('nse_%d',i);c.backward.scheme=schemes{i};c.backward.filtered=ismember(i,[4 5]);c.observe.full=i>=6;
 r=paperbfn.run(c);j=convert_config(c);j.nonlinear_cutoff=21;j.allow_restart_resample=true;
 cases{end+1}=struct('name',c.name,'config_json',jsonencode(j),'initial',r.initialReference,'terminal',terminal(c),'final',r.finalState,'errors',r.diagnostics(:,6:11),'endpoints',r.cycleErrors);
end
variants={'BFN','damped BFN','diffusive BFN','Voigt BFN','BFN','diffusive BFN','Voigt BFN'};
models={'damped','damped','damped','damped','viscous','viscous','viscous'};schemes={'standard','damping_reversal','diffusive','voigt','standard','diffusive','voigt'};
for i=1:7
 o=struct('model',models{i},'variant',variants{i},'N',128,'T',.0002,'cycles',2,'sampleEvery',2,'makePlot',false,'saveOutput',false);r=paperbfn.kdv_run(o);
 j=struct('name',sprintf('kdv_%d',i),'model',['kdv_' models{i}],'scheme',schemes{i},'n',128,'nonlinear_cutoff',42,'dt',1e-5,'window',.0002,'cycles',2,'spinup',0,'integrator','ifrk4','feedback_factor',true,'initialization','exact','observation_cutoff',15,'high_mode',20,'high_amplitude',.1,'transport_speed',0,'stride',2,'nu',0,'gamma',0);
 if i<=4,j.gamma=.5;j.alpha=1e-7;j.diffusion_correction=1e-5;else,j.nu=.001;j.alpha=.01;j.diffusion_correction=.001;end
 ep=table2array(r.endpoints(:,1:7));
 cases{end+1}=struct('name',j.name,'config_json',jsonencode(j),'initial',r.initialReferenceHat(:),'terminal',r.terminalReferenceHat(:),'final',r.lastEstimateHat(:),'errors',table2array(r.history(:,4:9)),'endpoints',ep);
end
cases=[cases{:}];save(fullfile(dest,'matlab_short_runs.mat'),'cases','-v7');
fprintf('Wrote %d MATLAB fixture cases\n',numel(cases));
end

function j=convert_config(c)
j=struct('name',c.name,'model',c.model,'profile',c.scenario,'n',c.N,'length',c.L,'dt',c.dt,'window',c.window,'cycles',c.cycles,'nu',c.nu,'gamma',c.gamma,'transport_speed',c.transportSpeed,'dispersion',c.dispersion,'integrator',c.integrator,'feedback_factor',c.feedbackFactor,'nonlinear_cutoff',c.nonlinearCutoff,'mu_forward',c.muForward,'mu_backward',c.muBackward,'observation_cutoff',c.observe.cutoff,'full_observations',c.observe.full,'components',c.observe.components,'fractions',c.observe.fractions,'schedule',c.observe.schedule,'high_mode',c.initial.highMode,'high_amplitude',c.initial.highAmplitude,'spinup',c.initial.spinup,'initialization',c.initial.guessType,'reference_initial',c.initial.reference(:)','guess_initial',c.initial.guess(:)','sigma',c.lorenz.sigma,'rho',c.lorenz.rho,'beta',c.lorenz.beta,'scheme',c.backward.scheme,'filtered',c.backward.filtered,'alpha',c.backward.alpha,'viscosity_factor',c.backward.viscosityFactor,'damping_factor',c.backward.dampingFactor,'viscosity_cutoff',c.backward.viscosityCutoff,'stride',c.output.diagnosticStride);
end

function u=terminal(c)
c=paperbfn.validate(c);m=paperbfn.model(c);t=paperbfn.reference(c,m);u=reshape(t(:,end),m.shape);
end
