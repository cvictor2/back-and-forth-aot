function summary=suite_1d(mode,folder,selection)
%SUITE_1D Transport/Burgers examples, optional heat, and seven KdV panels.
% Experiment selector, not the integrator. For transport/Burgers follow
% config -> model -> reference -> run -> step -> nonlinear -> errors.
% KdV instead calls kdv_run.m, which contains its full IF-RK4 implementation.
% "paper" uses the manuscript parameters. "smoke" takes the shorter settings
% shown here and is not a reproduction of the full paper trajectories.
if nargin<3,selection='all';end
selection=validatestring(selection,{'all','transport','burgers','kdv','heat'});
rows={};smoke=strcmp(mode,'smoke');
% Tags I/V mean inviscid/viscous, T/B mean transport/Burgers. The optional
% heat example uses HT. The three profiles map to high-only, low-only and
% mixed initial states, corresponding to none/full/partial observed signal.
models={'transport','transport','burgers','burgers'};tags={'IT','VT','IB','VB'};
profiles={'high','low','mixed'};suffix={'none','full','partial'};
if strcmp(selection,'heat'),models={'heat'};tags={'HT'};end
for i=1:numel(models)
 if ~ismember(selection,{'all','heat'}) && ~strcmp(selection,models{i}),continue;end
 for j=1:3
  c=paperbfn.config(models{i},profiles{j});c.name=[tags{i} '_' suffix{j}];
  % Viscous examples use nu=.05 and backward viscosity-only truncation.
  % The source VT high-only panel has amplitude 1 instead of .05. Preserve
  % this profile-specific amplitude when comparing the model families.
  if mod(i,2)==0,c.nu=.05;c.backward.scheme='cutoff';end
  if i==2 && j==1,c.initial.highAmplitude=1;end
  c.output.root=folder;
  % Reduce resolution/window/cycles for quick checks. These overrides remain
  % explicit so a smoke result cannot be mistaken for paper-resolution data.
  if smoke
   c.N=256;c.nonlinearCutoff=42;c.dt=1e-5;c.window=.0002;c.initial.spinup=.0001;c.cycles=2;
   c.output.diagnosticStride=2;
  end
  [r,path]=paperbfn.run(c);rows(end+1,:)={c.name,path,r.status};
  % Stable paper-asset names, separate from each run's raw result.
  for item={{'errors_1','error'},{'initial_state','soln'},{'initial_spectrum','spec'}}
   p=item{1};
   for ext={'.png','.pdf','.fig'}
    copyfile(fullfile(path,[p{1} ext{1}]),fullfile(folder,[tags{i} '_' p{2} '_' suffix{j} ext{1}]));
   end
  end
 end
end
% Seven KdV panels, with their forward physical model and backward variant.
% The diffusion, damping and alpha defaults are documented inside kdv_run.m.
% The viscous standard-BFN run is allowed to fail at paper resolution.
variants={'BFN','damped BFN','diffusive BFN','Voigt BFN','BFN','diffusive BFN','Voigt BFN'};
models={'damped','damped','damped','damped','viscous','viscous','viscous'};
names={'bfn_normal','bfn_damped','bfn_diffusive','bfn_VOIGT', ...
 'bfn_model_diffusive','bfn_viscous_model_diffusive','bfn_VOIGT_model_diffusive'};
for i=1:7
 if ~ismember(selection,{'all','kdv'}),continue;end
 o=struct('model',models{i},'variant',variants{i},'outputRoot',folder);
 if smoke,o.N=128;o.dt=1e-5;o.T=.0002;o.cycles=2;o.sampleEvery=2;end
 % Run the separate KdV driver, retain its failure status, and copy the error
 % exports under the stable asset name expected by the manuscript manifest.
 r=paperbfn.kdv_run(o);status='completed';if r.failed,status='failed';end
 rows(end+1,:)={names{i},r.outputFolder,status};
 for ext={'.png','.pdf','.fig'}
  copyfile(fullfile(r.outputFolder,['error' ext{1}]),fullfile(folder,[names{i} ext{1}]));
 end
end
% Return a catalog of individual experiment folders and outcomes. Raw states,
% settings and error tables live in those folders rather than only in images.
summary=cell2table(rows,'VariableNames',{'experiment','folder','status'});
writetable(summary,fullfile(folder,'runs.csv'));
end
