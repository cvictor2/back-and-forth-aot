function summary=suite_lorenz(mode,folder)
%SUITE_LORENZ Invariant, parameter and observation-window experiments.
% Run the invariant example, parameter sweeps and observation-window tests.
% This file selects experiments. The physical RK4 and explicit feedback are
% in lorenz_leg.m; the unnudged reference is generated through reference.m.
% Coordinates are shifted Lorenz variables, as explained beside rhs there.
historical=strcmp(mode,'historical');
smoke=strcmp(mode,'smoke');base=paperbfn.config('lorenz');
% Individual runs share a suite output root but have separate result folders.
% Plotting is deferred so related trajectories can be overlaid by compare.m.
% Resume checks the configuration, source hashes and any input-file hashes.
% A changed solver or changed input therefore requires a fresh calculation.
base.output.root=folder;base.output.plot=false;base.output.resume=true;
if smoke,base.window=.006;base.dt=1e-5;base.cycles=2;base.output.diagnosticStride=10;end
rows={};
% Invariant family, five cycles and the three amplitudes in Figure 1.
rs={};labels={};
% Reference and estimate lie in the x=y=0 invariant plane, differing initially
% only in the unobserved shifted z coordinate. Changing phi changes this initial
% error amplitude while leaving the model and observed component pair fixed.
for phi=[1 1e-5 1e-10]
 c=base;c.name=sprintf('invariant_%g',phi);c.cycles=5;
 if smoke,c.cycles=2;end
 c.initial.reference=[0;0;0];c.initial.guess=[0;0;phi];c.muForward=1000;c.muBackward=1000;
 [r,path]=paperbfn.run(c);rs{end+1}=r;labels{end+1}=sprintf('phi = %g',phi);rows(end+1,:)={c.name,path,r.status};
end
paperbfn.compare(rs,labels,folder,'lorenz_error');
% Parameter families. The rho and beta sweeps use sigma=20.
% Beta values are specified as numerators over 3, and their labels display
% that division explicitly. Other parameters stay at the base values.
for pair={{'sigma',[10 20 30 40 50 100]},{'rho',[10 20 28 30 40 50 100]},{'beta',[3 6 8 9 12 15 30]}}
 group=pair{1};field=group{1};rs={};labels={};
 for value=group{2}
  c=base;c.lorenz.(field)=value;
  if strcmp(field,'rho'),c.lorenz.sigma=20;end
  if strcmp(field,'beta'),c.lorenz.sigma=20;c.lorenz.beta=value/3;end
  c.name=sprintf('%s_%g',field,value);
  [r,path]=paperbfn.run(c);rs{end+1}=r;label=sprintf('%s = %g',field,value);
  if strcmp(field,'beta'),label=sprintf('b = %g/3',value);end
  labels{end+1}=label;rows(end+1,:)={c.name,path,r.status};
 end
 paperbfn.compare(rs,labels,folder,['Lorenz_' field]);
end
% Window experiments observe all three components, each in its own third
% of physical time. One family varies gain with full third-windows. The other
% varies observation fraction at fixed gain 1e4, including zero observations.
for family={'gain','fraction'}
 rs={};labels={};
 if strcmp(family{1},'gain'),values=[10 100 1000 10000 100000];else,values=[0 1/8 1/4 1/2 3/4 7/8 15/16 31/32 63/64 1];end
 for value=values
  c=base;c.cycles=10;if smoke,c.cycles=2;end
  c.observe.components=1:3;c.observe.schedule='thirds';
  % Publication mode uses current-time masks with half-open intervals.
  % The historical option selects closed-right observation windows and a
  % backward mask evaluated one time step before the current physical time.
  if historical,c.lorenz.windowRule='legacyClosed';c.lorenz.backwardObservationShift=-1;end
  c.muForward=1e4;c.muBackward=1e4;
  if strcmp(family{1},'gain'),c.muForward=value;c.muBackward=value;
  else,c.observe.fractions=value*[1 1 1];end
  c.name=sprintf('window_%s_%g',family{1},value);
  [r,path]=paperbfn.run(c);rs{end+1}=r;labels{end+1}=sprintf('%s = %g',family{1},value);rows(end+1,:)={c.name,path,r.status};
 end
 name='discrete_data';if strcmp(family{1},'fraction'),name='discrete_gaps';end
 paperbfn.compare(rs,labels,folder,name);
end
summary=cell2table(rows,'VariableNames',{'experiment','folder','status'});
writetable(summary,fullfile(folder,'runs.csv'));
end
