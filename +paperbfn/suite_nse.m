function summary=suite_nse(mode,folder)
%SUITE_NSE Three backward models, filtered variants and full observations.
% This file is the readable list of NSE experiments. It sets parameters and
% calls the solver; it does not contain the time integrator itself.
%
% To follow one run, open these files in order
%   config.m       default NSE grid, dt, viscosity and reference path
%   model.m        Fourier geometry, psi initial data and forcing
%   reference.m    unnudged forward target at every discrete physical time
%   run.m          forward/backward cycle loops
%   step.m         linear factor, explicit physical step and old-time feedback
%   nonlinear.m    velocities and Basdevant nonlinear streamfunction RHS
%   errors.m       velocity L2 and H1-seminorm diagnostics
%
% Paper mode uses the unchanged full settings. Smoke mode changes the grid,
% window and cycles only as shown below for a quick installation/review run.
smoke=strcmp(mode,'smoke');base=paperbfn.config('nse');base.output.root=folder;
if smoke
 % Spectral restriction of the supplied restart, not a synthetic target.
 a=load(base.initial.restartFile,'psi_hat_initial','f_hat');n=64;N=size(a.psi_hat_initial,1);
 % Retain the low positive/negative Fourier indices in each coordinate.
 % MATLAB FFT coefficients scale like grid point count, so restriction from
 % N-by-N to n-by-n needs (n/N)^2 for the same physical Fourier amplitudes.
 % This is an explicitly smaller test problem, not a fresh random target.
 take=[1:n/2 N-n/2+1:N];factor=(n/N)^2;
 psi_hat_initial=a.psi_hat_initial(take,take)*factor;f_hat=a.f_hat(take,take)*factor;
 smokeFile=fullfile(folder,'smoke_restart.mat');save(smokeFile,'psi_hat_initial','f_hat');
 base.N=n;base.initial.restartFile=smokeFile;base.observe.cutoff=12;
 base.window=.004;base.dt=.001;base.cycles=2;base.nonlinearCutoff=n/3;
end
% Eight configurations in a fixed order. First compare standard, diffusive
% and Voigt BFN on the observed band radius<50. Then add their filtered
% regularizations. Finally observe every nonzero mode for three variants.
schemes={'standard','diffusive','voigt','diffusive','voigt','standard','diffusive','voigt'};
labels={'BFN','dBFN','VBFN','filtered dBFN','filtered VBFN','full BFN','full dBFN','full VBFN'};
results={};rows={};endpoints=nan(8,6);
% Copy base settings for each experiment so a previous variant cannot leak
% into the next. alpha and viscosity come from config('nse'). Default dBFN
% reverses the signed-time viscous multiplier, giving backward diffusion.
for i=1:8
 c=base;c.name=strrep(labels{i},' ','_');c.backward.scheme=schemes{i};
 % Runs 4/5 regularize only unobserved modes. Runs 6/7/8 use full observations.
 % Initial projection follows the selected observation mask for EACH run.
 c.backward.filtered=ismember(i,[4 5]);c.observe.full=i>=6;
 % This is where integration starts. Set a breakpoint in step.m to inspect
 % an individual NSE update. The first calls build the unnudged reference;
 % subsequent calls advance the assimilated streamfunction.
 [r,path]=paperbfn.run(c);results{end+1}=r;rows(end+1,:)={c.name,path,r.status};
 if r.completedCycles>0,endpoints(i,:)=r.cycleErrors(end,2:7);end
end
% Plot error histories from raw numerical diagnostics. part=2 selects observed
% and part=3 unobserved errors. level=2 selects the velocity H1 seminorm.
% These comparisons do not rerun or modify the simulation.
for part=2:3
 suffix='observed';if part==3,suffix='unobserved';end
 paperbfn.compare(results(1:3),labels(1:3),folder,['nse_H1_' suffix],2,part);
 paperbfn.compare(results(1:5),labels(1:5),folder,['nse_filtered_H1_' suffix],2,part);
end
% Produce the three manuscript assets, including readable stacked H1 panels
% and the kinetic-energy spectrum. The endpoint table below includes all six
% errors for each method at the last completed backward endpoint.
paperbfn.plot_nse_paper(results,folder);
T=array2table(endpoints,'VariableNames',{'L2','observedL2','unobservedL2','H1','observedH1','unobservedH1'});
T=addvars(T,labels','Before',1,'NewVariableNames','method');writetable(T,fullfile(folder,'nse_endpoints.csv'));
summary=cell2table(rows,'VariableNames',{'experiment','folder','status'});
writetable(summary,fullfile(folder,'runs.csv'));
end
