function summary=run_1d(mode,selection)
%RUN_1D Transport, Burgers and KdV. Modes 'paper', 'smoke', 'archive'.
% Public selector for transport, Burgers, optional heat and all seven KdV panels.
% The experiment list and parameter overrides are in +paperbfn/suite_1d.m.
% Transport/Burgers/heat use run.m -> step.m -> nonlinear.m. KdV instead uses
% kdv_run.m, whose physicalStep helper contains the integrating-factor RK4.
%
% Examples from this repository folder
%   run_1d('smoke','burgers') quick inviscid/viscous three-profile checks
%   run_1d('paper','kdv')     the seven full-resolution KdV experiments
%   run_1d('paper')           all paper 1D families
%   run_1d('paper','heat')    optional numerical heat examples
% summary reports each experiment's output folder and status. "archive" redraws
% historical saved plot data and must not be used to validate the integrator.
if nargin<1,mode='paper';end
if nargin<2,selection='all';end
summary=paperbfn.entry('one_d',mode,selection);
end
