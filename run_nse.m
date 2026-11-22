function summary=run_nse(mode)
%RUN_NSE Partial, filtered and full observations. 'paper', 'smoke', 'archive'.
% Public NSE entry point. This file selects a mode and delegates the experiment
% list to +paperbfn/suite_nse.m. It contains no hidden integration formula.
%
% To inspect the mathematics, open these helpers in +paperbfn
%   model.m       Fourier multipliers, initial streamfunction and forcing
%   reference.m   physical forward trajectory without nudging
%   run.m         forward/backward time loops for one experiment
%   step.m        integrating-factor Euler and explicit feedback
%   nonlinear.m   velocities and the Basdevant NSE nonlinearity
%   errors.m      velocity L2 and H1-seminorm error calculations
%
% Examples from this repository folder
%   run_nse('smoke')   short/small runs for all eight variants
%   run_nse('paper')   full paper settings, the default
%   run_nse('archive') historical plot-data redraw, not a simulation
% summary lists output folders and statuses. To inspect a step, set a breakpoint
% in step.m before calling smoke mode. Initial calls generate the reference
% (truthOnly=true), followed by calls with assimilation feedback enabled.
if nargin<1,mode='paper';end
summary=paperbfn.entry('nse',mode);
end
