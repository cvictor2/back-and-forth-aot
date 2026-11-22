function summary=run_lorenz(mode)
%RUN_LORENZ Lorenz experiments. Modes 'paper' (default), 'smoke', 'archive'.
% Public selector for invariant, parameter-sweep and observation-window runs.
% Open +paperbfn/suite_lorenz.m to see each experiment's actual parameters.
% reference.m generates the unnudged target. run.m manages BFN cycles, while
% lorenz_leg.m contains the physical RK4 and separate explicit feedback loop.
% The variables use a shifted z coordinate, explained beside its rhs helper.
%
% Examples from this repository folder
%   run_lorenz('smoke')       short runs covering every experiment family
%   run_lorenz('paper')       full sweeps with half-open observation windows
%   run_lorenz('historical')  alternate closed-endpoint observation schedule
% summary lists output folders and statuses. Archive mode only redraws
% historical display coordinates, without evolving the Lorenz equations.
% Baseline sigma=10, rho=28, b=8/3, confirmed by the author.
if nargin<1,mode='paper';end
summary=paperbfn.entry('lorenz',mode);
end
