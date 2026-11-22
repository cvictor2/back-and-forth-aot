function summary=entry(group,mode,selection)
%ENTRY Each invocation has an isolated output folder and complete settings.
% Shared dispatcher used by the three public run_*.m entry points.
% It chooses an experiment suite and creates an isolated output directory.
% The numerical loop lives in run.m / lorenz_leg.m / kdv_run.m.
%
% paper       fresh full-resolution manuscript experiments
% smoke       fresh short/small experiments with every requested variant
% archive     redraw the supplemental display-coordinate data, no integration
% historical  Lorenz with the alternate closed-endpoint observation schedule
if nargin<3,selection='all';end
mode=validatestring(mode,{'paper','smoke','archive','historical'});
assert(~strcmp(mode,'historical') || strcmp(group,'lorenz'),'Historical mode is currently defined only for Lorenz.');
root=fileparts(fileparts(mfilename('fullpath')));
base=fullfile(root,'results');if ~isfolder(base),mkdir(base);end
folder=tempname(base);mkdir(folder);
fprintf('%s %s outputs: %s\n',group,mode,folder);
% Archive and fresh simulation paths are deliberately distinct. Archive plots
% cannot validate the solver because they do not evaluate any model equations.
if strcmp(mode,'archive')
 summary=paperbfn.replot_archive(group,folder);
else
 switch group
  case 'lorenz',summary=paperbfn.suite_lorenz(mode,folder);
  case 'one_d',summary=paperbfn.suite_1d(mode,folder,selection);
  case 'nse',summary=paperbfn.suite_nse(mode,folder);
 end
end
% Check each expected manuscript asset was exported in FIG, PDF and PNG.
% This verifies output coverage, not accuracy or stability of the model solver.
if ~strcmp(mode,'archive'),paperbfn.check_figure_coverage(group,folder,selection);end
fid=fopen(fullfile(folder,'run_mode.txt'),'w');fprintf(fid,'%s\n%s\n',mode,group);fclose(fid);
end
