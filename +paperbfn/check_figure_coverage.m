function coverage=check_figure_coverage(group,folder,selection)
%CHECK_FIGURE_COVERAGE Account for every manuscript asset expected from a suite.
% Compare the suite's root plot exports with FIGURE_MANIFEST.csv.
% This ensures every expected panel has editable FIG, PDF and PNG files.
% Numerical agreement and the validity of a reconstruction are assessed
% separately from saved states and error histories.
if nargin<3,selection='all';end
root=fileparts(fileparts(mfilename('fullpath')));
M=readtable(fullfile(root,'FIGURE_MANIFEST.csv'),'Delimiter',',', ...
 'ReadVariableNames',true,'TextType','string','VariableNamingRule','preserve');
M=M(strcmp(M.group,group),:);
% A selected 1D family is responsible only for its own panels. Optional heat
% examples have no manuscript assets, so their manifest subset is empty.
if strcmp(group,'one_d') && ~strcmp(selection,'all')
 switch selection
  case 'transport',keep=startsWith(M.asset,'IT_') | startsWith(M.asset,'VT_');
  case 'burgers',keep=startsWith(M.asset,'IB_') | startsWith(M.asset,'VB_');
  case 'kdv',keep=startsWith(M.asset,'bfn_');
  otherwise,keep=false(height(M),1); % optional heat has no manuscript panel
 end
 M=M(keep,:);
end
% Check formats independently so a missing editable FIG is not concealed by
% the existence of a PNG. Keep a CSV record of coverage for this invocation.
png=false(height(M),1);pdf=png;fig=png;
for i=1:height(M)
 png(i)=isfile(fullfile(folder,M.asset(i)+'.png'));
 pdf(i)=isfile(fullfile(folder,M.asset(i)+'.pdf'));
 fig(i)=isfile(fullfile(folder,M.asset(i)+'.fig'));
end
coverage=table(M.figure,M.asset,png,pdf,fig, ...
 'VariableNames',{'manuscriptFigure','asset','hasPNG','hasPDF','hasFIG'});
writetable(coverage,fullfile(folder,'figure_coverage.csv'));
assert(all(png & pdf & fig),'paperbfn:MissingFigure', ...
 'An expected manuscript figure export is missing. See figure_coverage.csv.');
fprintf('Verified %d manuscript plot assets in FIG, PDF and PNG formats.\n',height(M));
end
