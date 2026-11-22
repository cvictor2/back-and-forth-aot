function compare(results,labels,folder,name,level,part)
%COMPARE Overlay fresh numerical histories without replacing archival data.
% Overlay one diagnostic from several fresh result structures. Does not solve
% an equation or change any recorded numerical error. The suites use this for
% parameter families and observed/unobserved NSE comparisons.
%
% level=1 selects L2, level=2 selects H1 seminorm. part=1,2,3 selects total,
% observed,unobserved respectively. The defaults are total L2.
if nargin<5,level=1;end
if nargin<6,part=1;end
f=figure('Visible','off','Position',[100 100 1100 650]);hold on;
for i=1:numel(results)
 % Skip the first five time/pass bookkeeping columns to select the requested
 % norm. Hide nonpositive values only in the local log-plot copy.
 r=results{i};D=r.diagnostics;y=D(:,5+(level-1)*3+part);
 y(y<=0 | ~isfinite(y))=NaN;
 semilogy(D(:,1),y,'LineWidth',1.7,'DisplayName',labels{i});
end
set(gca,'YScale','log');grid on;xlabel('Completed-cycle coordinate');
ylabel('L^2 error');if strcmp(results{1}.config.model,'lorenz'),ylabel('Euclidean error');end
if level==2,ylabel('Velocity H^1 seminorm error');end
legend('Location','best','Interpreter','none');
paperbfn.save_figure(f,folder,name);
end
