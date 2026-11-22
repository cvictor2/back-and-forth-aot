function plot_result(r,folder)
%PLOT_RESULT Separate full-width error, state and spectrum plots.
% Read saved arrays and create error/state/spectrum figures. No PDE is solved.
% This is the generic plotting path for Lorenz, transport, Burgers and NSE.
% KdV has its own plotting block in kdv_run.m.
%
% r.reconstruction is the last completed backward endpoint at physical t=0.
% It can differ from r.finalState if integration failed partway through a pass.
% The reference used in state/spectrum plots is likewise the t=0 reference.
c=r.config;m=paperbfn.model(c,true);D=r.diagnostics;
% Columns 6:8 of diagnostics are L2 errors, and 9:11 are H1-seminorm errors.
% Only NSE makes an additional H1 figure. Within a triplet the order is total,
% observed, unobserved. Each figure uses completed-cycle coordinate on x.
levels=1;if strcmp(c.model,'nse'),levels=2;end
for level=1:levels
 f=figure('Visible','off','Position',[100 100 1100 600]);
 % Hide nonpositive/nonfinite values in the plotting COPY because log axes
 % cannot show them. Do not impose a floor on saved errors or endpoint tables.
 y=D(:,6+(level-1)*3:8+(level-1)*3);y(y<=0 | ~isfinite(y))=NaN;
 semilogy(D(:,1),y,'LineWidth',1.7);grid on;
 xlabel('Completed-cycle coordinate');ylabel('L^2 error');
 if strcmp(c.model,'lorenz'),ylabel('Euclidean error');end
 if level==2,ylabel('Velocity H^1 seminorm error');end
 legend('Total','Observed','Unobserved','Location','best');
 title(sprintf('%s, %s',c.name,r.status),'Interpreter','none');
 paperbfn.save_figure(f,folder,sprintf('errors_%d',level));
end
if strcmp(c.model,'lorenz'),return;end
% A 1D state is transformed back to physical space for direct comparison.
% The NSE state is streamfunction and has no analogous 1D physical-state panel.
u=r.initialReference;v=r.reconstruction;
if ~strcmp(c.model,'nse')
 f=figure('Visible','off','Position',[100 100 1100 600]);
 x=(0:c.N-1)'*c.L/c.N;
 plot(x,real(ifft(u)),'-',x,real(ifft(v)),'--','LineWidth',1.7);
 xlabel('x');ylabel('State at physical time zero');grid on;
 legend('Reference','Last completed backward endpoint','Location','best');
 paperbfn.save_figure(f,folder,'initial_state');
end
f=figure('Visible','off','Position',[100 100 1100 600]);
[k,pu]=spectral_power(u,m);[~,pv]=spectral_power(v,m);
loglog(k,positive(pu),'-',k,positive(pv),'--','LineWidth',1.7);grid on;
xlabel('Integer Fourier mode or radial shell');ylabel('Squared L^2 contribution');
legend('Reference','Last completed backward endpoint','Location','best');
xline(c.observe.cutoff-1,':','HandleVisibility','off');
if ~strcmp(c.model,'nse')
 xline(c.nonlinearCutoff,':','HandleVisibility','off');
 peak=max([pu(:);pv(:)]);if peak>0,ylim([max(peak*1e-18,1e-300),peak*4]);end
end
paperbfn.save_figure(f,folder,'initial_spectrum');
end
% Accumulate squared integral-L2 contributions by integer mode/radial shell.
% Both positive and negative 1D modes contribute to the same absolute shell.
% For NSE, convert streamfunction to velocity power before grouping by floor
% of radial mode magnitude. This plot is L2 power, not the .5 energy spectrum
% made specifically by plot_nse_paper.m.
function [k,p]=spectral_power(z,m)
q=abs(z).^2;
if strcmp(m.name,'nse'),q=(abs(m.dx).^2+abs(m.dy).^2).*q;end
shell=floor(m.radius(:));p=accumarray(shell+1,m.volume/m.count^2*q(:));
k=(1:numel(p)-1)';p=p(2:end);
end
% Display helper only. Missing log points become NaN instead of fake small
% positive data. The numerical arrays in result.mat remain untouched.
function y=positive(y)
y(y<=0 | ~isfinite(y))=NaN;
end
