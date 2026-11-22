function plot_nse_paper(rs,folder)
%PLOT_NSE_PAPER The three NSE assets used by the manuscript.
% Assemble the three NSE manuscript plot assets from completed numerical runs.
% This routine contains no integration and never estimates errors from images.
% rs is the ordered result list from suite_nse.m. Its first five elements are
% BFN, dBFN, VBFN, filtered dBFN and filtered VBFN.
labels={'BFN','dBFN','VBFN','filtered dBFN','filtered VBFN'};
colors=[.086 .502 .29;0 .447 .698;.835 .369 0;0 .447 .698;.835 .369 0];
% Two readable stacked-axis figures. The first compares the three base models;
% the second adds filtered variants, drawn dashed in their model's color.
% Upper/lower axes show observed/unobserved velocity H1-seminorm errors.
for count=[3 5]
 f=figure('Visible','off','Position',[100 100 1100 1100]);tiledlayout(2,1);
 for part=1:2
  nexttile;hold on;
  for i=1:count
   % The diagnostic columns are 10=observed H1 and 11=unobserved H1.
   % Values <=0 are omitted from the displayed log curve, not the saved data.
   D=rs{i}.diagnostics;y=D(:,9+part);y(y<=0)=NaN;style='-';if i>3,style='--';end
   semilogy(D(:,1),y,style,'Color',colors(i,:),'LineWidth',1.7,'DisplayName',labels{i});
  end
  set(gca,'YScale','log');grid on;xlabel('Completed-cycle coordinate');ylabel('Velocity H^1 seminorm error');
  if part==1,title('Observed modes');else,title('Unobserved modes');end;legend('Location','best');
 end
 name='nse_H1_histories';if count==5,name='nse_filtered_H1_histories';end
 paperbfn.save_figure(f,folder,name);
end
% Kinetic-energy spectrum of the INITIAL reference velocity. Start with
% streamfunction Fourier coefficients and multiply by |dx|^2+|dy|^2 to
% obtain velocity power. Parseval contributes volume/count^2. The factor .5
% is the kinetic-energy convention, not a normalization of the error norms.
r=rs{1};c=r.config;m=paperbfn.model(c,true);z=r.initialReference;
q=(abs(m.dx).^2+abs(m.dy).^2).*abs(z).^2;
% Bin by integer radial shells, summing all 2D Fourier directions in a shell.
% The plotted N/3 marker describes a per-coordinate nonlinear input cutoff;
% it is not a claim that the rectangular input mask has a circular boundary.
shell=floor(m.radius(:));power=accumarray(shell+1,.5*m.volume/m.count^2*q(:));
k=(0:numel(power)-1)';use=k>0 & power>0;
f=figure('Visible','off','Position',[100 100 1100 600]);loglog(k(use),power(use),'k','LineWidth',1.7);
xline(c.N/3,':','Per-coordinate cutoff');xlabel('Integer radial wave-number shell');ylabel('Kinetic energy per shell');grid on;
paperbfn.save_figure(f,folder,'nse_initial_spectrum');
end
