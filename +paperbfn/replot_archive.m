function summary=replot_archive(group,folder)
%REPLOT_ARCHIVE Redraw saved display coordinates without solving a PDE.
% Display-data viewer. This branch does not call the physical solver.
% It redraws the supplemental archive rather than the current paper results.
% Use paper mode for numerical integration and reference_figures/ for the
% manuscript plot set. Archive datasets can use different observation schedules.
% KdV entries in this archive are raster panels without trajectory arrays.
% A saved-image redraw does not provide a numerical test of the solver.
root=fileparts(fileparts(mfilename('fullpath')));
r=jsondecode(fileread(fullfile(root,'data','manuscript_plot_data.json')));
rows={};
for i=1:numel(r)
 name=r(i).name;
 switch group
  case 'lorenz',use=contains(lower(name),'lorenz') || startsWith(name,'discrete');
  case 'one_d',use=any(startsWith(name,{'IT_','VT_','IB_','VB_'}));
  case 'nse',use=startsWith(name,'nse_');
 end
 if ~use,continue;end
 f=figure('Visible','off','Position',[100 100 1100 650*numel(r(i).axes)]);
 tiledlayout(numel(r(i).axes),1,'TileSpacing','compact');
 % Restore saved axes, curve coordinates, scales and labels from the archive.
 % Vertical lines denote plot annotations rather than new simulation data.
 for j=1:numel(r(i).axes)
  a=r(i).axes(j);ax=nexttile;hold on;
  for q=reshape(a.lines,1,[])
   if q.vertical
    xline(q.x(1),'Color',q.color,'LineStyle',q.ls,'HandleVisibility','off');
   else
    plot(q.x,q.y,'Color',q.color,'LineStyle',q.ls,'LineWidth',max(q.lw,1.7), ...
      'DisplayName',q.label);
   end
  end
  set(ax,'XScale',a.xscale,'YScale',a.yscale);xlim(a.xlim);ylim(a.ylim);
  xticks(a.xticks);xticklabels(a.xticklabels);
  xlabel(a.xlabel,'Interpreter','latex');ylabel(a.ylabel,'Interpreter','latex');
  title(a.title,'Interpreter','latex');grid on;
  if a.legend,legend('Location','best','Interpreter','latex');end
 end
 [~,stem]=fileparts(name);paperbfn.save_figure(f,folder,stem);rows(end+1,:)={name,'saved display data'};
end
% The KdV archive contains PNG panels only. Copying a raster cannot
% supply hidden numerical values or an editable numerical FIG. New KdV runs
% produce their own full states, error histories and FIG exports separately.
if strcmp(group,'one_d')
 files=dir(fullfile(root,'data','kdv_original_panels','*.png'));
 for i=1:numel(files)
  copyfile(fullfile(files(i).folder,files(i).name),fullfile(folder,files(i).name));
  rows(end+1,:)={files(i).name,'original raster only'};
 end
end
summary=cell2table(rows,'VariableNames',{'figure','source'});
writetable(summary,fullfile(folder,'archive_manifest.csv'));
end
