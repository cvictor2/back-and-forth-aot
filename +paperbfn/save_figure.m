function save_figure(f,folder,name)
%SAVE_FIGURE Wide, readable plots in editable and publication formats.
% Shared export helper. PNG is a raster preview, PDF is vector publication
% output, and FIG retains MATLAB's editable curves and numerical plot data.
% All three exports describe the same figure. Close it afterward so full
% suites do not accumulate many open graphics windows.
if ~isfolder(folder),mkdir(folder);end
set(f,'Color','w');
% Apply readable axes formatting consistently. This changes presentation only,
% not stored states, error values or the numerical data behind plotted lines.
axs=findall(f,'Type','axes');set(axs,'FontSize',13,'LineWidth',0.8,'XMinorGrid','off','YMinorGrid','off');
exportgraphics(f,fullfile(folder,[name '.png']),'Resolution',250);
exportgraphics(f,fullfile(folder,[name '.pdf']),'ContentType','vector');
savefig(f,fullfile(folder,[name '.fig']));close(f);
end
