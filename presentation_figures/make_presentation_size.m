close all
fig_list={'atwood_design.fig','hist_sim_result.fig','linear_regression_histogram.fig', ...
          'modified_wynn.fig','nonlinear_regression_histogram_top.fig','samplegrid_fig.fig', ...
          'simulated_annealing_best.fig','wynn_convergence.fig'}
fig_list{6}
%%
for ii=[6]
    openfig(fig_list{ii})
    f=gcf;
    set(findall(gcf,'-property','FontSize'),'FontSize',45)
%     chil = get(gcf,'children');
%     h=chil(2).Children;
%     for i = 1:length(h)
%         if strcmp(get(h(i),'type'),'line')        % check if it's a line
%             set(h(i),'markersize',25);
%         end
%     end
    lines = findobj(gcf,'Type','ConstantLine');
    for i = 1:numel(lines)
        lines(i).LineWidth = 5;
    end
    lines = findobj(gcf,'Type','line');
    for i = 1:numel(lines)
        lines(i).LineWidth = 5;
        lines(i).MarkerSize = 25;
    end
    lines = findobj(gcf,'Type','scatter');
    for i = 1:numel(lines)
        lines(i).SizeData = 800;
    end
    plot_filename=strcat('pres_version_',fig_list{ii}(1:end-4));
   
    wd=6*5; % width
    ht=3*5; % height
    set(gcf,'PaperUnits','inches')
    set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
    print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
    clf
end
%%
for i = 1:numel(lines)
    lines(i).SizeData = 800;
end

%%
f=figure(2)