close all
fig_list={'atwood_design.fig','hist_sim_result.fig','linear_regression_histogram.fig', ...
          'modified_wynn.fig','nonlinear_regression_histogram_top.fig','samplegrid_fig.fig', ...
          'simulated_annealing_best.fig','wynn_convergence.fig'}

%%
for ii=1:numel(fig_list)
    openfig(fig_list{ii})
    f=gcf;
    set(findall(gcf,'-property','FontSize'),'FontSize',45/5)
%     chil = get(gcf,'children');
%     h=chil(2).Children;
%     for i = 1:length(h)
%         if strcmp(get(h(i),'type'),'line')        % check if it's a line
%             set(h(i),'markersize',25);
%         end
%     end
    lines = findobj(gcf,'Type','ConstantLine');
    for i = 1:numel(lines)
        lines(i).LineWidth = 5/5;
    end
    lines = findobj(gcf,'Type','line');
    for i = 1:numel(lines)
        lines(i).LineWidth = 5/5;
        lines(i).MarkerSize = 25/5;
    end
    lines = findobj(gcf,'Type','scatter');
    for i = 1:numel(lines)
        lines(i).SizeData = 800/5;
    end
    plot_filename=strcat('pres_version_',fig_list{ii}(1:end-4));
   
    wd=6*1; % width
    ht=3*1; % height
    if ii==5 || ii==6
        ht=ht/2;
    end
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
