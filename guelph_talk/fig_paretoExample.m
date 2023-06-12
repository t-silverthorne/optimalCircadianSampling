def_colours
close all
fname='paretoExample'
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
fig_width = 5; % width in inches
fig_height = 3; % height in inches
fig_aspect_ratio = fig_width / fig_height;
fig = figure;
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [fig_width fig_height], ...
         'PaperPosition', [0 0 fig_width fig_height], ...
         'PaperPositionMode', 'manual', ...
         'Units', 'inches', ...
         'Position', [0 0 fig_width fig_height]);

N=1e3;
X=1+randn(N,2);
hvis='on';
for ii=1:N
    if X(ii,1)*X(ii,2)>1 && X(ii,1)>0 && X(ii,2)>0
        plot(X(ii,1),X(ii,2),'.k','color',c3,HandleVisibility=hvis)
        if strcmp(hvis,'on')
            hvis='off';
        end
        hold on
    end
end
tvals=1e-3:1e-3:5;
plot(tvals,1./tvals,'--k','LineWidth',1)
xlim([0,5])
ylim([0,5])
xlabel('Criterion 1')
ylabel('Criterion 2')

legend({'compromise','Pareto frontier'}, ...
    'Location','southoutside','NumColumns',2)

print(gcf,strcat('/home/turner/research/overleaf/guelph_talk/figs_csc/',fname) ...
    ,'-dpng','-r600') % -r sets the resolution