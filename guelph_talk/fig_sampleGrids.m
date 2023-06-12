clear all
clf
close all
addpath('utils_core')
addpath('utils_cost_fun')

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
fig_width = 4; % width in inches
fig_height = 3; % height in inches
fig_aspect_ratio = fig_width / fig_height;
fig = figure;
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [fig_width fig_height], ...
         'PaperPosition', [0 0 fig_width fig_height], ...
         'PaperPositionMode', 'manual', ...
         'Units', 'inches', ...
         'Position', [0 0 fig_width fig_height]);

font_size = 10;
set(0, 'DefaultAxesFontSize', font_size);
set(0, 'DefaultTextFontSize', font_size);

def_colours

rng(2341)
p.Nmeas=8;
N=p.Nmeas;
tauL=0;
tauR=0.25;
[t_unif,t_nu]=getSamplingSchedules(4,4,tauL,tauR);
t_rand=sort(rand(1,N));

t_jit = t_unif + 2.5e-2*rand(1,N);

plot(t_unif,4,'ok','MarkerEdgeColor','none','MarkerFaceColor',c1)

hold on
plot(t_nu,3,'ok','MarkerEdgeColor','none','MarkerFaceColor',c2)
plot(t_rand,2,'ok','MarkerEdgeColor','none','MarkerFaceColor',c3)
plot(t_jit,1,'ok','MarkerEdgeColor','none','MarkerFaceColor',c4)

ylim([0.5 4.5])
xlabel('$t$ (hr)','interpreter','latex')
yticks([1 2 3 4])
yticklabels({'jittered','random','2-uniform','uniform'})
xticks([0 0.5 1])
xticklabels({'0','12','24'})

print(gcf,strcat('/home/turner/research/overleaf/guelph_talk/figs_csc/','sampleGrid.png'),'-dpng','-r600') % -r sets the resolution

