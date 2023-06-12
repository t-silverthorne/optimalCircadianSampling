close all
fname='unif_non_unif';
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
fig_width = 2.5; % width in inches
fig_height = 2.1; % height in inches
fig_aspect_ratio = fig_width / fig_height;
fig = figure;
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [fig_width fig_height], ...
         'PaperPosition', [0 0 fig_width fig_height], ...
         'PaperPositionMode', 'manual', ...
         'Units', 'inches', ...
         'Position', [0 0 fig_width fig_height]);

def_colours
addpath('utils_core/')

%%
xv=0:.1:24;
tiledlayout(2,1,'TileSpacing','tight')
nexttile(1)
f=@(xv) sin(2*pi*xv/24-pi/2)+0.5*exp(-(xv-20).^2/5).*cos(2*pi*xv/2);
plot(xv,f(xv),'-k')

hold on
c1="#648FFF"; % uniform colour
[xunif,x2u]=getSamplingSchedules(6,6,0.75,.9);
xunif=24*xunif;
x2u=24*x2u;
m1=plot(xunif,f(xunif),'.k','color',c1,'DisplayName','Uniform','MarkerSize',10);
xticks(0:12:24)
set(gca,'XTickLabel',[]);
xlim([0,24])
nexttile(2)


plot(xv,f(xv),'-k')
hold on
m2=plot(x2u,f(x2u),'.k','color',c2,'DisplayName','Non-uniform','MarkerSize',10);
xticks(0:12:24)
leg=legend([m1,m2],'NumColumns',2);
leg.Layout.Tile='south';
xlim([0,24])


print(gcf,strcat('/home/turner/research/overleaf/guelph_talk/figs_csc/',fname) ...
    ,'-dpng','-r600') % -r sets the resolution