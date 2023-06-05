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


xv=0:.1:36;
tiledlayout(2,1,'TileSpacing','tight')
nexttile(1)
f=@(xv) sin(2*pi*xv/24-pi/2)+0.5*exp(-(xv-20).^2/5).*cos(2*pi*xv/2);
plot(xv,f(xv),'-k')

hold on
c1="#648FFF"; % uniform colour
xunif=0:3:36;
m1=plot(xunif,f(xunif),'.k','color',c1,'DisplayName','Uniform','MarkerSize',10)
xticks(0:12:36)
set(gca,'XTickLabel',[]);
nexttile(2)

addpath('utils_core/')
[~,x2u]=getSamplingSchedules(6,7,0.5,0.6);
plot(xv,f(xv),'-k')
hold on
m2=plot(36*x2u,f(36*x2u),'.k','color',"#785EF0",'DisplayName','Non-uniform','MarkerSize',10)
xticks(0:12:36)
leg=legend([m1,m2],'NumColumns',2)
leg.Layout.Tile='south'



print(gcf,strcat('/home/turner/research/overleaf/guelph_talk/figs_csc/',fname) ...
    ,'-dpng','-r600') % -r sets the resolution