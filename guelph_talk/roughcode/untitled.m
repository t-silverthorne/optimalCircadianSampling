clear all
clf
addpath('../utils_core')
addpath('../utils_cost_fun')


c1="#648FFF"; % uniform colour
c2="#785EF0"; % 2-uniform colour
c3="#DC267F"; % Unif(0,1) colour
c4="#FE6100"; % jittered colour

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

font_size = 10;
set(0, 'DefaultAxesFontSize', font_size);
set(0, 'DefaultTextFontSize', font_size);



draw_lines('black','black','black','black')
yticks([])

plot_filename='waveforms_withuniform'
ht=2.6; % height
wd=5; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
fig=gcf;ax=fig.CurrentAxes;fig.Color='w';fig.OuterPosition=fig.InnerPosition;
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution

%%
clf
draw_lines(c1,'black','black','black')
yticks()
xlv=0:1/4:2; % uniform
xlv=xlv(1:end-1)
hx=xline(xlv,'-k','Color',c1,'LineWidth',2);

%%

clf
draw_lines('black','black','black','black')
xlv = xlv+.1*rand(1,8); % jittered
delete(hx)
hx=xline(xlv,'-k','Color',c4,'LineWidth',2);

%%
clf
draw_lines('black','black','black','black') % uniform
xlv = 2*rand(1,8); % random
delete(hx)
hx=xline(xlv,'-k','Color',c3,'LineWidth',2);


tauL=0.7; % two-uniform
tauR=0.8;
[~,xlv]=getSamplingSchedules(4,4,tauL,tauR);
xlv=xlv*2;
delete(hx)
hx=xline(xlv,'-k','Color',c2,'LineWidth',2);


%%

plot_filename='waveforms'
ht=2.6; % height
wd=5; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
fig=gcf;ax=fig.CurrentAxes;fig.Color='w';fig.OuterPosition=fig.InnerPosition;
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too



plot_filename='waveforms_withuniform'
ht=2.6; % height
wd=5; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
fig=gcf;ax=fig.CurrentAxes;fig.Color='w';fig.OuterPosition=fig.InnerPosition;
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution




plot_filename='waveforms_with2uniform'
ht=2.6; % height
wd=5; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
fig=gcf;ax=fig.CurrentAxes;fig.Color='w';fig.OuterPosition=fig.InnerPosition;
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too



function draw_lines(lc1,lc2,lc3,lc4)
t=0:.001:2;
y0=6;
y1=10;
y2=14;
y3=18;

plot(t,y3+cos(2*pi*t/2),'-k','Color',lc1)
hold on
plot(t,y2+cos(2*pi*t*4-0.1*pi),'-k','color',lc2)
plot(t,y1+sin(2*pi*t-0.4*pi)+.1*sin(2*pi*10*t),'-k','color',lc3)
plot(t,y0+sin(2*pi*t-pi)+(.02+1*exp(-(t-1.5).^2*1e2)).*sin(2*pi*30*t), ...
                          '-k','color',lc4)
end
