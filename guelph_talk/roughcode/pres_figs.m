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


rng(2345)
p.Nmeas=8;
N=p.Nmeas;
p.permMethod='naive_reuse_perms';
tauL=0.4;
tauR=0.6;
[t_unif,t_nu]=getSamplingSchedules(4,4,tauL,tauR);
t_rand=sort(rand(1,N));

t_jit = t_unif + 1e-2*rand(1,N);
c1="#648FFF"; % uniform colour
c2="#785EF0"; % 2-uniform colour
c3="#DC267F"; % Unif(0,1) colour
c4="#FE6100"; % jittered colour

plot(t_unif,4,'ok','MarkerEdgeColor','none','MarkerFaceColor',c1)

hold on
plot(t_nu,3,'ok','MarkerEdgeColor','none','MarkerFaceColor',c2)
plot(t_rand,2,'ok','MarkerEdgeColor','none','MarkerFaceColor',c3)
plot(t_jit,1,'ok','MarkerEdgeColor','none','MarkerFaceColor',c4)

yv=.2*linspace(-1,1,10);

plot(repmat(tauL,1,10),3+yv,'-k','Color',c2)
plot(repmat(tauR,1,10),3+yv,'-k','Color',c2)
ylim([0.5 4.5])
xlabel('$t$','interpreter','latex')
yticks([1 2 3 4])
yticklabels({'jittered','random','2-uniform','uniform'})

p.Nacro     = 64;
p.freq      = 3.8;
p.Amp       = 1.5;
p.noise     = .5;
p.Nbatch    = 20;
p.Nresidual = 1;
p.Nperm     = 1e2;

t_hires=linspace(0,1,100);
[t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);



plot_filename='sample_grids'
ht=2.6; % height
wd=5; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
fig=gcf;ax=fig.CurrentAxes;fig.Color='w';fig.OuterPosition=fig.InnerPosition;
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too



%%