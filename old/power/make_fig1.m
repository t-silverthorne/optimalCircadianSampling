close all
clear

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

simtype='medium';          % how many Monte Carlo samples to include
addpath('utils/')
checkUseGPU              % uses simType to construct param
param.NL=4;              % samples in left interval
param.NR=4;              % samples in right interval
Nmeastot=param.NL+param.NR;
param.freq_true=1;     % freq used in regression model
param.Amp=1;             % signal amplitude
param.noise=1;          % the noise actually used in the simulation
param.Nacro=32;

tiledlayout(2,2,'TileSpacing','tight','Padding','tight')
param.Amp=2.5;
fpv=[1.5,2.25,4.25,8.25];
interp=(fpv-1.5)'/(8.25-1.5);
cvals=[0.7 .11 .11].*(1-interp)+ [0 0.53 0.74].*interp;
for ff=1:numel(fpv)
    param.freq_true=fpv(ff);
    [acrovec,pwr]=simulatePWR(param,'uniform');
    plot(acrovec,pwr,'-k','LineWidth',2,'Color',cvals(ff,:))
    ylim([0,1])
    hold on
    pause(0.5)
end
%%
nexttile(2)
[t,~]=getSamplingSchedules(param.NL,param.NR,0,0.25);
STDmat=NaN(3,numel(acrovec));
for aa=1:numel(acrovec)
    paramloc=param;
    paramloc.acro=acrovec(aa);
    X=constructX(t,paramloc);
    STDmat(:,aa)=sqrt(diag(inv(X'*X)))
end
close all
plot(acrovec,STDmat(3,:))
%%

xlim([0,2*pi])
xlabel('$\phi$')
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$$0$$','$$\frac{\pi}{2}$$','$$\pi$$','$$\frac{3\pi}{2}$$','$$2\pi$$'})
ylabel('$\beta(\phi)$')
legend({'$f=1.5$','$f=2.25$','$f=4.25$','$f=8.25$'},'Location','southoutside','NumColumns',4)

nexttile(3)
load('results/Numfreq_10_min_1_max_10_Numamp_8_min_1_max_4.mat')
[X,Y]=meshgrid(ampvals,freqvals);

cmap=[0 0.2 .13].*linspace(1,0,100)'+[0.94 .88 .19].*linspace(0,1,100)';
colormap(cmap)

[M,c]=contourf(X,Y,unif_mat,10,'LineColor','none');
c.LineWidth=2;
caxis([0 1]);
xlabel('$A$','interpreter','latex')
ylabel('$f$','Interpreter','latex')

h=colorbar;
h.Label.Interpreter='latex'
set(h,'TickLabelInterpreter','latex')
h.Label.String='$\mathrm{min}_\phi\beta(\phi;A,f)$';

nexttile(4)
[M,c]=contourf(X,Y,nu_mat-unif_mat,10,'LineColor','none');
c.LineWidth=2;
xlabel('$A$','interpreter','latex')
set(gca,'YTickLabel',[]);
caxis([0 1]);

h=colorbar;
h.Label.Interpreter='latex'
set(h,'TickLabelInterpreter','latex')
h.Label.String='$\Delta \;\mathrm{min}_\phi\beta(\phi;A,f)$';

%%
set(findall(gcf,'-property','FontSize'),'FontSize',10)
plot_filename='figs/fig1_revised'
ht=3.5; % height
wd=5; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too

function CV = get_CV_exact(acro,param,t)
paramloc=param;
paramloc.acro=acro;
beta=[0; paramloc.Amp*sin(acro); paramloc.Amp*cos(acro)];
X=constructX(t,paramloc);
CV=abs(sqrt(diag(inv(X'*X)))./beta);
end