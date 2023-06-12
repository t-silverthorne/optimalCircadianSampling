% path and data
close all
fig=gcf

load('data/sweepNvals.mat')
addpath('utils_core')
addpath('utils_cost_fun')
%%
% aesthetics 
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
fig_width = 5; % width in inches
fig_height = 3; % height in inches
fig_aspect_ratio = fig_width / fig_height;
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [fig_width fig_height], ...
         'PaperPosition', [0 0 fig_width fig_height], ...
         'PaperPositionMode', 'manual', ...
         'Units', 'inches', ...
         'Position', [0 0 fig_width fig_height]);


clf
tiledlayout(2,2,'TileSpacing','tight')
nexttile(1,[2 1])
[M,c]=contourf(Agr,Fgr,Pmin,50,'LineColor','none');%,'FaceAlpha',0.5)
set(gca,'XScale','log')
%%
clim manual
clim([0,1])

cb = colorbar;
cbc
cb.Label.Interpreter='latex'
cb.Label.String='$\textrm{min}_\phi \gamma(\phi;A,f)$';
set(cb,'TickLabelInterpreter','latex')
cb.Label.Interpreter='latex'

p.Nbatch = 1;
p.Amp    = 2;
p.freq   = 4;

testing=false;
if testing
    p.Nresidual = 1e3;
    p.Nperm     = 1e2;
    p.Nacro=8;
else
    p.Nresidual = 1e3;
    p.Nperm     = 1e2;
    p.Nacro=32;
    p.Nbatch=1;

end

xl=xline(p.Amp,'--k','LineWidth',2,'color','white');
yl=yline(p.freq,'--k','LineWidth',2,'color','white');

xl.Color(4)=1;
yl.Color(4)=1;

[~,indA]=min(abs(Agr(1,:)-p.Amp));
[~,indf]=min(abs(Fgr(:,1)-p.freq));


nexttile(2)
[pwr,~,~]=getPowerBatch(t_unif,p)
pwr=mean(pwr,1);
acrovec=linspace(0,2*pi,p.Nacro+1);
acrovec=acrovec(1:end-1); % get acrophases
plot(acrovec,reshape(pwr,[1 p.Nacro 1 1]),'-ok')
hold on
yline(Pmin(indf,indA),'--k')
drawnow


nexttile(1)
xlabel('amplitude')
ylabel('frequency')

nexttile(2)
xlabel('acrophase','Interpreter','latex')
ylabel('power','Interpreter','latex')

ylim([0,1])
xlim([0 2*pi])
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$$0$$','$$\frac{\pi}{2}$$','$$\pi$$','$$\frac{3\pi}{2}$$','$$2\pi$$'})
xline(pi/2,'linewidth',2,'color',"#648FFF")

%%
nexttile(4)
tv=0:.01:1;
oscf=@(t) p.Amp*cos(2*pi*p.freq*t-pi/2);
plot(tv,oscf(tv),'-k','LineWidth',1)
hold on
plot(t_unif,oscf(t_unif),'.k','MarkerSize',20,'Color',"#648FFF")
xlabel('$t$')
ylabel('$\cos(2\pi f t - \pi/2)$','Interpreter','latex')


%%
print(gcf,strcat('/home/turner/research/overleaf/guelph_talk/figs_csc/','sweepNvals.png'),'-dpng','-r600') % -r sets the resolution
savefig(fname)


