% path and data
close all
tiledlayout(1,2)
fig=gcf

load('data/sweepNvals_test.mat')
addpath('utils_core')
addpath('utils_cost_fun')

% aesthetics 
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
fig_width = 10; % width in inches
fig_height = 5; % height in inches
fig_aspect_ratio = fig_width / fig_height;
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [fig_width fig_height], ...
         'PaperPosition', [0 0 fig_width fig_height], ...
         'PaperPositionMode', 'manual', ...
         'Units', 'inches', ...
         'Position', [0 0 fig_width fig_height]);

%%
clf

nexttile(1)
[M,c]=contourf(Agr,Fgr,Pmin,100,'LineColor','none');%,'FaceAlpha',0.5)
set(gca,'XScale','log')

clim manual
clim([0,1])

cb = colorbar;
cb.Layout.Tile = 'south';
cb.Label.Interpreter='latex'
cb.Label.String='$\textrm{min}_\phi \gamma(\phi)$';
set(cb,'TickLabelInterpreter','latex')
cb.Label.Interpreter='latex'

p.Nbatch=1;
p.Amp=3.5938;
p.freq=3 ;

testing=false;
if testing
    p.Nresidual = 1e3;
    p.Nperm     = 1e2;
    p.Nacro=8;
else
    p.Nresidual = 5e3;
    p.Nperm     = 2e2;
    p.Nacro=16;
    p.Nbatch=1;

end

xline(p.Amp,'-w')
yline(p.freq,'-w')

[~,indA]=min(abs(Agr(1,:)-p.Amp));
[~,indf]=min(abs(Fgr(:,1)-p.freq));


nexttile(2)
[pwr,~,~]=getPowerBatch(t_unif,p)
pwr=mean(pwr,1)
acrovec=linspace(0,2*pi,p.Nacro+1);
acrovec=acrovec(1:end-1); % get acrophases
plot(acrovec,reshape(pwr,[1 p.Nacro 1 1]),'-ok')
hold on
yline(Pmin(indf,indA),'--k')
drawnow
