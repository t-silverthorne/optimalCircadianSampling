% path and data
close all
fig=gcf

load('data/sweepNvals.mat')
addpath('utils_core')
addpath('utils_cost_fun')
def_colours
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
tiledlayout(1,2,'TileSpacing','loose')
nexttile(2)
[M,c]=contourf(Agr,Fgr,Pmin,50,'LineColor','none');%,'FaceAlpha',0.5)
set(gca,'XScale','log')
xlabel('amplitude')
ylabel('frequency/samling rate')
yticks([0 4 8 12 16])
yticklabels({'0','1/2','1','3/2','2'})

clim manual
clim([0,1])

cb = colorbar;
cb.Label.Interpreter='latex';
%cb.Label.String='$\textrm{min}_\phi \gamma(\phi;A,f)$';
cb.Label.String='worst-case power';
set(cb,'TickLabelInterpreter','latex')
cb.Label.Interpreter='latex';



p.Nbatch = 1;
p.Amp    = 2;
p.freq   = 4;

nexttile(2)
xl=xline(p.Amp,'--k','LineWidth',2,'color','white');
yl=yline(p.freq,'--k','LineWidth',2,'color','white');
xl.Color(4)=1;
yl.Color(4)=1;

[~,indA]=min(abs(Agr(1,:)-p.Amp));
[~,indf]=min(abs(Fgr(:,1)-p.freq));

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

nexttile(1)
tv=0:.01:1;
oscf=@(t) p.Amp*cos(2*pi*p.freq*t-pi/2);
plot(tv,oscf(tv),'-k','LineWidth',1)
hold on
plot(t_unif,oscf(t_unif),'.k','MarkerSize',20,'Color',c1)
xlabel('$t$ (hr)')
ylabel('$\cos(2\pi f t - \pi/2)$','Interpreter','latex')
xticks([0 0.5 1])
xticklabels({'0','12','24'})


print(gcf,strcat('/home/turner/research/overleaf/guelph_talk/figs_csc/','sweepNvals.png'),'-dpng','-r600') % -r sets the resolution
savefig(fname)


