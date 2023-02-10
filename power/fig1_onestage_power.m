close all
clear
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

freqvals=linspace(0.5,16,160);
SNRvals=logspace(-1,1,10);
cvals=[0 0 1].*linspace(0,1,numel(freqvals))' + [1 0 0].*(1-linspace(0,1,numel(freqvals))');;

% for ff=1:numel(freqvals)
%     param.freq_true=freqvals(ff);
%     pwr_vals=NaN(1,numel(SNRvals));
%     for ss=1:numel(SNRvals)
%         param.Amp=param.noise*SNRvals(ss);
%         [~,pwr]=simulatePWR(param,'uniform');
%         pwr_vals(ss)=min(pwr);
%     end
%     plot3(SNRvals,ones(1,numel(pwr_vals))*freqvals(ff),pwr_vals,'-k')
%     set(gca,'XScale','log')
%     set(gca,'YScale','log')
%     hold on
%     zlim([0,1])
%     drawnow
% end
tic
[S,F]=meshgrid(SNRvals,freqvals);
P=arrayfun(@(S,F) arrayfun_wrap_simulatePWR(S,F,param), S, F);
toc


%%
close all
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%%
close all
tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
nexttile(1)

cmap=[0 0.2 .13].*linspace(1,0,100)'+[0.94 .88 .19].*linspace(0,1,100)'
colormap(cmap)
[M,c]=contourf(S,F,P,50,'LineColor','none');%,'FaceAlpha',0.5)
set(gca,'XScale','log')
set(gca,'YScale','log')
c.LineWidth=2;

% fs=[1 2 4.5 8.5]';
% col_interp=(fs-1)/(fs(end)-1);
% col1=[0.9 0.17 0.31];
% col2=[0.6 .4 .8];
% cvals=col1.*col_interp+col2.*(1-col_interp);
% for ii=1:numel(fs)
%     plot3(0.1*ones(1,10),fs(ii)*ones(1,10),linspace(0,1,10),'-k','color',cvals(ii,:),'LineWidth',4)
% end
yticks([1 2 4 8 16])
xlabel('$\textrm{SNR}$','interpreter','latex')
ylabel('$f$','Interpreter','latex')
zlabel('$\beta$','Interpreter','latex')

h=colorbar

h.Ticks=0:.1:1
set(h,'TickLabelInterpreter','latex')
h.Label.Interpreter='latex'
h.Label.String='$\mathrm{min}_\phi\beta(\phi;A,f)$';

xline(1.5,'--w')
hold on
yline(1.5,'--w')
yline(2.25,'--w')
yline(4.25,'--w')
yline(8.25,'--w')


nexttile(2)
param.Amp=2.5;
fpv=[1.5,2.25,4.25,8.25];
interp=(fpv-1.5)'/(8.25-1.5);
cvals=[0.7 .11 .11].*(1-interp)+ [0 0.53 0.74].*interp 

for ff=1:numel(fpv)
    param.freq_true=fpv(ff);
    [acrovec,pwr]=simulatePWR(param,'uniform');
    plot(acrovec,pwr,'-k','LineWidth',2,'Color',cvals(ff,:))
    ylim([0,1])
    hold on
    pause(0.5)
end

xlim([0,2*pi])
xlabel('$\phi$')
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$$0$$','$$\frac{\pi}{2}$$','$$\pi$$','$$\frac{3\pi}{2}$$','$$2\pi$$'})
ylabel('$\beta(\phi)$')
legend({'$f=1.5$','$f=2.25$','$f=4.25$','$f=8.25$'},'Location','southoutside','NumColumns',4)
% cmap=cvals(1,:).*linspace(1,0,100)'+cvals(end,:).*linspace(0,1,100)'
% 
% colormap(cmap)
% h=colorbar
% caxis([min(fpv) max(fpv)]) % use this 
% clim([min(fpv) max(fpv)])
% h.Ticks=fpv
% set(h,'TickLabelInterpreter','latex')
% h.Label.Interpreter='latex'
% h.Label.String='$f$';

% add blue dots at intersections
%%
set(findall(gcf,'-property','FontSize'),'FontSize',10)
plot_filename='figs/fig1a'
ht=2.5; % height
wd=7; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too
function pout=arrayfun_wrap_simulatePWR(SNR,freq,param)
param.Amp=param.noise*SNR;
param.freq_true=freq;
[~,pwr]=simulatePWR(param,'uniform');
pout=min(pwr);
end