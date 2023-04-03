addpath('../utils/')

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

clf
param.NL=4;
param.NR=4;
param.freq_true=3;
param.Amp=2.5;
param.noise=1;
param.acrophase=pi/2;

tiledlayout(2,1,'TileSpacing','tight')
N=20;
map=[.36 .54 .66].*(1-((1:N)'-1)/(N-1)) + [.81 .1 .26].*((1:N)'-1)/(N-1);
jitter_scale=0.05;

for ii=1:N
    [tu,~]=getSamplingSchedules(param.NL,param.NR,0,0.5);
    tnu=tu+jitter_scale*(rand(1,numel(tu))-0.5);
    for jj=1:2
        if jj==1
            t=tu;
        else
            t=tnu;
        end
        Y=param.Amp*cos(2*pi*t*param.freq_true-param.acrophase)+param.noise*randn(1,numel(t));
        [pxx,f]=periodogram(Y,[],[],param.NL+param.NR);
        nexttile(jj)
        semilogx(f,pxx,'-k','color',map(ii,:),'HandleVisibility','off')
        hold on
        if ii==N
            xline(param.freq_true,'--k','DisplayName','true frequency')
        end
        hold on
        xlim([0,10])
    end
    drawnow   
end
%%
nexttile(1)
legend('location','northwest')

nvals_lab=5:5:N;
h=colorbar;
caxis([4 N])
h.Limits=[4 N]
h.Ticks=nvals_lab-.5;
set(h,'TickLabelInterpreter','latex')
h.TickLabels=string(nvals_lab);
h.Label.String='$N_L$';
h.Label.Interpreter='latex'
colormap(map)
h.Layout.Tile='east'

nexttile(1)
ylabel('psd')
set(gca,'XTickLabel',[]);

nexttile(2)
xlabel('frequency')
ylabel('psd')

nexttile(1)
title('($N_R=4,\tau=0.5$)')

plot_filename='biased_spectrum_example';
ht=4; % height
wd=6; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too
