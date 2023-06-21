addpath('utils_core/')
addpath('utils_cost_fun/')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
clf
p.Nmeas     = 30;
p.Nacro     = 1;
p.Nresidual = 1;
p.Nperm     = 1e2;
p.freq      = 5;
p.Amp       = 2.5;
p.noise     = 1;
p.Nbatch    = 1;

tiledlayout(3,1,'TileSpacing','tight')
Ncolour=20;
map=[.36 .54 .66].*(1-((1:Ncolour)'-1)/(Ncolour-1)) + [.81 .1 .26].*((1:Ncolour)'-1)/(Ncolour-1);
[t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);
[I3,I4]=constructUtilMats(p);
sc=4e-2
for ii=1:Ncolour
    for jj=1:3
        if jj==1
            t=t_unif;
        elseif jj==2
            t=sort(t_unif+sc*rand(1,p.Nmeas));
        else
            t=sort(rand(1,p.Nmeas));
        end
        Y=getSimulatedData(t,p);
        [pxx,f]=periodogram(Y,[],[],p.Nmeas);
        nexttile(jj)
        pl=semilogx(f,pxx,'-k','HandleVisibility','off');
        pl.Color=[pl.Color 0.4];
        hold on
        if ii==Ncolour
            xline(p.freq,'--k','DisplayName','true frequency')
        end
        hold on
        xlim([0,10])
    end
    drawnow   
end

nexttile(1)
legend('location','northwest')

% nvals_lab=5:5:Ncolour;
% h=colorbar;
% caxis([4 Ncolour])
% h.Limits=[4 Ncolour]
% h.Ticks=nvals_lab-.5;
% set(h,'TickLabelInterpreter','latex')
% h.TickLabels=string(nvals_lab);
% h.Label.String='$N_L$';
% h.Label.Interpreter='latex'
% colormap(map)
% h.Layout.Tile='east'

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
