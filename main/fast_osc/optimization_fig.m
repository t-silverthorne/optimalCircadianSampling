close all
clear
addpath('../utils/')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


fast_run = false;
Nmeas=8;
Amp=2;
freq=3.8;

tiledlayout(1,2,'TileSpacing','tight','Padding','none')
%% plot power oscillations
ax_osc=nexttile(2);
[mt,~]   = getSamplingSchedules(Nmeas,0,0,0);
acros   = linspace(0,2*pi,2^7);
freqs=3.5:.1:4.4;
oldf=[freqs(1) freqs(end)];
c1rgb=[1.0, 0.01, 0.2]; % slow color
c4rgb=[.6 .4 .8]; % fast color
newf=freqs;
colors=[c1rgb; c4rgb];
cloc=interp1(oldf,colors,newf);

for ii=1:length(freqs)
    powvals = arrayfun(@(x) getPowerExact(Amp,x,freqs(ii),mt),acros);
    plot(acros,powvals,'-k','Color',cloc(ii,:))
    hold on
end


%% plot worst case powers 
ax_sweep=nexttile(1);
numFgrid=2^7;
numAgrid=2^4;
if fast_run
    numFgrid=2^5;
    numAgrid=2^3;
end
freqvals        = linspace(1,18,numFgrid);
Ampvals         = logspace(-1,1,numAgrid);

[Agr,Fgr]=meshgrid(Ampvals,freqvals);
Pmin=arrayfun(@(Agr,Fgr) getMinPower(Agr,Fgr,mt),Agr, Fgr);
[M,c]=contourf(Agr,Fgr,Pmin,100,'LineColor','none');%,'FaceAlpha',0.5)
hold on
yline(freqs(end),'--k','Color',c4rgb,'LineWidth',2)
yline(freqs(1),'--w','Color',c1rgb,'LineWidth',2)

nexttile(1)
set(gca,'XScale','log')
xlabel('amplitude')
ylabel('frequency')



nexttile(2)
ylim([0,1])
xlim([0,2*pi])
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$$0$$','$$\frac{\pi}{2}$$','$$\pi$$','$$\frac{3\pi}{2}$$','$$2\pi$$'})
ylabel('$\gamma(\phi)$','interpreter','latex')
xlabel('$\phi$','interpreter','latex')


clim manual

% tweak color bars
cb_osc = colorbar(ax_osc);
cb_osc.Layout.Tile='east';
colormap(ax_osc,cloc)
clim(ax_osc,[freqs(1) freqs(end)])
cb_osc.Label.Interpreter='latex';
set(cb_osc,'TickLabelInterpreter','latex')
cb_osc.Label.String='frequency';


colormap(ax_sweep,'parula')
cb = colorbar(ax_sweep);
clim(ax_sweep,[0,1])
cb.Label.Interpreter='latex';
set(cb,'TickLabelInterpreter','latex')
cb.Label.String='$\textrm{min}_\phi \gamma(\phi;A,f)$';
cb.Label.Interpreter='latex';


%% savefig
fontsize(gcf,10,"points")

plot_filename='optimization_setup';
ht=2.5; % height
wd=6.5; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,strcat('~/research/overleaf/rate_limited_sampling/figures/',plot_filename), ...
    '-dpng','-r600') % -r sets the resolution