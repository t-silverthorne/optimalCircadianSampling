clf
clear
addpath('../utils/')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% aesthetic changes
tiledlayout(2,2,'TileSpacing','tight','Padding','none')
freqs=3.5:.1:4.4;           
oldf=[freqs(1) freqs(end)];
c1rgb=[1.0, 0.01, 0.2]; % slow color
c4rgb=[.6 .4 .8];       % fast color
newf=freqs;
colors=[c1rgb; c4rgb];
cloc=interp1(oldf,colors,newf);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET PARAMETERS AND CONSTANTS 
Nmeas    = 8;
Amp      = 2;
freq     = 3.8;
[mt,~]   = getSamplingSchedules(Nmeas,0,0,0);
fast_run = false;

numFgrid=2^7; % for contour plots
numAgrid=2^4;
if fast_run
    numFgrid=2^5;
    numAgrid=2^3;
end

freq_for_opt = 3.8;
maxit        = 1e3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOP LEFT 
cont(1)=nexttile(1);
freqvals  = linspace(1,18,numFgrid);
Ampvals   = logspace(-1,1,numAgrid);

[Agr,Fgr] = meshgrid(Ampvals,freqvals);
Pmin      = arrayfun(@(Agr,Fgr) getMinPower(Agr,Fgr,mt),Agr, Fgr);
[M,c]     = contourf(Agr,Fgr,Pmin,100,'LineColor','none');%,'FaceAlpha',0.5)

% aesthetics and other lines
hold on % lines for indicating cross sections
yline(freqs(1),'--w','Color',c1rgb,'LineWidth',2)
yline(freqs(end),'--k','Color',c4rgb,'LineWidth',2)
ylim([freqvals(1), freqvals(end)])
set(gca,'XScale','log') % rescale x axis
xlabel('amplitude') % xlabel
ylabel('frequency') % ylabel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOP RIGHT 
cont(2) = nexttile(2);

% run simulated annealing to maximize min power 
opts = optimoptions(@simulannealbnd,'Display','iter', ...
    'MaxIterations',maxit,'DisplayInterval',100,'ReannealInterval',50);
x0=rand(1,Nmeas); % initial guess
xopt = simulannealbnd(@(t) -getMinPower(Amp,freq_for_opt,t), ...
                   x0,zeros(1,Nmeas),ones(1,Nmeas),opts);

% make contour plot
[Agr,Fgr]=meshgrid(Ampvals,freqvals);
Pmin=arrayfun(@(Agr,Fgr) getMinPower(Agr,Fgr,xopt),Agr, Fgr);
[M,c]=contourf(Agr,Fgr,Pmin,100,'LineColor','none');%,'FaceAlpha',0.
hold on % lines for indicating cross sections
yline(freqs(1),'--w','Color',c1rgb,'LineWidth',2)
yline(freqs(end),'--k','Color',c4rgb,'LineWidth',2)
ylim([freqvals(1), freqvals(end)])
set(gca,'XScale','log') % rescale x axis
xlabel('amplitude') % xlabel
set(gca,'YTickLabel',[]);


% add colorbar shared by first two plots
set(cont, 'Colormap',parula, 'CLim', [0 1])
% assign color bar to one tile 
cbtop = colorbar(cont(end));
cbtop.Label.Interpreter='latex';
set(cbtop,'TickLabelInterpreter','latex')
cbtop.Label.String='$\textrm{min}_\phi \gamma(\phi;A,f)$';
cbtop.Label.Interpreter='latex';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOTTOM LEFT 
osc_row(1)=nexttile(3);
acros   = linspace(0,2*pi,2^7);
for ii=1:length(freqs)
    powvals = arrayfun(@(x) getPowerExact(Amp,x,freqs(ii),mt),acros);
    plot(acros,powvals,'-k','Color',cloc(ii,:))
    hold on
end
xlabel('$\phi$','interpreter','latex') % xlabel
ylabel('$\gamma(\phi)$','interpreter','latex') % ylabel
ylim([0,1])
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$$0$$','$$\frac{\pi}{2}$$','$$\pi$$','$$\frac{3\pi}{2}$$','$$2\pi$$'})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOTTOM RIGHT 
osc_row(2)=nexttile(4);
for ii=1:length(freqs)
    powvals = arrayfun(@(x) getPowerExact(Amp,xopt,freqs(ii),mt),acros);
    plot(acros,powvals,'-k','Color',cloc(ii,:))
    hold on
end

ylim([0,1])
xlim([0,2*pi])
set(gca,'YTickLabel',[]);
xlabel('$\phi$','interpreter','latex') % xlabel
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$$0$$','$$\frac{\pi}{2}$$','$$\pi$$','$$\frac{3\pi}{2}$$','$$2\pi$$'})

% add colorbar shared by bottom two plots
set(osc_row, 'Colormap', cloc , 'CLim', [freqs(1) freqs(end)])
cb_osc = colorbar(osc_row(end));
cb_osc.Label.Interpreter='latex';
set(cb_osc,'TickLabelInterpreter','latex')
cb_osc.Label.String='frequency';
cbtop.Label.Interpreter='latex';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPORT 
plot_filename='single_freq_optimization_combined';
ht=2.5*2; % height
wd=6.5; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,strcat('~/research/overleaf/rate_limited_sampling/figures/',plot_filename), ...
    '-dpng','-r600') % -r sets the resolution
