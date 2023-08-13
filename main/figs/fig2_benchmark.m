close all
clear
addpath('../utils/')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

tiledlayout(1,2,'TileSpacing','tight')
nreps    = 5;
Amp      = 2;
fast_run = true;

nvals  = 5:15;
c1rgb  = [1.0, 0.44, 0.37]; % sparse color
c4rgb  = [.24   .17   .12];       % dense color
oldz   = [nvals(1) nvals(end)];
newz   = nvals;
colors = [c1rgb; c4rgb];
cloc   = interp1(oldz,colors,newz);


numFgrid=2^7; % for contour plots
numAgrid=2^4;
if fast_run
    numFgrid=2^5;
    numAgrid=2^3;
end
freqvals  = linspace(1,18,numFgrid);

for nn=1:length(nvals)
    Nmeas   = nvals(nn);
    [mt,~]  = getSamplingSchedules(Nmeas,0,0,0);


    powvals = NaN(1,numFgrid);
    lvals   = NaN(1,numFgrid);
    uvals   = NaN(1,numFgrid);
    tvals   = NaN(1,numFgrid);
    for ii=1:length(freqvals)
        fv=freqvals(ii);
    
        U=cell2mat(arrayfun( @(ind) directOptimize(Amp,fv,Nmeas), 1:nreps,'UniformOutput',false));
        
        powvals(ii) = mean(U(1,:));
        tvals(ii)   = mean(U(2,:));
        qv          = quantile(U(2,:),[.025 .975]);
        lvals(ii)   = qv(1);
        uvals(ii)   = qv(2);
    end
    
    
    sh(1)=nexttile(1);
    ee=errorbar(freqvals,1e3*tvals,lvals,uvals,'color',cloc(nn,:));
    
    hold on

    sh(2)=nexttile(2);
    plot(freqvals,powvals,'color',cloc(nn,:))
    hold on
    drawnow
end
%%
set(sh, 'Colormap', cloc , 'CLim', [nvals(1)-.5 nvals(end)+.5])
cb = colorbar(sh(end));
cb.Label.Interpreter='latex';
set(cb,'TickLabelInterpreter','latex')
cb.Label.String='$N$';
cbtop.Label.Interpreter='latex';

%%
nexttile(1)
ylabel('wall time (ms)')
xlabel('frequency $f$')

nexttile(2)
ylabel('$\textrm{min}_\phi \gamma(\phi;A,f)$')
xlabel('frequency $f$')
ylim([0,1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPORT for paper
plot_filename='/Users/turnersilverthorn/research/overleaf/samplingPaper/figures/fig2.png';
ht=2.5; % height
wd=6.5; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename, ...
    '-dpng','-r600') % -r sets the resolution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPORT for presentation
set(findall(gcf,'-property','FontSize'),'FontSize',7)
plot_filename='/Users/turnersilverthorn/research/overleaf/samplingTalk/figsAMMCS/fig2.png';
% Full slide dimensions
ht=2.9; % height
wd=4; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename, ...
    '-dpng','-r600') % -r sets the resolution



function power_time = directOptimize(Amp,freq,Nmeas)
options=optimoptions('fmincon',...
  'CheckGradients',false,...
  'SpecifyObjectiveGradient',true,'Display','none');

[mt,~]   = getSamplingSchedules(Nmeas,0,0,0);
tic
[x,~]    = fmincon(@(t) get_MINUS_MinLambdaAndDiff(t,freq),mt,[],[],[],[],zeros(1,Nmeas),ones(1,Nmeas),[],options);
time     = toc;
%Pow      = getMinPower(Amp,freq,x)/getMinPower(Amp,freq,mt);
Pow      = getMinPower(Amp,freq,x);
power_time = [Pow,time]';
end



