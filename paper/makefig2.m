clear all; close all;
testing=false;

ny1=4*3;ny2=2*3;ny3=1*3;ny4=1*3;
nx1=2*3;nx2=2*3;
tileind=@(row,col) sub2ind([nx1+nx2,ny1+ny2+ny3+ny4],col,row);

tiledlayout(ny1+ny2+ny3+ny4,nx1+nx2,'TileSpacing','tight','Padding','tight')

if testing
    numFgrid=2;
    numAgrid=2;
else
    numFgrid=10;
    numAgrid=32;
end

Nmeasvals=[8,16,32];

for ii=1:3
    nexttile(tileind(1,1+(ii-1)*(nx1+nx2)/3),[ny1,(nx1+nx2)/3])
    p.Nmeas     = Nmeasvals(ii);
    if testing
        p.Nacro     = 4;
        p.Nresidual = 1e1;
        p.Nperm     = 1e1;
        p.noise     = 1;
        p.Nbatch    = 1;
    else
        p.Nacro     = 32;
        p.Nresidual = 1e3;
        p.Nperm     = 1e2;
        p.noise     = 1;
        p.Nbatch    = 1;
    end
    p.permMethod='fast';
    freqvals=linspace(1,32,numFgrid);
    Ampvals =logspace(-1,1,numAgrid);
    
    [t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);
    [I3,I4]=constructUtilMats(p);
    tic
    [Agr,Fgr]=meshgrid(Ampvals,freqvals);
    Pmin=arrayfun(@(Agr,Fgr) arrayfun_wrap_getPowerBatch(Agr,Fgr,t_unif,p,I3,I4), ...
                            Agr, Fgr);
    toc
    
    [M,c]=contourf(Agr,Fgr,Pmin,50,'LineColor','none');%,'FaceAlpha',0.5)
    set(gca,'XScale','log')
    drawnow
end


% aesthetics
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

nexttile(1)
ylabel('frequency')
yticks(0:8:32)
xticks([10^(-1), 10^0, 10^1])
xlabel('amplitude')
title('($N=8$)','Interpreter','latex')
clim manual
clim([0,1])

nexttile(1+(nx1+nx2)/3)
yticks(0:8:32)
yticklabels('')
xticks([10^(-1), 10^0, 10^1])
xlabel('amplitude')
title('($N=16$)','Interpreter','latex')
clim manual
clim([0,1])

nexttile(1+2*(nx1+nx2)/3)
yticks(0:8:32)
yticklabels('')
xticks([10^(-1), 10^0, 10^1])
xlabel('amplitude')
title('($N=32$)','Interpreter','latex')
clim manual
clim([0,1])

cb = colorbar;
cb.Layout.Tile = 'south';
cb.Label.Interpreter='latex'
cb.Label.String='$\textrm{min}_\phi \gamma(\phi)$';
set(cb,'TickLabelInterpreter','latex')
cb.Label.Interpreter='latex'

% plot_filename='figs/fig2'
% ht=3; % height
% wd=6; % width
% set(gcf,'PaperUnits','inches')
% set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
% print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
% savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too

%%

if testing
    p.Nmeas     = 8;
    p.Nacro     = 4;
    p.Nresidual = 1e1;
    p.Nperm     = 1e1;
    p.noise     = 1;
    p.Nbatch    = 1;
    p.Amp       = 2;
    p.freq      =3.8;
else
    p.Nmeas     = 8;
    p.Nacro     = 16;
    p.Nresidual = 1e3;
    p.Nperm     = 1e2;
    p.noise     = 1;
    p.Nbatch    = 10;
    p.Amp       = 2;
    p.freq      =3.8;
end
addpath('utils_core/')
addpath('utils_cost_fun/')
[I3,I4]=constructUtilMats(p);

% p.permMethod='fast';
% %% particle swarm
% opts = optimoptions(@particleswarm, ... %'HybridFcn',@fminsearch,...
%                                        'Display','iter', ...    
%                                        'SwarmSize',24, ...
%                                        'MaxIter',10, ... % TODO make realistic
%                                        'FunctionTolerance',1e-3, ... % might be too high
%                                        'UseParallel',true);%,'PlotFcn',@pswmyfun); % TODO make true
% 
% [t_opt,~]=particleswarm( @(t) -costfun_wrap_getPowerBatch(t,p,I3,I4),8,zeros(1,8),ones(1,8),opts)
% %% fmincon
eps_cstr=5e-3;
Aineq=eye(p.Nmeas-1,p.Nmeas); 

for ii=1:p.Nmeas-1
    Aineq(ii,ii+1)=-1;
end
bineq=-eps_cstr*ones(p.Nmeas-1,1);
[t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);
% [pwr,~,~]=getPowerBatch(t_unif,p,I3,I4)
% min(mean(pwr,1))
% 
% opts=optimoptions(@fmincon,'Display','iter','UseParallel',true,'MaxIter',10);
% x = fmincon(@(t) -costfun_wrap_getPowerBatch(t,p,I3,I4),t_unif,Aineq,bineq,[],[],zeros(1,8),ones(1,8),[],opts);
% 
% %% fminsearch
% opts=optimset('Display','iter','MaxIter',10);
% x=fminsearch(@(t) -fminsearch_wrap_getPowerBatch(t,p,I3,I4),log10(t_unif),opts);

% patternsearch
if testing
    opts=optimoptions(@patternsearch,'Display','iter','maxiter',2);
else
    opts=optimoptions(@patternsearch,'Display','iter','maxiter',50);
end
x = patternsearch(@(t) -costfun_wrap_getPowerBatch(t,p,I3,I4),t_unif,Aineq,bineq,[],[],zeros(1,8),ones(1,8),[],opts) ;

nexttile(tileind(ny1+1,1),[ny2,nx1+nx2])
[pwr,~,~]=getPowerBatch(t_unif,p,I3,I4);
pwr_unif=mean(pwr,1);
acrovec=linspace(0,2*pi,p.Nacro+1);
acrovec=acrovec(1:end-1);
plot(acrovec,pwr_unif)
hold on

[pwr,~,~]=getPowerBatch(x,p,I3,I4);
pwr_opt=mean(pwr,1);
plot(acrovec,pwr_opt)
legend({'uniform','optimal'},'Interpreter','latex','location','southoutside','NumColumns',2)
%% pareto plot
% just worry about power and acro variance
if testing
    popu_size=1;
    opts_pareto=optimoptions('paretosearch','Display','iter',...
                      'UseParallel',false, 'InitialPoints',repmat(t_unif,popu_size,1),...
                      'MeshTolerance',1e-2,'MaxIterations',3);
else
    popu_size=24;
    opts_pareto=optimoptions('paretosearch','Display','iter',...
                      'UseParallel',true, 'InitialPoints',repmat(t_unif,popu_size,1),...
                      'MeshTolerance',1e-2,'MaxIterations',3);
end
nexttile(tileind(ny1+ny2+1,1),[ny3+ny4,nx1+nx2])
Npar=5;
for ii=1:Npar
    [xopt_pareto,fopt_pareto] = paretosearch(@(t) wrap_getCostFun(t,p,I3,I4,[4,5]),p.Nmeas,Aineq,bineq,[],[],zeros(1,p.Nmeas),ones(1,p.Nmeas),[],opts_pareto);
    for ii=1:size(fopt_pareto,1)
        plot(fopt_pareto(ii,:),'.k')
        pause(0.2)
    end
    hold on
end
savefig('figs/fig2.fig')

% nexttile(tileind(ny1+ny2+1,1),[ny3,nx1])
% plot(1,1)
% nexttile(tileind(ny1+ny2+ny3+1,1),[ny3,nx1])
% plot(1,1)
% %% sweeps
% nexttile(tileind(ny1+ny2+1,nx1+1),[ny3,nx2])
% plot(1,1)
% nexttile(tileind(ny1+ny2+ny3+1,nx1+1),[ny3,nx2])
% plot(1,1)
% 

%% global search
% problem = createOptimProblem('fmincon',...
%     'objective',@(t) -costfun_wrap_getPowerBatch(t,p,I3,I4),...
%     'x0',t_unif,'options',...
%     optimoptions(@fmincon,'Display','iter','UseParallel',false,'MaxIter',10),);


function pmin=arrayfun_wrap_getPowerBatch(Amp,freq,t,p,I3,I4)
p.Amp=p.noise*Amp;
p.freq=freq;
[pwr,~,~]=getPowerBatch(t,p,I3,I4);
pmin=min(mean(pwr,1));
end

function pmin=costfun_wrap_getPowerBatch(t,p,I3,I4)
[pwr,~,~]=getPowerBatch(t,p,I3,I4);
pmin=min(mean(pwr,1));
end

function pmin=fminsearch_wrap_getPowerBatch(t,p,I3,I4)
t=sort(0.5*(1+tanh(10.^t)));
pmin=costfun_wrap_getPowerBatch(t,p,I3,I4);
end
