clear all; close all;
testing=true;

ny1=4*3;ny2=2*3;ny3=1*3;ny4=1*3;
nx1=2*3;nx2=2*3;
tileind=@(row,col) sub2ind([nx1+nx2,ny1+ny2+ny3+ny4],col,row);

tiledlayout(1,3,'TileSpacing','tight','Padding','tight')

if testing
    numFgrid=2;
    numAgrid=2;
else
    numFgrid=10;
    numAgrid=32;
end

Nmeasvals=[8,16,32];

for ii=1:3
    nexttile(ii)
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

nexttile(2)
yticks(0:8:32)
yticklabels('')
xticks([10^(-1), 10^0, 10^1])
xlabel('amplitude')
title('($N=16$)','Interpreter','latex')
clim manual
clim([0,1])

nexttile(3)
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
savefig('figs/fig2top.fig')

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
