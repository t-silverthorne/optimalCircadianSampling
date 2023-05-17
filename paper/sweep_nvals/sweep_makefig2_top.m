testing=false;
addpath('../utils_core/')
addpath('../utils_cost_fun/')
tiledlayout(2,4,'TileSpacing','tight','Padding','tight')

if testing
    numFgrid=2;
    numAgrid=2;
    Nmeasvals=8:15; 
    partname='test'
else
    numFgrid=32;
    numAgrid=10;
end

for ii=1:length(Nmeasvals)
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
        p.Nperm     = 2e2;
        p.noise     = 1;
        p.Nbatch    = 1;
    end
    p.permMethod       = 'naive_reuse_perms'; %p.permActionMethod = 'index'; % options index or matrix for 'naive_make_perms_first'
    freqvals=linspace(1,32,numFgrid);
    Ampvals =logspace(-1,1,numAgrid);
    
    [t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);
    tic
    [Agr,Fgr]=meshgrid(Ampvals,freqvals);
    Pmin=arrayfun(@(Agr,Fgr) arrayfun_wrap_getPowerBatch(Agr,Fgr,t_unif,p), ...
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

for ii=1:length(Nmeasvals)
    nexttile(ii)
    yticks(0:8:32)
    xticks([10^(-1), 10^0, 10^1])
    xlabel('amplitude')
    if (ii==1 | ii==5)
        ylabel('frequency')
    else
        yticklabels('')
    end

    title(strcat('($N=',int2str(Nmeasvals(ii)),'$)'),'Interpreter','latex')
    clim manual
    clim([0,1])
end
cb = colorbar;
cb.Layout.Tile = 'south';
cb.Label.Interpreter='latex';
cb.Label.String='$\textrm{min}_\phi \gamma(\phi)$';
set(cb,'TickLabelInterpreter','latex')
cb.Label.Interpreter='latex';
savefig(strcat('fig2topsweep_',string(partname),'.fig'))


function pmin=arrayfun_wrap_getPowerBatch(Amp,freq,t,p)
p.Amp=p.noise*Amp;
p.freq=freq;
[pwr,~,~]=getPowerBatch(t,p);
pmin=min(mean(pwr,1));
end

function pmin=costfun_wrap_getPowerBatch(t,p)
[pwr,~,~]=getPowerBatch(t,p);
pmin=min(mean(pwr,1));
end
