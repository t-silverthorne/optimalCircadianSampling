testing=false;
fname='sweepNvals_largerN';
addpath('utils_core')
addpath('utils_cost_fun/')
tiledlayout(2,4,'TileSpacing','tight','Padding','tight')

if testing
    numFgrid=2;
    numAgrid=2;
    Nmeasvals=8:15; 
    partname='test';
else
    numFgrid=32*4;
    numAgrid=20;
end

if testing
    p.Nresidual = 50;
    p.Nperm     = 10;
    p.Nacro     = 8; % choose random acrophase
else
    p.Nresidual = 5e3;
    p.Nperm     = 2e2;
    p.Nacro     = 16; % choose random acrophase
end
p.Nmeas     = 8;  % eventually loop over this

p.freq      = 3.7;
p.Amp       = 2;
p.noise     = 1;
p.Nbatch    = 1;


p.permMethod       = 'FY_double_for_loop'; %p.permActionMethod = 'index'; % options index or matrix for 'naive_make_perms_first'
p.permActionMethod = 'index';
freqvals        = linspace(1,18,numFgrid);
Ampvals         = logspace(-1,1,numAgrid);

[t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);

tic
[Agr,Fgr]=meshgrid(Ampvals,freqvals);
Pmin=arrayfun(@(Agr,Fgr) arrayfun_wrap_getPowerBatch(Agr,Fgr,t_unif,p), ...
                        Agr, Fgr);
toc

save(fname)
%[M,c]=contourf(Agr,Fgr,Pmin,50,'LineColor','none');%,'FaceAlpha',0.5)
%set(gca,'XScale','log')
%drawnow


function pmin=arrayfun_wrap_getPowerBatch(Amp,freq,t,p)
p.Amp=Amp;
p.freq=freq;
[pwr,~,~]=getPowerBatch(t,p);
pmin=min(mean(pwr,1));
end

