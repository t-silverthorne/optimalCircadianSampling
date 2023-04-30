addpath('utils_core/')
addpath('utils_cost_fun/')
tiledlayout(2,5,'TileSpacing','tight','Padding','tight')

p.permMethod='fast';

p.Nmeas     = 8;
p.Nacro     = 32;
p.Nresidual = 1e3;
p.Nperm     = 1e2;
p.noise     = 1;
p.Nbatch    = 2;
p.Amp       = 2;
p.freq      = 2.8;

[t,~]=getSamplingSchedules(p.Nmeas,0,0,0);
[I3,I4]=constructUtilMats(p);
pwr_record=NaN(i    i,1);
parfor ii=1:1e3
    [pwr,~,~]=getPowerBatch(t,p,I3,I4);
    pwr=mean(pwr,1);
    pwr_record(ii) = mean(pwr(1));
end
histogram(pwr_record)