% test of wrapper for cost function
p.Nmeas     = 10;
p.Nacro     = 16;
p.Nresidual = 5e2;
p.Nperm     = 2e2;
p.freq      = 2;
p.Amp       = 1;
p.noise     = 1;
p.Nbatch    = 1;
addpath('../utils_core')
[t,~]=getSamplingSchedules(p.Nmeas,0,0,0);
[X,I3,I4]=constructUtilMats(t,p);
tic
wrap_getCostFun(t,p,X,I3,I4)
toc
%%
tic
getPowerBatch(t,p,X,I3,I4);
toc