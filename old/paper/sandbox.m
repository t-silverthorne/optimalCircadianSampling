addpath('utils_core/')
addpath('utils_cost_fun/')

p.Nmeas     = 8;
p.Nacro     = 4;
p.Nresidual = 1e1;
p.Nperm     = 1e1;
p.noise     = 1;
p.Nbatch    = 1;
p.Amp       = 2;
p.freq      = 3.8;
p.permMethod='naive_reuse_perms';

[t,~]=getSamplingSchedules(p.Nmeas,0,0,0)


wrap_getCostFun(t,p,[4,5])