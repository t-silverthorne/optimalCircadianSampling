% test of wrapper for cost function
p.Nmeas     = 8;
p.Nacro     = 16;
p.Nresidual = 1e3;
p.Nperm     = 1e1;
p.freq      = 7.9;
p.Amp       = 2.5;
p.noise     = 1;
p.Nbatch    = 1;
addpath('../utils_core')

clf
[t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);
[I3,I4]=constructUtilMats(p);
N=54;
sc=5e-2;    
wrap_getCostFun(t_unif,p,I3,I4)