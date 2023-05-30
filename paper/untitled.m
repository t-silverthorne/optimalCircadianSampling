% test of wrapper for cost function
addpath('utils_core')
addpath('utils_cost_fun')
p.Nmeas     = 8;
p.Nacro     = 16;
p.Nresidual = 2e3;
p.Nperm     = 200;
p.freq      = 3.8;
p.Amp       = 1.5;
p.noise     = 1.5;
p.Nbatch    = 2;

p.permMethod = 'FY_double_for_loop';
p.permActionMethod='index'
p.perm_start_ind = 1;
p.perm_stop_ind  = factorial(p.Nmeas);
[t,~]=getSamplingSchedules(p.Nmeas,0,0,0);

tic
[pwr,~,~]=getPowerBatch(t,p);
toc

clf
if p.Nbatch>1
    plot(mean(pwr),'-o')
else
    plot(pwr,'-o')
end