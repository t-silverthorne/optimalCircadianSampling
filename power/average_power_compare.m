addpath('utils/')
clear
simtype='verylong';
checkUseGPU
param.useGPU=true;
param.NL=4;
param.NR=4;
param.Nmeas=param.NL+param.NR;
param.freq_true=7.6;
param.Amp=1;
param.noise=.5;
param.Nacro=32;
nodes='uniform';

disp('hi')
tic
[acrovec,pwr]=wrap_simulatePWR_matperm_fv(param,nodes);
toc
close all
plot(acrovec,pwr)
%%
Navg=1e2;
p1_master = zeros(1,param.Nacro);
p2_master = zeros(1,param.Nacro);
for ii=1:Navg
    disp(ii)
    [acrovec,pwr_slow]=simulatePWR(param,nodes);
    [~,pwr_matperm_fv]=simulatePWR_matperm_fv(param,nodes);
    p1_master=p1_master+pwr_slow;
    p2_master=p2_master+pwr_matperm_fv;
end
p1_master=p1_master/Navg;
p2_master=p2_master/Navg;

close all
plot(acrovec,p1_master)
hold on
plot(acrovec,p2_master)