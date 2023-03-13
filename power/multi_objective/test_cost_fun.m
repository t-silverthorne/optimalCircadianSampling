addpath('utils/')
simtype='medium';
checkUseGPU
param.useGPU=true;
param.NL=4;
param.NR=4;
param.Nmeas=param.NL+param.NR;
param.freq_true=2.4;
param.Amp=1;
param.noise=.6;
param.Nacro=32;
nodes='uniform';

v=0;
nrun=10;
for i=1:nrun
    disp(i)
    v=v+costfun_power_bias_var(param,nodes);
end
v/nrun

% 0.0647    0.0233    0.3136    0.9441    0.4969

g%%

%costfun = @(t) out2_wrap_simulatePWR_matperm_fv(param,t);

costfun = @(t) costfun_power_bias_var(param,t);

costfun('uniform')
%%