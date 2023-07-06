clear all
testing=false;
addpath('utils_core')
addpath('utils_cost_fun/')

numOpt=5;
p.positive_cost_fun=true;

p.Nmeas     = 5;  % eventually loop over this
if testing
    p.Nresidual = 50;
    p.Nperm     = 10;
    p.Nacro     = 8; % choose random acrophase
    nouter      = 5; % number of outer optimzation loops
    ninner      = 10; % number of inner optimization loops
    parpool_size=2;
else
    p.Nresidual = 1e3;
    p.Nperm     = factorial(p.Nmeas);
    p.Nacro     = 5; % choose random acrophase
end
nMC=4;

p.freq      = 1.3;
p.Amp       = 4.4;
p.noise     = 1;
p.Nbatch    = 1;

%p.permMethod       = ''; %p.permActionMethod = 'index'; % options index or matrix for 'naive_make_perms_first'
p.perm_start_ind=1;
p.perm_stop_ind=p.Nperm;
p.permMethod = 'exact_in_batches';%'index';

[t,~]=getSamplingSchedules(p.Nmeas,0,0,0);
close all
pwr=0;
for ii=1:nMC
    pwr=getPower(t,p);

    %pwr=pwr/nMC;
    pwr=reshape(pwr,1,length(pwr));
    ac=linspace(0,2*pi,p.Nacro+1);
    ac=ac(1:end-1);
    
    
    plot(ac,pwr,'-ok')
    hold on
    ylim([0,1])
    drawnow
end
pwr_exact = arrayfun(@(phi) getPowerExact(p.Amp,phi,p.freq,t),ac)
hold on
plot(ac,pwr_exact,'-ob')
function beta=getPowerExact(Amp,acro,freq,mt)
Nmeas    = length(mt); % num samples
alpha    = 0.05; % type I error rate

csq    = cos(2*pi*freq*mt-acro);
lambda =  Amp^2*csq*csq'; % non-centrality

beta=1-ncfcdf( finv(1-alpha,2,Nmeas-3) ,2,Nmeas-3,lambda);
end