%% setup parameters
testing=false;
fname = 'test_amp4';
addpath('utils_core')
addpath('utils_cost_fun/')

numOpt=5;
p.positive_cost_fun=true;

if testing
    p.Nresidual = 50;
    p.Nperm     = 10;
    p.Nacro     = 8; % choose random acrophase
    nouter      = 5; % number of outer optimzation loops
    ninner      = 10; % number of inner optimization loops
    parpool_size=2;
else
    p.Nresidual = 1e3;
    p.Nperm     = 1e2;
    p.Nacro     = 16; % choose random acrophase
    nouter      = 30; % number of outer optimzation loops
    ninner      = 1000; % number of inner optimization loops
%    ninner      = 100; % number of inner optimization loops
    parpool_size=30;
end
parpool(parpool_size)

p.Nmeas     = 8;  % eventually loop over this
p.freq      = 3.8;
p.Amp       = 2;
p.noise     = 1;
p.Nbatch    = 1;

p.permMethod       = 'FY_double_for_loop'; %p.permActionMethod = 'index'; % options index or matrix for 'naive_make_perms_first'
p.permActionMethod = 'index';

%% get scalarization samples


lambdaMaster=[];                 % Sample positive orthant of n sphere
while size(lambdaMaster,1)<nouter
    lambdaMat=randn(nouter,numOpt);
    lambdaMat=lambdaMat./sqrt(sum(lambdaMat.^2,2));
    lambdaMaster = [lambdaMaster; lambdaMat(sum(lambdaMat>0,2)==numOpt,:)];
end

%% run optimization
active_inds = 1:5;
parfor ii=1:nouter
    opts = optimoptions(@simulannealbnd,'Display','iter', ...
            'MaxIterations',ninner,'DisplayInterval',1,'ReannealInterval',50,'TemperatureFcn','temperaturefast');
    x0=rand(1,p.Nmeas);
    xopt = simulannealbnd(@(t) hvol_scalar_costfun_v2(t,p,active_inds,lambdaMaster(ii,:)), ...
			       x0,zeros(1,p.Nmeas),ones(1,p.Nmeas),opts);
    xmaster(ii,:)=xopt;
end



save(fname)


function J=hvol_scalar_costfun_v2(t,p,active_inds,lambdaloc)
t=sort(t);
J=max(wrap_getCostFun(t,p,active_inds)./lambdaloc);
end










