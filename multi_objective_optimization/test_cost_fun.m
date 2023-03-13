clear
poolobj = gcp('nocreate');
delete(poolobj)
gpuDevice(1)
addpath('utils/')
simtype='long';
checkUseGPU
param.useGPU=true;
param.useParallel=false;
if param.useParallel
    param.poolobj=parpool;
end
param.NL=3;
param.NR=3;
param.Nmeas=param.NL+param.NR;
param.freq_true=7.6;
param.Amp=1;
param.noise=.5;
param.Nacro=16;
nodes='uniform';
Nmeastot=param.NL+param.NR;
[t_unif,~]=getSamplingSchedules(param.NL,param.NR,0,0.5); % initial guess for sampling

param.Pmat=construct_all_perms(param.Nmeas);

% multi objective func to be minimized
costfun = @(t) costfun_power_bias_var(param,t);

%costfun(rand(1,Nmeastot))-costfun(t_unif)

[]=simulatePWR_matperm_fv(param,'uniform')
%% Run constrained optimization


Aineq=eye(Nmeastot-1,Nmeastot); % inequality constraints
for ii=1:Nmeastot-1
    Aineq(ii,ii+1)=-1;
end
bineq=ones(Nmeastot-1,1);

%% paretosearch (pattern search for multiobjective

opts=optimoptions('paretosearch','Display','iter',...
                  'InitialPoints',t_unif,'UseParallel',param.useParallel, ...
                  'MeshTolerance',1e-2 );
[xopt,fopt] = paretosearch(costfun,Nmeastot,Aineq,bineq,[],[],zeros(Nmeastot,1),ones(Nmeastot,1),[],opts);

%% gamultiobj (genetic algorithm for multiobjective)


param.useParallel=true;
Nmeastot=param.NL+param.NR;
opts=optimoptions('gamultiobj','Display','iter','InitialPopulationMatrix',t_unif,'UseParallel',param.useParallel,...
                      'PopulationSize',50);
[xopt,fopt] = gamultiobj(costfun,Nmeastot,Aineq,bineq,[],[],zeros(Nmeastot,1),ones(Nmeastot,1),[],opts);

%% Run reparameterized optimization

%% Penalty methods

