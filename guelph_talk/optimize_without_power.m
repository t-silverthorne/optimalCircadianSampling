addpath('utils_core')
addpath('utils_cost_fun')
numOpt=4;
nouter=100;
ninner=5e4;
%p.Nmeas     = 8;
%p.freq      = 3.8;
%p.Amp       = 2;
fname='data/ignore_power1.mat'

[t,~]=getSamplingSchedules(p.Nmeas,0,0,0);
p.Nresidual = 1e3;
p.Nperm     = 1;

p.noise     = 1;
p.Nacro     = 32; % choose random acrophase
p.Nbatch    = 1;
p.positive_cost_fun=true;
getCostFunNoPower(t,p)
%%
lambdaMaster=[];                 % Sample positive orthant of n sphere
while size(lambdaMaster,1)<nouter
    lambdaMat=randn(nouter,numOpt);
    lambdaMat=lambdaMat./sqrt(sum(lambdaMat.^2,2));
    lambdaMaster = [lambdaMaster; lambdaMat(sum(lambdaMat>0,2)==numOpt,:)];
end



eps_cstr    = 1e-3;
active_inds = 1:5;
tic
parfor ii=1:nouter
    opts = optimoptions(@simulannealbnd,'Display','none', ...
            'MaxIterations',ninner,'DisplayInterval',100,'ReannealInterval',100, ...
            'TemperatureFcn','temperaturefast');
    x0=rand(1,p.Nmeas);
    xopt = simulannealbnd(@(t) getScalarizedCfun(t,p,lambdaMaster(ii,:)), ...
			       x0,eps_cstr*ones(1,p.Nmeas),ones(1,p.Nmeas),opts);
    xmaster(ii,:)=xopt;
end
toc

save(fname)

function Jvec = getScalarizedCfun(t,p,lambda)
t=sort(t);
Jvec=max(getCostFunNoPower(t,p)./lambda);
end


function Jvec = getCostFunNoPower(t,p)
Y=getSimulatedData(t,p);
X=constructX(t,p);

betas_obs = pagemldivide(X'*X,pagemtimes(X',pagetranspose(Y))); % observed error
acro_est  = atan2(betas_obs(2,:,:,:),betas_obs(3,:,:,:));
amp_est   = sqrt(betas_obs(2,:,:,:).^2 + betas_obs(3,:,:,:).^2);

amp_est   = reshape(amp_est,[p.Nresidual,p.Nacro]);
acro_est  = reshape(acro_est,[p.Nresidual,p.Nacro]);

acrovec   = linspace(0,2*pi,p.Nacro+1);
acrovec   = acrovec(1:end-1);

[acro_bias,acro_var] = getAcroStats(acro_est,acrovec);
[amp_bias,amp_var]   = getAmpStats(amp_est,p.Amp);

Jvec(1)=max(abs(amp_bias)); 
Jvec(2)=max(abs(acro_bias));
Jvec(3)=max(amp_var);  
Jvec(4)=-min(acro_var);       

if p.positive_cost_fun % useful for scalarization
    Jvec(4)=Jvec(4)+1;
end
end
