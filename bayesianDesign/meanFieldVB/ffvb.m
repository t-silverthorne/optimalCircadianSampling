clf
clear all
close all
rng(1)
%% simulated data
A0=2.5;
B0=2.5;
f0=2.1;
fmax=12;
Lcut=fmax;
T=1000;
Nsamp=1e4;
Nmeas=8;

t_obs=linspace(0,1,Nmeas+1);
t_obs=t_obs(1:end-1);
y_obs=A0*cos(2*pi*f0*t_obs)+B0*sin(2*pi*f0*t_obs) + randn(1,numel(t_obs));

%% priors
settings.muAprior  =0;
settings.sigAprior =2;
settings.muBprior  =0;
settings.sigBprior =2;
settings.muTprior  =0;
settings.sigTprior =1;
settings.fmax=fmax;
settings.y_obs=y_obs;
settings.t_obs=t_obs;

muAprior  =  settings.muAprior;
sigAprior =  settings.sigAprior;
muBprior  =  settings.muBprior;
sigBprior =  settings.sigBprior;
muTprior  =  settings.muTprior;
sigTprior =  settings.sigTprior;

%% MCMC
tic
X=randn(Nsamp,3);
X(:,1)=X(:,1)*sigAprior+muAprior; X(:,2)=X(:,2)*sigBprior+muBprior;
X(:,1)=X(:,1)*sigTprior+muTprior;

pfun=@(par) -sum((y_obs - par(:,1).*cos(2*pi*0.5*(fmax+fmax*tanh(par(:,3)))*t_obs) ... % likelihood
    - par(:,2).*sin(2*pi*0.5*(fmax +fmax*tanh(par(:,3)))*t_obs)  ).^2/2,2) - ...
    log(sqrt(2*pi)) - ... % normalization of likelihood (sigma=1)
    (muAprior-par(:,1)).^2/2/sigAprior^2 -  ... % prior
    (muBprior-par(:,2)).^2/2/sigBprior^2 - ...  
    (muTprior-par(:,3)).^2/2/sigTprior^2;

Y=randn(Nsamp,3,T);
Y(:,1)=Y(:,1)*sig1; Y(:,2)=Y(:,2)*sig2; Y(:,3)=Y(:,3)*sig3;
for t=2:T
    Z=Y(:,:,t)+X;
    diffp=pfun(Z)-pfun(X); 
    logdiff=log(prod(normpdf(X,Z,[sig1 sig2 sig3]),2))- ...
        log(prod(normpdf(Z,X,[sig1 sig2 sig3]),2));
    alphamat=min(diffp+logdiff,0);
    rmat=log(rand(Nsamp,1));
    X=(rmat<alphamat).*Y(:,:,t) + X;
end
toc

tiledlayout(3,1)
ptrue=[A0 B0 f0];
for ind=1:3
    nexttile(ind)
    histogram(X(:,ind),max(120,floor(sqrt(Nsamp))),'normalization','pdf','EdgeColor','none')
    if ind==3
        histogram(fmax/2 + fmax/2*tanh(X(:,ind)),100,'normalization','pdf','EdgeColor','none')
    end
    hold on
    xline(ptrue(ind))
end

%% sampling from multigaussians for each component of posterior
% build param struct using random initial guess of parameters 
NgaussA=2;
% param.muvecA = rand(1,NgaussA);
% param.sdvecA = rand(1,NgaussA);
% param.wvecA  = ones(1,NgaussA)/NgaussA;
NgaussB=2;
% param.muvecB = rand(1,NgaussB);
% param.sdvecB = rand(1,NgaussB);
% param.wvecB  = ones(1,NgaussB)/NgaussB;
NgaussT=4;
% param.muvecT = rand(1,NgaussT);
% param.sdvecT = rand(1,NgaussT);
% param.wvecT  = ones(1,NgaussT)/NgaussT;
% 
% 
% xA=[param.muvecA param.sdvecA param.wvecA];
% xB=[param.muvecB param.sdvecB param.wvecB];
% xT=[param.muvecT param.sdvecT param.wvecT];
% 
% costFunOneCosinorFFVB(xA,xB,xT,NgaussA,NgaussB,NgaussT,settings,Nsamp)

%% try optimizing
opts = optimoptions(@particleswarm,'HybridFcn',@fminsearch,...
                                   'Display','iter', ...    
                                   'SwarmSize',300, ...
                                   'InitialSwarmSpan',15, ...
                                   'UseParallel',true);%,'PlotFcn',@pswmyfun); % TODO make true

dist_param=particleswarm(@(x) -costFunOneCosinorFFVB(x(1:3*NgaussA), ...
                                         x(3*NgaussA+1:3*(NgaussA+NgaussB)), ...
                                         x(3*(NgaussA+NgaussB)+1:end),NgaussA,NgaussB,NgaussT,settings,Nsamp), ...
                                    3*(NgaussA+NgaussB+NgaussT),[ ],[ ],opts);
%%
xA=[dist_param(1:3*NgaussA)];
xB=[dist_param(3*NgaussA+1:3*(NgaussA+NgaussB))];
xT=[dist_param(3*(NgaussA+NgaussB)+1:end)];
param = getParamFromVec(xA,xB,xT,NgaussA,NgaussB,NgaussT);
muvecA=param.muvecA;
sdvecA=param.sdvecA;
wvecA =param.wvecA;
muvecB=param.muvecB;
sdvecB=param.sdvecB;
wvecB =param.wvecB;
muvecT=param.muvecT;
sdvecT=param.sdvecT;
wvecT =param.wvecT;

multinormpdf   = @(x,muvec,sdvec,wvec) ...
       normpdf(repmat(x,1,length(muvec)), ...
               repmat(muvec,length(x),1), ...
               repmat(sdvec,length(x),1))*wvec';

for ind=1:3
    t=nexttile(ind);
    xv=t.XLim(1):.01:t.XLim(2);
    xv=xv';
    switch ind
        case 1
            plot(xv,multinormpdf(xv,muvecA,sdvecA,wvecA))
        case 2
            plot(xv,multinormpdf(xv,muvecB,sdvecB,wvecB))
        case 3
            xv=-100:.01:100;
            xv=xv'
            plot(fmax/2*(1+tanh(xv)),multinormpdf(xv,muvecT,sdvecT,wvecT))
            xlim([0,fmax])
    end
end


%%
function param = getParamFromVec(xA,xB,xT,NgaussA,NgaussB,NgaussT)
param.muvecA=xA(1:NgaussA);
param.sdvecA=abs(xA(NgaussA+1:2*NgaussA));
param.wvecA =abs(xA(2*NgaussA+1:end))/sum(abs(xA(2*NgaussA+1:end)));
param.muvecB=xB(1:NgaussB);
param.sdvecB=abs(xB(NgaussB+1:2*NgaussB));
param.wvecB =abs(xB(2*NgaussB+1:end))/sum(abs(xB(2*NgaussB+1:end)));
param.muvecT=xT(1:NgaussT);
param.sdvecT=abs(xT(NgaussT+1:2*NgaussT));
param.wvecT =abs(xT(2*NgaussT+1:end))/sum(abs(xT(2*NgaussT+1:end)));
end

function LB = costFunOneCosinorFFVB(xA,xB,xT,NgaussA,NgaussB,NgaussT,settings,Nsamp)
% unpack param struct for evaluating cost function
param = getParamFromVec(xA,xB,xT,NgaussA,NgaussB,NgaussT);

muvecA=param.muvecA;
sdvecA=param.sdvecA;
wvecA =param.wvecA;
muvecB=param.muvecB;
sdvecB=param.sdvecB;
wvecB =param.wvecB;
muvecT=param.muvecT;
sdvecT=param.sdvecT;
wvecT =param.wvecT;

% unpack other variables
muAprior  =  settings.muAprior;
sigAprior =  settings.sigAprior;
muBprior  =  settings.muBprior;
sigBprior =  settings.sigBprior;
muTprior  =  settings.muTprior;
sigTprior =  settings.sigTprior;
fmax      =  settings.fmax;
y_obs     =  settings.y_obs;
t_obs     =  settings.t_obs;

hlambda=@(par) -sum((y_obs - par(:,1).*cos(2*pi*0.5*(fmax +fmax*tanh(par(:,3)))*t_obs) ... % likelihood
    - par(:,2).*sin(2*pi*0.5*(fmax +fmax*tanh(par(:,3)))*t_obs)  ).^2/2,2) - ...
    log(sqrt(2*pi)) - ... % normalization of likelihood (sigma=1)
    (muAprior-par(:,1)).^2/2/sigAprior^2 -  ... % prior
    (muBprior-par(:,2)).^2/2/sigBprior^2 - ...  
    (muTprior-par(:,3)).^2/2/sigTprior^2 - ...
    log(sqrt(2*pi)*sigAprior) - ... ; % normalization from prior % TODO Check normalization
    log(sqrt(2*pi)*sigBprior) - ... 
    log(sqrt(2*pi)*sigTprior) - ...
    ( log((  exp(-(par(:,1)-muvecA).^2/2./sdvecA.^2)./(sqrt(2*pi)*sdvecA)  )*wvecA') + ...
    log((  exp(-(par(:,2)-muvecB).^2/2./sdvecB.^2)./(sqrt(2*pi)*sdvecB)  )*wvecB') + ...
    log((  exp(-(par(:,3)-muvecT).^2/2./sdvecT.^2)./(sqrt(2*pi)*sdvecT)  )*wvecT') );
%    (  -(par(:,1)-muvecA).^2/2./sdvecA.^2 -log((sqrt(2*pi)*sdvecA))  +log(wvecA') + ...
%       -(par(:,2)-muvecB).^2/2./sdvecB.^2 -log((sqrt(2*pi)*sdvecB))  +log(wvecB') + ...
%       -(par(:,3)-muvecT).^2/2./sdvecT.^2 -log((sqrt(2*pi)*sdvecT))  +log(wvecT'));

% qlambda=@(par) log((  exp(-(par(:,1)-muvecA).^2/2./sdvecA.^2)./(sqrt(2*pi)*sdvecA)  )*wvecA') + ...
%     log((  exp(-(par(:,2)-muvecB).^2/2./sdvecB.^2)./(sqrt(2*pi)*sdvecB)  )*wvecB') + ...
%     log((  exp(-(par(:,3)-muvecT).^2/2./sdvecT.^2)./(sqrt(2*pi)*sdvecT)  )*wvecT');
par=sampleCosinorFFVB(Nsamp,param);
LB=sum(hlambda(par))/size(par,1);
end

function S=sampleCosinorFFVB(Nsamp,param)
% all three factors of posterior are modeled as multi-Gaussians
muvecA=param.muvecA; % cosine amp
sdvecA=param.sdvecA;
wvecA =param.wvecA;
muvecB=param.muvecB; % sine coeff
sdvecB=param.sdvecB;
wvecB =param.wvecB;
muvecT=param.muvecT; % period 
sdvecT=param.sdvecT;
wvecT =param.wvecT;

mu_indsA =randsample(1:length(muvecA),Nsamp,'true',wvecA);
mu_indsB =randsample(1:length(muvecB),Nsamp,'true',wvecB);
mu_indsT =randsample(1:length(muvecT),Nsamp,'true',wvecT);
S       =randn([Nsamp,3]);
if length(sdvecA)>1
    S(:,1)       =(S(:,1).*sdvecA(mu_indsA)'+muvecA(mu_indsA)'); 
else
    S(:,1)       =(S(:,1).*sdvecA(mu_indsA)+muvecA(mu_indsA)); 
end
if length(sdvecB)>1
    S(:,2)       =(S(:,2).*sdvecB(mu_indsB)'+muvecB(mu_indsB)'); 
else
    S(:,2)       =(S(:,2).*sdvecB(mu_indsB)+muvecB(mu_indsB)); 
end
if length(sdvecT)>1
    S(:,3)       =(S(:,3).*sdvecT(mu_indsT)'+muvecT(mu_indsT)'); 
else 
    S(:,3)       =(S(:,3).*sdvecT(mu_indsT)+muvecT(mu_indsT)); 
end
end

function S=sampleNnormal(Nsamp,muvec,sdvec,wvec)
mu_inds =randsample(1:length(muvec),Nsamp,'true',wvec);
S       =randn([Nsamp,1]);
S       =(S.*sdvec+muvec(mu_inds)');
end

%%

% 
% %%
% % assumes x is a column vector and muvec,sdvec, wvec are all rows
% multinormcdf   = @(x,muvec,sdvec,wvec) ...
%        normcdf(repmat(x,1,length(muvec)), ...
%                repmat(muvec,length(x),1), ...
%                repmat(sdvec,length(x),1))*wvec';
% Nsamp=1e5;
%  
% S=rand(Nsamp,1);
% S0=S;
% [~,ind]=min(abs(S0-cumsum(wvec))');
% S=muvec(ind)';
% for i=1:5
%     S=S-(multinormcdf(S,muvec,sdvec,wvec)-S0)./multinormpdf(S,muvec,sdvec,wvec); 
% end
% close all
% histogram(S,'Normalization','pdf')
%multinormcdf(a,muvec,sdvec,wvec)
