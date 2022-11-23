clear
A0=2.5;
B0=2.5;
f0=2.1;
fmax=12;
Lcut=fmax;
T=1000;
Nsamp=3;
Nmeas=8;

t_obs=linspace(0,1,Nmeas+1);
t_obs=t_obs(1:end-1);
y_obs=A0*cos(2*pi*f0*t_obs)+B0*sin(2*pi*f0*t_obs) + randn(1,numel(t_obs));
%%
param.NgaussA=1;
param.NgaussB=2;
param.NgaussT=4;

%% pre compute symbolic gradient
syms S1 S2 S3
muvecA=sym('muA',[1 param.NgaussA]);
muvecB=sym('muB',[1 param.NgaussB]);
muvecT=sym('muT',[1 param.NgaussT]);
sdvecA=sym('sdA',[1 param.NgaussA]);
sdvecB=sym('sdB',[1 param.NgaussB]);
sdvecT=sym('sdT',[1 param.NgaussT]);
wvecA=sym('wA',[1 param.NgaussA]);
wvecB=sym('wB',[1 param.NgaussB]);
wvecT=sym('wT',[1 param.NgaussT]);
assume([muvecA muvecB muvecT ...
        sdvecA sdvecB sdvecT ...
        wvecA  wvecB  wvecT ...
        S1     S2     S3],'real')

logq=( log((  exp(-(S1-muvecA).^2/2./sdvecA.^2)./(sqrt(2*pi)*sdvecA)  )*wvecA') + ...
log((  exp(-(S2-muvecB).^2/2./sdvecB.^2)./(sqrt(2*pi)*sdvecB)  )*wvecB') + ...
log((  exp(-(S3-muvecT).^2/2./sdvecT.^2)./(sqrt(2*pi)*sdvecT)  )*wvecT') );

% gsym=gradient(logq,[muvecA, sdvecA, wvecA, ...
%             muvecB, sdvecB, wvecB, ...
%             muvecT, sdvecT, wvecT])'
% 
% muvecAval=rand(1,numel(muvecA)); sdvecAval=rand(1,numel(muvecA)); wvecAval=rand(1,numel(wvecA));
% muvecBval=rand(1,numel(muvecB)); sdvecBval=rand(1,numel(muvecB)); wvecBval=rand(1,numel(wvecB));
% muvecTval=rand(1,numel(muvecT)); sdvecTval=rand(1,numel(muvecT)); wvecTval=rand(1,numel(wvecT));
% 
% subs(gsym,[muvecA, sdvecA, wvecA, ...
%             muvecB, sdvecB, wvecB, ...
%             muvecT, sdvecT, wvecT], ...
%             [muvecAval, sdvecAval, wvecAval, ...
%             muvecBval, sdvecBval, wvecBval, ...
%             muvecTval, sdvecTval, wvecTval])
% gradlogq=matlabFunction(subs(gsym,[muvecA, sdvecA, wvecA, ...
%             muvecB, sdvecB, wvecB, ...
%             muvecT, sdvecT, wvecT], ...
%             [muvecAval, sdvecAval, wvecAval, ...
%             muvecBval, sdvecBval, wvecBval, ...
%             muvecTval, sdvecTval, wvecTval]));

%%
gradlogq=matlabFunction(gradient(logq,[muvecA, sdvecA, wvecA, ...
            muvecB, sdvecB, wvecB, ...
            muvecT, sdvecT, wvecT])','Vars',[...
            muvecA muvecB muvecT ...
            sdvecA sdvecB sdvecT ...
            wvecA  wvecB  wvecT ...
            S1     S2     S3]);
rr=rand(1,24);
rrc=num2cell(rr)
size(gradlogq(rrc{1:end-1},[1]))
%%



param.fmax=10;
param.y_obs=y_obs;
param.t_obs=t_obs;
param.muAprior  =0;
param.sigAprior =2;
param.muBprior  =0;
param.sigBprior =2;
param.muTprior  =0;
param.sigTprior =1;

lambda=rand(1,3*(param.NgaussA+param.NgaussB+param.NgaussT));

% func start
% unpack
fmax      =  param.fmax;     % priors, settings, observations, Ngausses
y_obs     =  param.y_obs;
t_obs     =  param.t_obs;
NgaussA   =  param.NgaussA;
NgaussB   =  param.NgaussB;
NgaussT   =  param.NgaussT;
muAprior  =  param.muAprior;
sigAprior =  param.sigAprior;
muBprior  =  param.muBprior;
sigBprior =  param.sigBprior;
muTprior  =  param.muTprior;
sigTprior =  param.sigTprior;
% keep unpacking
xA=lambda(1:3*NgaussA);
xB=lambda(3*NgaussA+1:3*(NgaussA+NgaussB));
xT=lambda(3*(NgaussA+NgaussB)+1:3*(NgaussA+NgaussB+NgaussT));
muvecA=xA(1:NgaussA);
sdvecA=abs(xA(NgaussA+1:2*NgaussA));
wvecA =abs(xA(2*NgaussA+1:end))/sum(abs(xA(2*NgaussA+1:end)));
muvecB=xB(1:NgaussB);
sdvecB=abs(xB(NgaussB+1:2*NgaussB));
wvecB =abs(xB(2*NgaussB+1:end))/sum(abs(xB(2*NgaussB+1:end)));
muvecT=xT(1:NgaussT);
sdvecT=abs(xT(NgaussT+1:2*NgaussT));
wvecT =abs(xT(2*NgaussT+1:end))/sum(abs(xT(2*NgaussT+1:end)));

% generate random samples from qlambda(theta)
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

% compute LB from S
LB=-sum((y_obs - S(:,1).*cos(2*pi*0.5*(fmax +fmax*tanh(S(:,3)))*t_obs) ... % likelihood
    - S(:,2).*sin(2*pi*0.5*(fmax +fmax*tanh(S(:,3)))*t_obs)  ).^2/2,2) - ...
    log(sqrt(2*pi)) - ... % normalization of likelihood (sigma=1)
    (muAprior-S(:,1)).^2/2/sigAprior^2 -  ... % prior
    (muBprior-S(:,2)).^2/2/sigBprior^2 - ...  
    (muTprior-S(:,3)).^2/2/sigTprior^2 - ...
    log(sqrt(2*pi)*sigAprior) - ... ; % normalization from prior % TODO Check normalization
    log(sqrt(2*pi)*sigBprior) - ... 
    log(sqrt(2*pi)*sigTprior) - ...
    ( log((  exp(-(S(:,1)-muvecA).^2/2./sdvecA.^2)./(sqrt(2*pi)*sdvecA)  )*wvecA') + ...
    log((  exp(-(S(:,2)-muvecB).^2/2./sdvecB.^2)./(sqrt(2*pi)*sdvecB)  )*wvecB') + ...
    log((  exp(-(S(:,3)-muvecT).^2/2./sdvecT.^2)./(sqrt(2*pi)*sdvecT)  )*wvecT') );

% symbolic gradient of log(qlambda(x))

%%
muvecAc=num2cell(muvecA);  muvecBc=num2cell(muvecB); muvecTc=num2cell(muvecT);
sdvecAc=num2cell(sdvecA);  sdvecBc=num2cell(sdvecB); sdvecTc=num2cell(sdvecT);
wvecAc=num2cell(wvecA);    wvecBc=num2cell(wvecB);    wvecTc=num2cell(wvecT);

arrayfun(@(ind) gradlogq(   muvecAc{:},  muvecBc{:}, muvecTc{:}, ...
            sdvecAc{:},  sdvecBc{:}, sdvecTc{:}, ...
            wvecAc{:},   wvecBc{:},  wvecTc{:}, ...
            S(ind,1),  S(ind,2), S(ind,3)), (1:size(S,1))','UniformOutput',false)

%LB=@(lambda) getLB(lambda,param)
%gradLB=@(lambda)

