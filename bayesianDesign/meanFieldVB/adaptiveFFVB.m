figure(1)
clf
figure(2)
clf
rng(12)
%% simulated data

sqrw =@(w) w.^2/sum(w.^2);
sqrw2=@(w) w.^2/sum(w.^2)^2;
A0=2.5;
B0=2.5;
f0=2.1;
fmax=16;
T=1000;
Nsamp=1e4;
Nmeas=8;

t_obs=linspace(0,1,Nmeas+1);
t_obs=t_obs(1:end-1);
y_obs=A0*cos(2*pi*f0*t_obs)+B0*sin(2*pi*f0*t_obs) + randn(1,numel(t_obs));

NgaussA=1;
NgaussB=2;
NgaussT=3;
gradDim=3*(NgaussA+NgaussB+NgaussT);

muAprior  =0;
sigAprior =2;
muBprior  =0;
sigBprior =2;
muTprior  =0;
sigTprior =2;


% derivatives wrt mu param
muvecgrad = @(S,muvec,sdvec,wvec)( (2.^(1./2).*exp(-(S - muvec).^2./(2.*sdvec.^2)).*(2.*S - 2.*muvec))./(4.*pi.^(1./2).*sdvec.^2.*abs(sdvec)) ).*sqrw(wvec)./ ...
           (  (2.^(1./2).*exp(-(S - muvec).^2./(2.*sdvec.^2)))./(2.*pi.^(1./2).*abs(sdvec))*sqrw(wvec)'  );
% derivatives wrt sdev param
sdvecgrad = @(S,muvec,sdvec,wvec)( -(2.^(1./2).*exp(-(S - muvec).^2./(2.*sdvec.^2)).*sign(sdvec).*(2.*S.*muvec + ... 
                 sdvec.^2.*sign(sdvec).^2 - S.^2 - muvec.^2))./(2.*pi.^(1./2).*sdvec.^4) ).*sqrw(wvec)./ ...
           (  (2.^(1./2).*exp(-(S - muvec).^2./(2.*sdvec.^2)))./(2.*pi.^(1./2).*abs(sdvec))*sqrw(wvec)'  );
% derivatives wrt weights
wvecgrad  = @(S,muvec,sdvec,wvec)( -2*wvec.*( (2.^(1./2).*exp(-(S - muvec).^2./(2.*sdvec.^2)))./(2.*pi.^(1./2).*abs(sdvec))*sqrw2(wvec)'  ) + ... 
                                      2*wvec.*exp(-(muvec-S).^2/2./sdvec.^2)/sqrt(2*pi)./abs(sdvec)/sum(wvec.^2) )./ ...
           (  (2.^(1./2).*exp(-(S - muvec).^2./(2.*sdvec.^2)))./(2.*pi.^(1./2).*abs(sdvec))*sqrw(wvec)'  );
logqgrad=@(S,muvecA,sdvecA,wvecA, ...
             muvecB,sdvecB,wvecB, ...
             muvecT,sdvecT,wvecT) ...
             [muvecgrad(S(:,1),muvecA,sdvecA,wvecA) sdvecgrad(S(:,1),muvecA,sdvecA,wvecA) wvecgrad(S(:,1),muvecA,sdvecA,wvecA) ...
              muvecgrad(S(:,2),muvecB,sdvecB,wvecB) sdvecgrad(S(:,2),muvecB,sdvecB,wvecB) wvecgrad(S(:,2),muvecB,sdvecB,wvecB) ...
              muvecgrad(S(:,3),muvecT,sdvecT,wvecT) sdvecgrad(S(:,3),muvecT,sdvecT,wvecT) wvecgrad(S(:,3),muvecT,sdvecT,wvecT)];


muvecA    =  linspace(-5,5,NgaussA);
sdvecA    =  1*ones(1,NgaussA);
wvecA     =  sqrw(ones(1,NgaussA));
muvecB    =  linspace(-5,5,NgaussB);
sdvecB    =  1*ones(1,NgaussB);
wvecB     =  sqrw(ones(1,NgaussB));
muvecT    =  atanh(2*[2 2 6]/fmax-1);%atanh(linspace(1/2,fmax,NgaussT)/fmax-1);
sdvecT    =  .1*ones(1,NgaussT);
wvecT     =  sqrw(ones(1,NgaussT));

lambda=[muvecA sdvecA wvecA ...
        muvecB sdvecB wvecB ...
        muvecT sdvecT wvecT];
%%

% INITIALIZATION
S=sample_from_q(Nsamp,muvecA,sdvecA,sqrw(wvecA),muvecB, ...
                        sdvecB,sqrw(wvecB),muvecT,sdvecT,sqrw(wvecT));

%%
% UNBIASED ESTIMATE OF LB GRADIENT
hlambda=@(S,muvecA,sdvecA,wvecA, ...
             muvecB,sdvecB,wvecB, ...
             muvecT,sdvecT,wvecT) ...
             -sum((y_obs - S(:,1).*cos(2*pi*0.5*(fmax +fmax*tanh(S(:,3)))*t_obs) ... % likelihood
    - S(:,2).*sin(2*pi*0.5*(fmax +fmax*tanh(S(:,3)))*t_obs)  ).^2/2,2) - ...
    log(sqrt(2*pi)) - ... % normalization of likelihood (sigma=1)
    (muAprior-S(:,1)).^2/2/sigAprior^2 -  ... % prior
    (muBprior-S(:,2)).^2/2/sigBprior^2 - ...  
    (muTprior-S(:,3)).^2/2/sigTprior^2 - ...
    log(sqrt(2*pi)*sigAprior) - ... ; % normalization from prior % TODO Check normalization
    log(sqrt(2*pi)*sigBprior) - ... 
    log(sqrt(2*pi)*sigTprior) - ...
    (log((  exp(-(S(:,1)-muvecA).^2/2./sdvecA.^2)./(sqrt(2*pi)*abs(sdvecA))  )*sqrw(wvecA)') + ...
     log((  exp(-(S(:,2)-muvecB).^2/2./sdvecB.^2)./(sqrt(2*pi)*abs(sdvecB))  )*sqrw(wvecB)') + ...
     log((  exp(-(S(:,3)-muvecT).^2/2./sdvecT.^2)./(sqrt(2*pi)*abs(sdvecT))  )*sqrw(wvecT)') );
hlambdaval=hlambda(S,muvecA,sdvecA,wvecA, ...
             muvecB,sdvecB,wvecB, ...
             muvecT,sdvecT,wvecT);
logqgradval = logqgrad(S,muvecA,sdvecA,wvecA, ...
             muvecB,sdvecB,wvecB, ...
             muvecT,sdvecT,wvecT);
gradLBhat=sum(logqgradval.*hlambdaval,1)/Nsamp;

g0=gradLBhat;
v0=g0.^2;
gbar=g0;
vbar=v0;

% learning params
beta1=.1; beta2=0.1; eps0=1e-3;tau=2000;tW=30;maxPat=1000;
t=0;patience=0;stop=false;
LBvec=[];
LBwindowvec=[];
gtvec=[];
close all
tiledlayout(3,1)
while stop==false && t<500

    if visualize
        % visualization
        multinormpdf   = @(x,muvec,sdvec,wvec) ...
               normpdf(repmat(x,1,length(muvec)), ...
                       repmat(muvec,length(x),1), ...
                       repmat(abs(sdvec),length(x),1))*sqrw(wvec)';
        for ind=1:3
            nexttile(ind)
            xv=-10:.1:10;
            xv=xv';
            switch ind
                case 1
                    plot(xv,multinormpdf(xv,muvecA,sdvecA,wvecA))
                case 2
                    plot(xv,multinormpdf(xv,muvecB,sdvecB,wvecB))
                case 3
                    xv=0:.01:fmax;
                    xv=xv';
                    plot(xv,2/fmax*multinormpdf(atanh(2*xv/fmax-1),muvecT,sdvecT,wvecT)./abs((( 2*xv/fmax).^2 -1)))
                    xlim([0,fmax])
            end
        end
        drawnow
    end

    fprintf('%d\n',t)
    S=sample_from_q(Nsamp,muvecA,sdvecA,sqrw(wvecA),muvecB, ...
                        sdvecB,sqrw(wvecB),muvecT,sdvecT,sqrw(wvecT));
    hlambdaval=hlambda(S,muvecA,sdvecA,wvecA, ...
             muvecB,sdvecB,wvecB, ...
             muvecT,sdvecT,wvecT);
    logqgradval = logqgrad(S,muvecA,sdvecA,wvecA, ...
             muvecB,sdvecB,wvecB, ...
             muvecT,sdvecT,wvecT);
    
    Sc=sample_from_q(Nsamp,muvecA,sdvecA,sqrw(wvecA),muvecB, ...
                        sdvecB,sqrw(wvecB),muvecT,sdvecT,sqrw(wvecT));
    hlambdavalc=hlambda(Sc,muvecA,sdvecA,wvecA, ...
             muvecB,sdvecB,wvecB, ...
             muvecT,sdvecT,wvecT);
    logqgradvalc = logqgrad(Sc,muvecA,sdvecA,wvecA, ...
             muvecB,sdvecB,wvecB, ...
             muvecT,sdvecT,wvecT);
    gr1=logqgradvalc.*hlambdavalc;
    gr2=logqgradvalc;
    c=arrayfun(@(ind) [1 0]*cov(gr1(:,ind),gr2(:,ind))*[0;1]/var(gr2(:,ind)),1:gradDim); %control variate
    c(abs(c)>1e5)=0; % ad hoc way of dealing with low covariance terms
    gt=sum(logqgradval.*(hlambdaval-c),1)/Nsamp;
    if norm(gt)>1e8
        disp('huge gradient, possibly unstable')
    end
    gt(isnan(gt))=0; % ad hoc
    %disp(gt)
    %c(isnan(c))=0; % ad hoc fix for guys with zero variance
    vt=gt.^2;
    gbar=beta1*gbar+(1-beta1)*gt;
    vbar=beta2*vbar+(1-beta2)*vt;
    
    alphat=min(eps0,eps0*tau/t);
    if sum(isnan(gt) | isnan(vt))
        disp('error imaginary learning velocity')
    end
    lambda=lambda+alphat*gbar./sqrt(vbar);
    %lambda(isnan(lambda))=1; % ad hoc
    
    xA        =  lambda(1:3*NgaussA);
    xB        =  lambda(3*NgaussA+1:3*(NgaussA+NgaussB));
    xT        =  lambda(3*(NgaussA+NgaussB)+1:3*(NgaussA+NgaussB+NgaussT));
    muvecA    =  xA(1:NgaussA);
    sdvecA    =  xA(NgaussA+1:2*NgaussA);
    wvecA     =  xA(2*NgaussA+1:end);
    muvecB    =  xB(1:NgaussB);
    sdvecB    =  xB(NgaussB+1:2*NgaussB);
    wvecB     =  xB(2*NgaussB+1:end);
    muvecT    =  xT(1:NgaussT);
    sdvecT    =  xT(NgaussT+1:2*NgaussT);
    wvecT     =  xT(2*NgaussT+1:end);


    % compute lower bound estimate
    hlambdaval=hlambda(S,muvecA,sdvecA,wvecA, ...
                 muvecB,sdvecB,wvecB, ...
                 muvecT,sdvecT,wvecT);
    %LBvec(end+1)=sum(hlambdaval(hlambdaval<Inf))/Nsamp;%adhoc

    LBvec(end+1)=sum(hlambdaval)/Nsamp;
    gtvec(end+1)=norm(gt);
    if t>=1e3 % tW ad hoc
        LBwindow= sum(LBvec(end-tW+1:end))/tW;
        LBwindowvec(end+1)=LBwindow;
        if LBwindow>max(LBwindowvec(1:end-1))
            patience=0;
        else
            patience=patience+1;
        end
    end
    if patience>maxPat
        stop=true;
    end
    t=t+1;


end

%%
disp('continued')
close all
figure(2)
tiledlayout(3,1)
nexttile
plot(1:length(LBvec),LBvec)
nexttile
plot(LBwindowvec)
nexttile
semilogy(gtvec)
figure(1)
%% MCMC
tic
sig1=2;sig2=2;sig3=2;
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

%% comparison

multinormpdf   = @(x,muvec,sdvec,wvec) ...
       normpdf(repmat(x,1,length(muvec)), ...
               repmat(muvec,length(x),1), ...
               repmat(abs(sdvec),length(x),1))*sqrw(wvec)';
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
            xv=0:.01:fmax;
            xv=xv';
            plot(xv,2/fmax*multinormpdf(atanh(2*xv/fmax-1),muvecT,sdvecT,wvecT)./abs((( 2*xv/fmax).^2 -1)))
            xlim([0,fmax])
    end
end


function S=sample_from_q(Nsamp,muvecA,sdvecA,wvecA,muvecB, ...
                        sdvecB,wvecB,muvecT,sdvecT,wvecT)
    % generate random samples from qlambda(theta)
    mu_indsA =randsample(1:length(muvecA),Nsamp,'true',wvecA);
    mu_indsB =randsample(1:length(muvecB),Nsamp,'true',wvecB);
    mu_indsT =randsample(1:length(muvecT),Nsamp,'true',wvecT);
    
    % INITIALIZATION: sample from qlambda
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

