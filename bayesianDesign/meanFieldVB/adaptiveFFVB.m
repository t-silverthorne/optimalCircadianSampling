clear
rng(12)
%% simulated data

sqrw=@(w) w.^2/norm(w.^2);% DEFINED IN OTHER FUNCTION TOO

A0=2.5;
B0=2.5;
f0=2.1;
fmax=15;
T=1000;
Nsamp=1e5;
Nmeas=8;

t_obs=linspace(0,1,Nmeas+1);
t_obs=t_obs(1:end-1);
y_obs=A0*cos(2*pi*f0*t_obs)+B0*sin(2*pi*f0*t_obs) + randn(1,numel(t_obs));

NgaussA=2;
NgaussB=2;
NgaussT=3;
gradDim=3*(NgaussA+NgaussB+NgaussT);

muAprior  =0;
sigAprior =2;
muBprior  =0;
sigBprior =2;
muTprior  =0;
sigTprior =2;

% muvecA=rand(1,NgaussA);
% sdvecA=rand(1,NgaussA);
% wvecA=rand(1,NgaussA);
% 
% muvecB=rand(1,NgaussB);
% sdvecB=rand(1,NgaussB);
% wvecB=rand(1,NgaussB);
% 
% muvecT=rand(1,NgaussT);
% sdvecT=rand(1,NgaussT);
% wvecT=rand(1,NgaussT);

% want the gradient of this
% logq=( log((  exp(-(S1-muvecA).^2/2./sdvecA.^2)./(sqrt(2*pi)*sdvecA)  )*wvecA' ) + ...
% log((  exp(-(S2-muvecB).^2/2./sdvecB.^2)./(sqrt(2*pi)*sdvecB)  )*wvecB') + ...
% log((  exp(-(S3-muvecT).^2/2./sdvecT.^2)./(sqrt(2*pi)*sdvecT)  )*wvecT') );

% this is the gradient
muvecAgrad = @(S,muvecA,sdvecA,wvecA) ...
           ( (2^(1/2)*exp(-(S(:,1) - muvecA).^2./(2*sdvecA.^2)).*(2*S(:,1) - ...
           2*muvecA))./(4*pi^(1/2)*abs(sdvecA).^3).*sqrw(wvecA)  )./...
           ( exp(-(S(:,1)-muvecA).^2/2./sdvecA.^2)./(sqrt(2*pi)*abs(sdvecA))*sqrw(wvecA)'  ) ;
sdvecAgrad = @(S,muvecA,sdvecA,wvecA) -...
           (2.^(1./2).*exp(-(S(:,1) - muvecA).^2./(2.*sdvecA.^2)).*sign(sdvecA).*(2.*S(:,1).*muvecA + ...
           sdvecA.^2.*sign(sdvecA).^2 - S(:,1).^2 - muvecA.^2))./(2.*pi.^(1./2).*sdvecA.^4)./...
           ( exp(-(S(:,1)-muvecA).^2/2./sdvecA.^2)./(sqrt(2*pi)*abs(sdvecA))*sqrw(wvecA)'  );
wvecAgrad  = @(S,muvecA,sdvecA,wvecA) ...
           (2*wvecA/norm(wvecA.^2) - 2*wvecA.^5/norm(wvecA.^2)^3 ).*exp(-(S(:,1) - muvecA).^2./(2*sdvecA.^2))./...
           ( exp(-(S(:,1)-muvecA).^2/2./sdvecA.^2)./(sqrt(2*pi)*abs(sdvecA))*sqrw(wvecA)'  );

muvecBgrad = @(S,muvecB,sdvecB,wvecB) ... 
           ( (2^(1/2)*exp(-(S(:,2) - muvecB).^2./(2*sdvecB.^2)).*(2*S(:,2) - ...
           2*muvecB))./(4*pi^(1/2)*abs(sdvecB).^3).*sqrw(wvecB)  )./...
           ( exp(-(S(:,2)-muvecB).^2/2./sdvecB.^2)./(sqrt(2*pi)*abs(sdvecB))*sqrw(wvecB)'  ) ;
sdvecBgrad = @(S,muvecB,sdvecB,wvecB) -...
           (2.^(1./2).*exp(-(S(:,2) - muvecB).^2./(2.*sdvecB.^2)).*sign(sdvecB).*(2.*S(:,2).*muvecB + ...
           sdvecB.^2.*sign(sdvecB).^2 - S(:,2).^2 - muvecB.^2))./(2.*pi.^(1./2).*sdvecB.^4)./...
           ( exp(-(S(:,2)-muvecB).^2/2./sdvecB.^2)./(sqrt(2*pi)*abs(sdvecB))*sqrw(wvecB)'  );
wvecBgrad  = @(S,muvecB,sdvecB,wvecB) ...
           (2*wvecB/norm(wvecB.^2) - 2*wvecB.^5/norm(wvecB.^2)^3 ).*exp(-(S(:,2) - muvecB).^2./(2*sdvecB.^2))./...
           ( exp(-(S(:,2)-muvecB).^2/2./sdvecB.^2)./(sqrt(2*pi)*abs(sdvecB))*sqrw(wvecB)'  );

muvecTgrad = @(S,muvecT,sdvecT,wvecT) ... 
           ( (2^(1/2)*exp(-(S(:,3) - muvecT).^2./(2*sdvecT.^2)).*(2*S(:,3) - ...
           2*muvecT))./(4*pi^(1/2)*abs(sdvecT).^3).*sqrw(wvecT)  )./...
           ( exp(-(S(:,3)-muvecT).^2/2./sdvecT.^2)./(sqrt(2*pi)*abs(sdvecT))*sqrw(wvecT)'  ) ;
sdvecTgrad = @(S,muvecT,sdvecT,wvecT) -...
           (2.^(1./2).*exp(-(S(:,3) - muvecT).^2./(2.*sdvecT.^2)).*sign(sdvecT).*(2.*S(:,3).*muvecT + ...
           sdvecT.^2.*sign(sdvecT).^2 - S(:,3).^2 - muvecT.^2))./(2.*pi.^(1./2).*sdvecT.^4)./...
           ( exp(-(S(:,3)-muvecT).^2/2./sdvecT.^2)./(sqrt(2*pi)*abs(sdvecT))*sqrw(wvecT)'  );
wvecTgrad  = @(S,muvecT,sdvecT,wvecT) ...
           (2*wvecT/norm(wvecT.^2) - 2*wvecT.^5/norm(wvecT.^2)^3 ).*exp(-(S(:,3) - muvecT).^2./(2*sdvecT.^2))./...
           ( exp(-(S(:,3)-muvecT).^2/2./sdvecT.^2)./(sqrt(2*pi)*abs(sdvecT))*sqrw(wvecT)'  );
logqgrad=@(S,muvecA,sdvecA,wvecA, ...
             muvecB,sdvecB,wvecB, ...
             muvecT,sdvecT,wvecT) ...
             [muvecAgrad(S,muvecA,sdvecA,wvecA) sdvecAgrad(S,muvecA,sdvecA,wvecA) wvecAgrad(S,muvecA,sdvecA,wvecA) ...
              muvecBgrad(S,muvecB,sdvecB,wvecB) sdvecBgrad(S,muvecB,sdvecB,wvecB) wvecBgrad(S,muvecB,sdvecB,wvecB) ...
              muvecTgrad(S,muvecT,sdvecT,wvecT) sdvecTgrad(S,muvecT,sdvecT,wvecT) wvecTgrad(S,muvecT,sdvecT,wvecT)];


% evaluate
% S=rand(10,3);
% logqgrad(S,muvecA,sdvecA,wvecA, ...
%              muvecB,sdvecB,wvecB, ...
%              muvecT,sdvecT,wvecT)        
%%
% syms sig mu S pi
% %simplify(diff(1/sqrt(2*pi)/sig*exp(-(mu-S)^2/2/sig^2),mu))
% %(2^(1/2)*exp(-(S - mu)^2/(2*sig^2))*(2*S - 2*mu))/(4*pi^(1/2)*sig^3)
% 
% vectorize(simplify(diff(1/sqrt(2*pi)/abs(sig)*exp(-(mu-S)^2/2/sig^2),sig)))
%-(2^(1/2)*exp(-(S - mu)^2/(2*sig^2))*sign(sig)*(2*S*mu + sig^2*sign(sig)^2 - S^2 - mu^2))/(2*pi^(1/2)*sig^4)
%%
% random starting lambda
%lambda=rand(1,3*(NgaussA+NgaussB+NgaussT));

% func start

% unpack global variables
%xA        =  lambda(1:3*NgaussA);
%xB        =  lambda(3*NgaussA+1:3*(NgaussA+NgaussB));
%xT        =  lambda(3*(NgaussA+NgaussB)+1:3*(NgaussA+NgaussB+NgaussT));
sc=.1;
muvecA    =  zeros(1,NgaussA)+sc*rand(1,NgaussA);
sdvecA    =  ones(1,NgaussA)+sc*rand(1,NgaussA);
wvecA     =  sqrw(ones(1,NgaussA));
muvecB    =  zeros(1,NgaussB)+sc*rand(1,NgaussB);
sdvecB    =  ones(1,NgaussB)+sc*rand(1,NgaussB);
wvecB     =  sqrw(ones(1,NgaussB));
muvecT    =  atanh(linspace(1/2,fmax,NgaussT)/fmax-1);
sdvecT    =  2*ones(1,NgaussT)+sc*rand(1,NgaussT);
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
beta1=.5; beta2=0.5; eps0=5;tau=500;tW=30;maxPat=80;
t=0;patience=0;stop=false;
LBvec=[];
LBwindowvec=[];

while stop==false && t<2000
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
    hlambdavalc=hlambda(S,muvecA,sdvecA,wvecA, ...
             muvecB,sdvecB,wvecB, ...
             muvecT,sdvecT,wvecT);
    logqgradvalc = logqgrad(S,muvecA,sdvecA,wvecA, ...
             muvecB,sdvecB,wvecB, ...
             muvecT,sdvecT,wvecT);
    gr1=logqgradvalc.*hlambdavalc;
    gr2=logqgradvalc;
    c=arrayfun(@(ind) [1 0]*cov(gr1(:,ind),gr2(:,ind))*[0;1]/var(gr2(:,ind)),1:gradDim); %control variate
    gt=sum(logqgradval.*(hlambdaval-c),1)/Nsamp;
    gt(isnan(gt))=0;
    %disp(gt)
    c(isnan(c))=0; % ad hoc fix for guys with zero variance
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
    wvecA     =  xA(2*NgaussA+1:end)/sum(abs(xA(2*NgaussA+1:end)));
    muvecB    =  xB(1:NgaussB);
    sdvecB    =  xB(NgaussB+1:2*NgaussB);
    wvecB     =  xB(2*NgaussB+1:end)/sum(abs(xB(2*NgaussB+1:end)));
    muvecT    =  xT(1:NgaussT);
    sdvecT    =  xT(NgaussT+1:2*NgaussT);
    wvecT     =  xT(2*NgaussT+1:end)/sum(abs(xT(2*NgaussT+1:end)));


    % compute lower bound estimate
    hlambdaval=hlambda(S,muvecA,sdvecA,wvecA, ...
                 muvecB,sdvecB,wvecB, ...
                 muvecT,sdvecT,wvecT);
    LBvec(end+1)=sum(hlambdaval(hlambdaval<Inf))/Nsamp;%adhoc
    if t>=tW
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
close all;plot(LBwindowvec)
% %%
% close all
% plot(1:length(LBvec),LBvec)
% %%
% %% MCMC
% tic
% sig1=2;sig2=2;sig3=2;
% X=randn(Nsamp,3);
% X(:,1)=X(:,1)*sigAprior+muAprior; X(:,2)=X(:,2)*sigBprior+muBprior;
% X(:,1)=X(:,1)*sigTprior+muTprior;
% 
% pfun=@(par) -sum((y_obs - par(:,1).*cos(2*pi*0.5*(fmax+fmax*tanh(par(:,3)))*t_obs) ... % likelihood
%     - par(:,2).*sin(2*pi*0.5*(fmax +fmax*tanh(par(:,3)))*t_obs)  ).^2/2,2) - ...
%     log(sqrt(2*pi)) - ... % normalization of likelihood (sigma=1)
%     (muAprior-par(:,1)).^2/2/sigAprior^2 -  ... % prior
%     (muBprior-par(:,2)).^2/2/sigBprior^2 - ...  
%     (muTprior-par(:,3)).^2/2/sigTprior^2;
% 
% Y=randn(Nsamp,3,T);
% Y(:,1)=Y(:,1)*sig1; Y(:,2)=Y(:,2)*sig2; Y(:,3)=Y(:,3)*sig3;
% for t=2:T
%     Z=Y(:,:,t)+X;
%     diffp=pfun(Z)-pfun(X); 
%     logdiff=log(prod(normpdf(X,Z,[sig1 sig2 sig3]),2))- ...
%         log(prod(normpdf(Z,X,[sig1 sig2 sig3]),2));
%     alphamat=min(diffp+logdiff,0);
%     rmat=log(rand(Nsamp,1));
%     X=(rmat<alphamat).*Y(:,:,t) + X;
% end
% toc
% 
% tiledlayout(3,1)
% ptrue=[A0 B0 f0];
% for ind=1:3
%     nexttile(ind)
%     histogram(X(:,ind),max(120,floor(sqrt(Nsamp))),'normalization','pdf','EdgeColor','none')
%     if ind==3
%         histogram(fmax/2 + fmax/2*tanh(X(:,ind)),100,'normalization','pdf','EdgeColor','none')
%     end
%     hold on
%     xline(ptrue(ind))
% end
% 
% %% comparison
% 
% multinormpdf   = @(x,muvec,sdvec,wvec) ...
%        normpdf(repmat(x,1,length(muvec)), ...
%                repmat(muvec,length(x),1), ...
%                repmat(abs(sdvec),length(x),1))*sqrw(wvec)';
% figure(1)
% for ind=1:3
%     t=nexttile(ind);
%     xv=t.XLim(1):.01:t.XLim(2);
%     xv=xv';
%     switch ind
%         case 1
%             plot(xv,multinormpdf(xv,muvecA,sdvecA,wvecA))
%         case 2
%             plot(xv,multinormpdf(xv,muvecB,sdvecB,wvecB))
%         case 3
%             xv=0:.01:fmax;
%             xv=xv';
%             plot(xv,2/fmax*multinormpdf(atanh(2*xv/fmax-1),muvecT,sdvecT,wvecT)./abs((( 2*xv/fmax).^2 -1)))
%             xlim([0,fmax])
%     end
% end


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

