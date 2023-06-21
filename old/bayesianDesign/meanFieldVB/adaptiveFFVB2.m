%syms sdvec muvec S pi wvec
vectorize(simplify(diff(1/sqrt(2*pi)/abs(sdvec)*exp(-(muvec-S)^2/2/sdvec^2),muvec)))
vectorize(simplify(diff(1/sqrt(2*pi)/abs(sdvec)*exp(-(muvec-S)^2/2/sdvec^2),sdvec)))

%%
clear
sqrw =@(w) w.^2/sum(w.^2);
sqrw2=@(w) w.^2/sum(w.^2)^2;

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
% derivatives wrt mu param
muvecgrad = @(S,muvec,sdvec,wvec)( (2.^(1./2).*exp(-(S - muvec).^2./(2.*sdvec.^2)).*(2.*S - 2.*muvec))./(4.*pi.^(1./2).*sdvec.^2.*abs(sdvec)) ).*sqrw(wvec)./ ...
           (  (2.^(1./2).*exp(-(S - muvec).^2./(2.*sdvec.^2)))./(2.*pi.^(1./2).*abs(sdvec))*sqrw(wvec)'  );
% derivatives wrt sdev param
sdvecgrad = @(S,muvec,sdvec,wvec)( -(2.^(1./2).*exp(-(S - muvec).^2./(2.*sdvec.^2)).*sign(sdvec).*(2.*S.*muvec + ... 
                 sdvec.^2.*sign(sdvec).^2 - S.^2 - muvec.^2))./(2.*pi.^(1./2).*sdvec.^4) ).*sqrw(wvec)./ ...
           (  (2.^(1./2).*exp(-(S - muvec).^2./(2.*sdvec.^2)))./(2.*pi.^(1./2).*abs(sdvec))*sqrw(wvec)'  );
% derivatives wrt weights
wvecgrad  = @(S,muvec,sdvec,wvec)( -2*wvec.*( (2.^(1./2).*exp(-(S - muvec).^2./(2.*sdvec.^2)))./(2.*pi.^(1./2).*abs(sdvec))*sqrw2(wvec)'  ) + ... 
                                      2*wvec.*exp(-(muvec-S).^2/2/sdvec.^2)/sqrt(2*pi)./abs(sdvec) )./ ...
           (  (2.^(1./2).*exp(-(S - muvec).^2./(2.*sdvec.^2)))./(2.*pi.^(1./2).*abs(sdvec))*sqrw(wvec)'  );
logqgrad=@(S,muvecA,sdvecA,wvecA, ...
             muvecB,sdvecB,wvecB, ...
             muvecT,sdvecT,wvecT) ...
             [muvecgrad(S,muvecA,sdvecA,wvecA) sdvecgrad(S,muvecA,sdvecA,wvecA) wvecgrad(S,muvecA,sdvecA,wvecA) ...
              muvecgrad(S,muvecB,sdvecB,wvecB) sdvecgrad(S,muvecB,sdvecB,wvecB) wvecgrad(S,muvecB,sdvecB,wvecB) ...
              muvecgrad(S,muvecT,sdvecT,wvecT) sdvecgrad(S,muvecT,sdvecT,wvecT) wvecgrad(S,muvecT,sdvecT,wvecT)];

