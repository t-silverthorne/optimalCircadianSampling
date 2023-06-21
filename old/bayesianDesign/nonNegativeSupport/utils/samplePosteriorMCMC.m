function X = samplePosteriorMCMC(N,yvec,tvec,settings)
mu1=settings.mu1;
mu2=settings.mu2;
mu3=settings.mu3;
sig1=settings.sig1;sig2=settings.sig2;sig3=settings.sig3;
T=settings.T;
% TODO: is it ok that non-neg only comes at the end
logp=@(pmat) -sum((yvec-pmat(:,1).*cos(2*pi*tvec.*pmat(:,3)-pmat(:,2))).^2/2,2) - ...
	(mu1-pmat(:,1)).^2/2/sig1^2 - ...
	(mu2-pmat(:,2)).^2/2/sig2^2 - ...
	(mu3-pmat(:,3)).^2/2/sig3^2 ;

X=sampleTruncatedPrior(N,settings);
Y=randn([N,3,T]).*[sig1 sig2 sig3];

if settings.run_gpu
    X=gpuArray(X);
    Y=gpuArray(Y);
end

for t=2:T
    Z=Y(:,:,t)+X;
    diffp=logp(Z)-logp(X); 
    logdiff=log(prod(normpdf(X,Z,[sig1 sig2 sig3]),2))- ...
        log(prod(normpdf(Z,X,[sig1 sig2 sig3]),2));
    alphamat=min(diffp+logdiff,0);
    rmat=log(rand(N,1));
    X=(rmat<alphamat).*(Z(:,1)>0).*(Z(:,2)>0).*(2*pi>Z(:,2)).*(Z(:,3)>0).*Y(:,:,t) + X;
end

if settings.run_gpu
    X=gather(X);
end

end

