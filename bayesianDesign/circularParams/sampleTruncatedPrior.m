function X = sampleTruncatedPrior(N,settings)
if nargin<2
	settings.model='cosinorOneFreq';
    settings.run_gpu=false;
	settings.Tprior=100;
    settings.sig1=1;
    settings.sig2=settings.sig1;
    settings.sig3=settings.sig1;
    settings.mu1=1;
    settings.mu2=0;
    settings.mu3=2;
end
T=settings.Tprior;

sig1=settings.sig1;sig2=settings.sig2;sig3=settings.sig3;
mu1=settings.mu1;mu2=settings.mu2;mu3=settings.mu3;

pfun=@(p) (p(:,1)>0).*(p(:,2)>0).*(2*pi>p(:,2)).*(p(:,3)>0).* ...
	normpdf(p,[mu1,mu2,mu3],[sig1,sig2,sig3]); % TODO check this evaluates Gaussian columnwise with correct mean

X=[mu1 mu2 mu3]+zeros(N,3,1); % initial state
Y=randn([N,3,T]).*[sig1 sig2 sig3]; % noise for Brownian motion


if settings.run_gpu
    X=gpuArray(X);
    Y=gpuArray(Y);
end

for t=2:T
    Z=Y(:,:,t)+X;
    diffp=log(pfun(Z))-log(pfun(X));
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
