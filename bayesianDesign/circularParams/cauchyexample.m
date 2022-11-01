% GPU
tic
N=1e4; % desired number of samples
T=50; % length of Markov chain
k=3;
sig=2;
p=@(x) prod(1./(1+((x)/sig).^2),2);
X=randn([N,k,T])+3;
Y=randn([N,k,T]);

run_gpu=true;

if run_gpu
    X=gpuArray(X);
    Y=gpuArray(Y);
end

for t=2:T
    Z=Y(:,:,t)+X(:,:,t-1);
    diffp=log(p(Z))-log(p(X(:,:,t-1)));
    logdiff=log(prod(normpdf(X(:,:,t-1),Z),2))- ...
        log(prod(normpdf(Z,X(:,:,t-1)),2));
    alphamat=min(diffp+logdiff,0);
    rmat=log(rand(N,1));
    X(:,:,t)= (rmat<alphamat).*Y(:,:,t) + X(:,:,t-1);
end
if run_gpu
    X=gather(X);
end
toc
clf

histogram(X(:,1,end),80,'Normalization','pdf')
xv=-8:.01:8;
hold on
plot(xv,1./(1+(xv/sig).^2)/pi/sig,'linewidth',3)
ylim([0 0.5])
xlim([-20 20])