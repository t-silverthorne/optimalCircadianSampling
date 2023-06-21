Npar=1e7;
pmat=rand(Npar,3);
Nmeas=20;
tobs=rand(1,Nmeas);
tobsg=gpuArray(tobs);
yobs=rand(1,Nmeas);
yobsg=gpuArray(yobs);
mutheta1=0;
sigtheta1=1;
muamp1=10;
sigamp1=2;
muT1=0.5;
sigT1=0.01;

pmatg=gpuArray(pmat);
% vectorize over observations, arrayfun over pmat
logpV1=@(pmat,tvec,yvec,muamp1,sigamp1, ...
    mutheta1,sigtheta1, ...
    muT1,sigT1) -sum((yvec-pmat(1).*cos(2*pi*tvec./pmat(3)-pmat(2))).^2/2) + ...
         (muamp1-pmat(1)).^2/2/sigamp1^2 + ...
             (mutheta1-pmat(2)).^2/2/sigtheta1^2 + ...
                (muT1-pmat(3)).^2/2/sigT1^2 ;

% vectorize over pmat, arrayfun over observations
logpV2=@(pmat,tvec,yvec,muamp1,sigamp1, ...
    mutheta1,sigtheta1, ...
    muT1,sigT1) -sum((yvec-pmat(:,1).*cos(2*pi*tvec./pmat(:,3)-pmat(:,2))).^2/2,2) + ...
         (muamp1-pmat(:,1)).^2/2/sigamp1^2 + ...
            (mutheta1-pmat(:,2)).^2/2/sigtheta1^2 + ...
              (muT1-pmat(:,3)).^2/2/sigT1^2 ;

%% bad vectorize
% tic
% res1=gather(arrayfun(@(ind) logpV1(pmat(ind,:),tobs,yobs,muamp1,sigamp1,mutheta1,sigtheta1,muT1,sigT1), ...
%             1:size(pmat,1)));
% toc

% good vectorize no GPU
tic
res2=gather(logpV2(pmat,tobs,yobs,muamp1,sigamp1, ...
    mutheta1,sigtheta1, ...
    muT1,sigT1));
toc

% make pmat a gpu array
tic
res2=gather(logpV2(pmatg,tobs,yobs,muamp1,sigamp1, ...
    mutheta1,sigtheta1, ...
    muT1,sigT1));
toc
% make observation times also gpu array
tic
res2=gather(logpV2(pmatg,tobsg,yobsg,muamp1,sigamp1, ...
    mutheta1,sigtheta1, ...
    muT1,sigT1));
toc

%B=logpV2(pmat,tobs(1),yobs(1),muamp1,sigamp1,mutheta1,sigtheta1,muT1,sigT1) 
%%
A=rand(size(pmat,2),1e3);
Ag=gpuArray(A);
tic
pmat*A;
toc
tic
pmat*Ag;
toc