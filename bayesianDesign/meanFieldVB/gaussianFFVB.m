% assume posterior is sum of mutlivariate normal distributions
% syms d N
% simplify(solve(N-d*(d+1)/2,d))
% f=@(N) (8*N + 1)^(1/2)/2 - 1/2
% f()
clear 
% GPU random variable generation
NS=1e1;% number of samples
d=3;
Ngauss=2;
useGPU=false;
Lvec = rand(d*(d+1)/2,Ngauss);
muvec = rand(d,Ngauss);
wvec = rand(Ngauss,1); wvec = wvec/sum(wvec);
fprintf('running\n')
%%
S=generateSamples(muvec,Lvec,wvec,NS,d,Ngauss,useGPU);
S=reshape(S,[d,1,NS]);
%%
fmax=10;
sigA1=1;
sigB1=1;
sigT1=1;
Sigmat_prior=diag([sigA1,sigB1,sigT1]);
muvec_prior=[0;0;0];

%% prior and likelihood (vectorized to act on samples)
% prior is product of Gaussians
prior=@(A1,B1,T1) exp(-0.5*([A1;B1;T1] - muvec_prior)'*(Sigmat_prior\([A1;B1;T1] - muvec_prior)))*(2*pi)^(-d/2)/sqrt(det(Sigmat_prior));
lhood=@(tobs,yobs,A1,B1,T1) exp(-0.5*sum(( A1*cos(pi*fmax*(1+tanh(T1)*tobs))+B1*sin(pi*fmax*(1+tanh(T1))*tobs) -yobs).^2))/sqrt(2*pi);
% vectorized version of prior (capable of using higher dimensional arrays)
priorvec=@(S) exp(pagemtimes(-0.5*pagetranspose(S - muvec_prior), ...
                             (pagemldivide(Sigmat_prior,(S - muvec_prior)))))*...
                  (2*pi)^(-d/2)/sqrt(det(Sigmat_prior));

eta=@(S,tobs) S(1,:,:).*cos(pi*fmax*(1+tanh(S(3,:,:)).*tobs))+...
    S(2,:,:).*sin(pi*fmax*(1+tanh(S(3,:,:))).*tobs);
lhoodvec=@(S,tobs,yobs) exp(-0.5*sum( (eta(S,tobs)-yobs).^2,1))/sqrt(2*pi);
%lhood([1;2;3],[4;5;6],1,1,2)
%lhoodvec(S,[1;2;3],[4;5;6])
%% derivatives of prior and likelihood function
tobs=[1;2;3];
yobs=[4;5;6];
jac1=@(S,tobs) cos(pi*fmax*(pagemtimes(tobs,tanh(S(3,:,:))) + 1));
jac2=@(S,tobs) sin(pi*fmax*pagemtimes(tobs,(tanh(S(3,:,:))) + 1));
jac3=@(S,tobs) S(1,:,:).*fmax.*tobs.*pi.*sin(pi*fmax*(tobs.*tanh(S(3,:,:)) + 1)).*(tanh(S(3,:,:)).^2 - 1)- ...
               S(2,:,:).*fmax.*tobs.*pi.*cos(pi*fmax*tobs.*(tanh(S(3,:,:)) + 1)).*(tanh(S(3,:,:)).^2 - 1);
jaceta=@(S,tobs)[jac1(S,tobs) jac2(S,tobs) jac3(S,tobs)];

difflogprior=@(S) -pagemldivide(Sigmat_prior,S-muvec_prior);
diffloglhood=@(S,tobs,yobs) pagemtimes(jaceta(S,tobs),yobs-eta(S,tobs));

%diffprior(S)
%difflhood(S,tobs,yobs)
%% derivative of hlambda and g
beta1=.1; beta2=0.1; eps0=1e-3;tau=2000;tW=30;maxPat=1000;

niter=1;

Sigm1=ivh(Lvec(:,1))*ivh(Lvec(:,1))';
Sigm2=ivh(Lvec(:,2))*ivh(Lvec(:,2))';
%Lmat2=ivh(Lvec(:,3));

qlambda=    @(S) wvec(1)*(sqrt(det(Sigm1))^(-1)*(2*pi)^(-d/2)*...
                                   exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,1)),pagemldivide(Sigm1,S-muvec(:,1))) )) + ...
                 wvec(2)*(sqrt(det(Sigm2))^(-1)*(2*pi)^(-d/2)*...
                                   exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,2)),pagemldivide(Sigm2,S-muvec(:,2))) ));


diffqlambda=@(S) wvec(1)*(sqrt(det(Sigm1))^(-1)*(2*pi)^(-d/2)*...
                                   pagemldivide(Sigm1,(muvec(:,1)-S)).*...
                                   exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,1)),pagemldivide(Sigm1,S-muvec(:,1))) )) + ...
                 wvec(2)*(sqrt(det(Sigm2))^(-1)*(2*pi)^(-d/2)*...
                                   pagemldivide(Sigm2,(muvec(:,2)-S)).*...
                                   exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,2)),pagemldivide(Sigm2,S-muvec(:,2))) ));


difflogqlambda=@(S) diffqlambda(S)./qlambda(S);

difflogqlambda(S)
diffhlambda=@(S) difflogprior(S)+diffloglhood(S,tobs,yobs)-difflogqlambda(S);
diffhlambda(S)

switch useGPU
    case true

        u    = rand(1,NS,'gpuArray');
    case false

        u    = rand(1,NS);
end
avec = cumsum(wvec);
avecm1= [0; avec(1:end-1)];
xi=avecm1<u & u<avec;
counts=gather(sum(xi,2));
cc=cumsum(counts);

% generate samples
switch useGPU
    case true
        eps=randn(d,1,NS,'gpuArray');
    case false
        eps=randn(d,1,NS);
end

dMU_mean_mat= NaN(d,Ngauss);
dL_mean_mat = NaN(d*(d+1)/2,Ngauss);
dA_mean_mat= zeros(Ngauss-1,1);
for ii=1:Ngauss
    numU=cc(ii);
    
    switch useGPU
        case true
            epsloc    = randn(d,1,numU,'gpuArray');
        case false
            epsloc    = randn(d,1,numU);
    end
    Sloc      = muvec(:,ii)+pagemtimes(ivh(Lvec(:,ii)),eps);
    dMU       = diffhlambda(epsloc); % derivative wrt mu
    dMU_mean  = sum(reshape(dMU,[d,numU]),2)/numU;
    
    Lder1     = pagemtimes(dMU,pagetranspose(epsloc));
    Lder2     = pagefun(@tril,Lder1);
    dL_d2     = reshape(Lder2,[d^2,numU]);
    dL_d2     = sum(dL_d2,2)/numU;
    dL_mean   = vech(reshape(dL_d2,[d d]));
        
    fgauss1=    (sqrt(det(Sigm1))^(-1)*(2*pi)^(-d/2)*...
                            exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,1)),pagemldivide(Sigm1,S-muvec(:,1))) ));
    fgauss2=    (sqrt(det(Sigm2))^(-1)*(2*pi)^(-d/2)*...
                            exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,2)),pagemldivide(Sigm2,S-muvec(:,2))) ));
    dA_mean   =  sum(fgauss2-fgauss1)/NS;
    dMU_mean_mat(:,ii)= dMU_mean;
    dL_mean_mat(:,ii) = dL_mean;
    dA_mean_mat      = dA_mean_mat+dA_mean;
end
dMU_mean_mat

vMU= dMU_mean_mat.^2;
vL = dL_mean_mat.^2;
vA = dA_mean_mat.^2;

if niter==1
    dMUbar=dMU_mean_mat;
    dLbar =dL_mean_mat;
    dAbar =dA_mean_mat;
    vMUbar=vMU;vLbar=vL;vAbar=vA;
else
    dMUbar= beta1*dMUbar+(1-beta1)*dMU_mean_mat;
    dLbar = beta1*dLbar+(1-beta1)*dL_mean_mat;
    dAbar = beta1*dAbar+(1-beta1)*dA_mean_mat;
    vMUbar= beta1*vMUbar+(1-beta1)*vMU;
    vLbar = beta1*vLbar+(1-beta1)*vL;
    vAbar = beta1*vAbar+(1-beta1)*vA;
end

muvec=muvec+dMUbar./sqrt(vMUbar);
Lvec=Lvec+dLbar./sqrt(vLbar);
avec=avec+dAbar./sqrt(vAbar);
%%

%%




%%
for ii=2:Ngauss
    i0=cc(ii-1)+1;
    Lmatnow=ivh(Lvec(:,ii));
    switch useGPU
        case true
            S(:,i0:cc(ii))=muvec(:,ii)+Lmatnow*randn(d,cc(ii)+1-i0,'gpuArray');
        case false
            S(:,i0:cc(ii))=muvec(:,ii)+Lmatnow*randn(d,cc(ii)+1-i0);
    end
end

%Sloc=


%pagemtimes(diffhlambda(S),pagetranspose(eps))

%%
% syms S1 S2 S3 fmax tobs
% size(jacobian(S1.*cos(pi*fmax*(1+tanh(S3).*tobs))+...
%     S2.*sin(pi*fmax*(1+tanh(S3)).*tobs),[S1; S2; S3]))
% clear tobs
% tobs=[1;2;3]
% [cos(pi*fmax*(pagemtimes(tobs,tanh(S(3,:,:))) + 1)), ...
%     sin(pi*fmax*pagemtimes(tobs,(tanh(S(3,:,:))) + 1)), ...
%     S(1,:,:)*fmax*tobs*pi*sin(pi*fmax*(tobs*tanh(S(3,:,:)) + 1))*(tanh(S(3,:,:))^2 - 1) - S(2,:,:)*fmax*tobs*pi*cos(pi*fmax*tobs*(tanh(S(3,:,:)) + 1))*(tanh(S(3,:,:))^2 - 1)]
%  

function S=generateSamples(muvec,Lvec,wvec,NS,d,Ngauss,useGPU)
switch useGPU
    case true
        u    = rand(1,NS,'gpuArray');
    case false
        u    = rand(1,NS);
end
avec = cumsum(wvec);
avecm1= [0; avec(1:end-1)];
xi=avecm1<u & u<avec;
counts=gather(sum(xi,2));
cc=cumsum(counts);

switch useGPU
    case true
        S(:,1:cc(1))=muvec(:,1)+ivh(Lvec(:,1))*randn(d,cc(1),'gpuArray');
    case false
        S(:,1:cc(1))=muvec(:,1)+ivh(Lvec(:,1))*randn(d,cc(1));
end

for ii=2:Ngauss
    i0=cc(ii-1)+1;
    Lmatnow=ivh(Lvec(:,ii));
    switch useGPU
        case true
            S(:,i0:cc(ii))=muvec(:,ii)+Lmatnow*randn(d,cc(ii)+1-i0,'gpuArray');
        case false
            S(:,i0:cc(ii))=muvec(:,ii)+Lmatnow*randn(d,cc(ii)+1-i0);
    end
end
end
   
function Avech=vech(A)
d=size(A,1);
Avech=[];
for ii=1:size(A,2)
    for jj=ii:d
        Avech(end+1)=A(jj,ii);
    end
end
Avech=reshape(Avech,[d*(d+1)/2,1]);
end


function Lmat = ivh(Lvec)
% inverse of vech map
d= (8*length(Lvec) + 1)^(1/2)/2 - 1/2;
Lmat = zeros(d,d);
i0=1;
for i=1:d
    Lmat(i:d,i)=Lvec(i0:i0+d-i);
    i0=i0+d+1-i;
end
end

function Ltall = ivhmat(Lmat)
% apply inverse of vech map to each column then take transpose of each lower
% triangular matrix and stack the results vertically
d= (8*size(Lmat,1) + 1)^(1/2)/2 - 1/2;
Ltall=NaN(d*size(Lmat,2),d);
for ii=1:size(Lmat,2)
    Ltall( (1+d*(ii-1)):d*ii ,:)=ivh(Lmat(:,ii));
end

end