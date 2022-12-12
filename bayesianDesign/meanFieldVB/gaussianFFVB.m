% assume posterior is sum of mutlivariate normal distributions
% syms d N
% simplify(solve(N-d*(d+1)/2,d))
% f=@(N) (8*N + 1)^(1/2)/2 - 1/2
% f()
clear
rng(2)
% GPU random variable generation
NS=1e3;% number of samples
d=3;
Ngauss=2;

A0=2.5;
B0=2.5;
f0=2.1;
fmax=15;
Nmeas=12;

tobs=linspace(0,1,Nmeas+1);
tobs=tobs(1:end-1);
yobs=A0*cos(2*pi*f0*tobs)+B0*sin(2*pi*f0*tobs) + randn(1,numel(tobs));
tobs=tobs';
yobs=yobs';
useGPU=true;
Lvec = rand(d*(d+1)/2,Ngauss);
muvec =rand(d,Ngauss);
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



jac1=@(S,tobs) cos(pi*fmax*(pagemtimes(tobs,tanh(S(3,:,:))) + 1));
jac2=@(S,tobs) sin(pi*fmax*pagemtimes(tobs,(tanh(S(3,:,:))) + 1));
jac3=@(S,tobs) S(1,:,:).*fmax.*tobs.*pi.*sin(pi*fmax*(tobs.*tanh(S(3,:,:)) + 1)).*(tanh(S(3,:,:)).^2 - 1)- ...
               S(2,:,:).*fmax.*tobs.*pi.*cos(pi*fmax*tobs.*(tanh(S(3,:,:)) + 1)).*(tanh(S(3,:,:)).^2 - 1);
jaceta=@(S,tobs)[jac1(S,tobs) jac2(S,tobs) jac3(S,tobs)];

difflogprior=@(S) -pagemldivide(Sigmat_prior,S-muvec_prior);
diffloglhood=@(S,tobs,yobs) pagemtimes(pagetranspose(jaceta(S,tobs)),yobs-eta(S,tobs));

logprior=@(S) -d/2*log(2*pi)-0.5*log(det(Sigmat_prior)) -0.5*pagemtimes(pagetranspose(S-muvec_prior),pagemldivide(Sigmat_prior,S-muvec_prior));
loglhood=@(S,tobs,yobs)  -0.5*log(2*pi) - 0.5*pagemtimes(pagetranspose(yobs-eta(S,tobs)),yobs-eta(S,tobs));
%% derivative of hlambda and g
beta1=.1; beta2=0.1; eps0=.1;tau=2000;tW=30;maxPat=1000;

niter=1e3;

for itervar=1:niter

if mod(itervar,100)==0
    disp(itervar)
    fprintf('%d\n',LBnow)
end
avec = cumsum(wvec.^2/sum(wvec.^2));
avecm1= [0; avec(1:end-1)];


Sigm1=ivh(Lvec(:,1))*ivh(Lvec(:,1))';
Sigm2=ivh(Lvec(:,2))*ivh(Lvec(:,2))';
%Lmat2=ivh(Lvec(:,3));

qlambda=    @(S) wvec(1).^2/sum(wvec.^2)*(sqrt(det(Sigm1))^(-1)*(2*pi)^(-d/2)*...
                                   exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,1)),pagemldivide(Sigm1,S-muvec(:,1))) )) + ...
                 wvec(2).^2/sum(wvec.^2)*(sqrt(det(Sigm2))^(-1)*(2*pi)^(-d/2)*...
                                   exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,2)),pagemldivide(Sigm2,S-muvec(:,2))) ));


diffqlambda=@(S) wvec(1).^2/sum(wvec.^2)*(sqrt(det(Sigm1))^(-1)*(2*pi)^(-d/2)*...
                                   pagemldivide(Sigm1,(muvec(:,1)-S)).*...
                                   exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,1)),pagemldivide(Sigm1,S-muvec(:,1))) )) + ...
                 wvec(2).^2/sum(wvec.^2)*(sqrt(det(Sigm2))^(-1)*(2*pi)^(-d/2)*...
                                   pagemldivide(Sigm2,(muvec(:,2)-S)).*...
                                   exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,2)),pagemldivide(Sigm2,S-muvec(:,2))) ));


difflogqlambda=@(S) diffqlambda(S)./qlambda(S);

diffhlambda=@(S) difflogprior(S)+diffloglhood(S,tobs,yobs)-difflogqlambda(S);

LB=@(S,tobs,yobs) logprior(S) + loglhood(S,tobs,yobs) - log(qlambda(S));

switch useGPU
    case true
        u    = rand(1,NS,'gpuArray');
    case false
        u    = rand(1,NS);
end

xi=avecm1<u & u<avec;
counts=gather(sum(xi,2));

cc=counts;
% generate samples
% switch useGPU
%     case true
%         eps=randn(d,1,NS,'gpuArray');
%     case false
%         eps=randn(d,1,NS);
% end

dMU_mean_mat= zeros(d,Ngauss);
dL_mean_mat = zeros(d*(d+1)/2,Ngauss);
dW_mean_mat = zeros(Ngauss,1);


dxi_dw=@(i,l) get_dadw(wvec,i-1,l) - get_dadw(wvec,i,l) ;

LBnow=0;
for ii=1:Ngauss
    numU=cc(ii);
    if numU>0
        dxidw_l=arrayfun(@(ind) dxi_dw(ii,ind),1:Ngauss);
    
        switch useGPU
            case true
                epsloc    = randn(d,1,numU,'gpuArray');
            case false
                epsloc    = randn(d,1,numU);
        end
        Sloc      = muvec(:,ii)+pagemtimes(ivh(Lvec(:,ii)),epsloc);
        diffhlambda_epsloc=diffhlambda(epsloc);
        dMU       = diffhlambda_epsloc; % derivative wrt mu
        dMU_mean  = sum(reshape(dMU,[d,numU]),2)/numU;
        
        Lder1     = pagemtimes(dMU,pagetranspose(epsloc));
        Lder2     = pagefun(@tril,Lder1);
        dL_d2     = reshape(Lder2,[d^2,numU]);
        dL_d2     = sum(dL_d2,2)/numU;
        dL_mean   = vech(reshape(dL_d2,[d d]));
            
        dW_mean=sum(pagemtimes(pagetranspose(dxidw_l.*Sloc),diffhlambda_epsloc),3)/numU;
        if isnan(dW_mean)
            itervar
            disp('bad dW_mean')
        end
        dW_mean_mat   =  dW_mean_mat+dW_mean;
        dMU_mean_mat(:,ii)= dMU_mean;
        dL_mean_mat(:,ii) = dL_mean;

        LBnow=LBnow+sum(LB(Sloc,tobs,yobs))/numU;
    end
end

vMU= dMU_mean_mat.^2;
vL = dL_mean_mat.^2;
vW = dW_mean_mat.^2;

if itervar==1
    dMUbar=dMU_mean_mat;
    dLbar =dL_mean_mat;
    dWbar =dW_mean_mat;
    vMUbar=vMU;vLbar=vL;vWbar=vW;
else
    dMUbar= beta1*dMUbar+(1-beta1)*dMU_mean_mat;
    dLbar = beta1*dLbar+(1-beta1)*dL_mean_mat;
    dWbar = beta1*dWbar+(1-beta1)*dW_mean_mat;
    vMUbar= beta2*vMUbar+(1-beta2)*vMU;
    vLbar = beta2*vLbar+(1-beta2)*vL;
    vWbar = beta2*vWbar+(1-beta2)*vW;
end

alphat=eps0;
muvec=muvec-alphat*dMUbar./sqrt(vMUbar);
Lvec=Lvec-alphat*dLbar./sqrt(vLbar);
wvec=wvec-alphat*dWbar./sqrt(vWbar);
if sum(isnan(wvec))>0
    disp('bad wvec')
end



end
%%
muvec
Lvec
wvec
%%

%%

%%
% for ii=2:Ngauss
%     i0=cc(ii-1)+1;
%     Lmatnow=ivh(Lvec(:,ii));
%     switch useGPU
%         case true
%             S(:,i0:cc(ii))=muvec(:,ii)+Lmatnow*randn(d,cc(ii)+1-i0,'gpuArray');
%         case false
%             S(:,i0:cc(ii))=muvec(:,ii)+Lmatnow*randn(d,cc(ii)+1-i0);
%     end
% end

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
function da_idw_l = get_dadw(wvec,i,l)
if i==0 || i==length(wvec)
    da_idw_l=0;
elseif l<=i
    da_idw_l=2*wvec(l)*(1/sum(wvec.^2) - sum(wvec(1:i).^2)/sum(wvec.^2));
else
    da_idw_l=-2*wvec(l)*sum(wvec(1:i).^2)/sum(wvec.^2);
end
end

function S=generateSamples(muvec,Lvec,wvec,NS,d,Ngauss,useGPU)
switch useGPU
    case true
        u    = rand(1,NS,'gpuArray');
    case false
        u    = rand(1,NS);
end
avec = cumsum(wvec.^2/sum(wvec.^2));
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

% function Ltall = ivhmat(Lmat)
% % apply inverse of vech map to each column then take transpose of each lower
% % triangular matrix and stack the results vertically
% d= (8*size(Lmat,1) + 1)^(1/2)/2 - 1/2;
% Ltall=NaN(d*size(Lmat,2),d);
% for ii=1:size(Lmat,2)
%     Ltall( (1+d*(ii-1)):d*ii ,:)=ivh(Lmat(:,ii));
% end
% 
% end