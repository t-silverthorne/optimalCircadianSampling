clear
rng(2)
use_tolerance=true;

% options for variational inference
NS    = 1e4; % number of samples
d     = 3;   % dimension of param space (3*num harmonics in cosinor)
NG    = 2;   % number of Gaussians to use in posterior

% specify true parameters for model
A0    = 3.5;% cos coeff
B0    = 3.5;% sin coeff
f0    = 1.8;% frequency
fmax  = 10; % max freq in posterior (used in repar)

% generate data from model
Nmeas = 10; % number of measurements
tobs  = linspace(0,1,Nmeas+1); % measurement times
tobs  = tobs(1:end-1);
yobs  = A0*cos(2*pi*f0*tobs)+B0*sin(2*pi*f0*tobs)+randn(1,numel(tobs)); % observations
tobs  = tobs'; 
yobs  = yobs';

% parameters for prior
sigA1 = 5;
sigB1 = 5;
sigT1 = 5;
Sigmat_prior=diag([sigA1,sigB1,sigT1]);
muvec_prior=[0;0;atanh(2/fmax-1)];

% randomly generated parameters for posterior
useGPU=true; % use GPU for generating random numbers and linear algebra
Lvec  = repmat(vech(eye(d)),1,NG); % vector of entries for Cholesky factors of covariance matrices
muvec = [rand(1,NG);rand(1,NG); atanh(2*linspace(1e-1,5,NG)/fmax-1)]

%rand(d,NG); % means of Gaussians in posterior  repmat([A0; B0; atanh(2*f0/fmax-1)],1,NG)
wvec  = 0.5*ones(NG,1);%rand(NG,1); wvec = wvec.^2/sum(wvec.^2); % weigths of Gaussians in posterior


% derivatives of prior and likelihood function

% syms S1 S2 S3 fmax tobs 
% 
% eta=S1*cos(tobs.*pi*fmax*(1+tanh(S3)))+S2*sin(tobs.*pi*fmax*(1+tanh(S3)))
% diff(eta,S3)
%%
% PENALTY METHOD 10\tanh\left(\frac{x}{10}\right)^{8}-10
jac1=@(S,tobs) cos(pi*fmax*pagemtimes(tobs,tanh(S(3,:,:)) + 1));
jac2=@(S,tobs) sin(pi*fmax*pagemtimes(tobs,tanh(S(3,:,:)) + 1));
%jac3OLD=@(S,tobs) S(1,:,:).*fmax.*tobs.*pi.*sin(pi*fmax*(tobs.*tanh(S(3,:,:)) + 1)).*(tanh(S(3,:,:)).^2 - 1)- ...
%               S(2,:,:).*fmax.*tobs.*pi.*cos(pi*fmax*tobs.*(tanh(S(3,:,:)) + 1)).*(tanh(S(3,:,:)).^2 - 1);
jac3=@(S,tobs) S(1,:,:)*fmax.*tobs.*pi.*sin(pi*fmax*tobs.*(tanh(S(3,:,:)) + 1)).*(tanh(S(3,:,:)).^2 - 1) - ...
            S(2,:,:)*fmax.*tobs.*pi.*cos(pi*fmax*tobs.*(tanh(S(3,:,:)) + 1)).*(tanh(S(3,:,:)).^2 - 1);

eta=@(S,tobs) S(1,:,:).*cos(pi*fmax*(1+tanh(S(3,:,:)).*tobs))+...
    S(2,:,:).*sin(pi*fmax*(1+tanh(S(3,:,:))).*tobs);
jaceta=@(S,tobs)[jac1(S,tobs) jac2(S,tobs) jac3(S,tobs)];

difflogprior=@(S) -pagemldivide(Sigmat_prior,S-muvec_prior);

diffloglhood=@(S,tobs,yobs) pagetranspose(sum((yobs-eta(S,tobs)).*jaceta(S,tobs),1));

logprior=@(S) -d/2*log(2*pi)-0.5*log(det(Sigmat_prior)) -0.5*pagemtimes(pagetranspose(S-muvec_prior),pagemldivide(Sigmat_prior,S-muvec_prior));
loglhood=@(S,tobs,yobs)  -0.5*log(2*pi) - 0.5*pagemtimes(pagetranspose(yobs-eta(S,tobs)),yobs-eta(S,tobs));

% stochastic gradient ascent for parameter optimization
beta1=.8; beta2=0.8; eps0=1e-2;tau=2000;tW=30;maxPat=1000;
niter=5e3;
LBnow=0;
for itervar=1:niter
if mod(itervar,50)==0
    disp(itervar)
    fprintf('%d\n',LBnow)
end

avec = cumsum(wvec.^2/sum(wvec.^2));
avecm1= [0; avec(1:end-1)];

switch NG
    case 2
        Sigm1=ivh(Lvec(:,1))*ivh(Lvec(:,1))';
        if isnan(Lvec)
            disp('sigm1')
        end
        Sigm2=ivh(Lvec(:,2))*ivh(Lvec(:,2))';
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


    case 3
        Sigm1=ivh(Lvec(:,1))*ivh(Lvec(:,1))';
        Sigm2=ivh(Lvec(:,2))*ivh(Lvec(:,2))';
        Sigm3=ivh(Lvec(:,3))*ivh(Lvec(:,3))';
        qlambda=    @(S) wvec(1).^2/sum(wvec.^2)*(sqrt(det(Sigm1))^(-1)*(2*pi)^(-d/2)*...
                                                   exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,1)),pagemldivide(Sigm1,S-muvec(:,1))) )) + ...
                                 wvec(2).^2/sum(wvec.^2)*(sqrt(det(Sigm2))^(-1)*(2*pi)^(-d/2)*...
                                                   exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,2)),pagemldivide(Sigm2,S-muvec(:,2))) )) + ...
                                 wvec(3).^2/sum(wvec.^2)*(sqrt(det(Sigm3))^(-1)*(2*pi)^(-d/2)*...
                                                   exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,3)),pagemldivide(Sigm3,S-muvec(:,3))) ));
        diffqlambda=@(S) wvec(1).^2/sum(wvec.^2)*(sqrt(det(Sigm1))^(-1)*(2*pi)^(-d/2)*...
                                           pagemldivide(Sigm1,(muvec(:,1)-S)).*...
                                           exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,1)),pagemldivide(Sigm1,S-muvec(:,1))) )) + ...
                         wvec(2).^2/sum(wvec.^2)*(sqrt(det(Sigm2))^(-1)*(2*pi)^(-d/2)*...
                                           pagemldivide(Sigm2,(muvec(:,2)-S)).*...
                                           exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,2)),pagemldivide(Sigm2,S-muvec(:,2))) )) + ...
                         wvec(3).^2/sum(wvec.^2)*(sqrt(det(Sigm3))^(-1)*(2*pi)^(-d/2)*...
                                           pagemldivide(Sigm3,(muvec(:,3)-S)).*...
                                           exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,3)),pagemldivide(Sigm3,S-muvec(:,3))) ));
    case 4
        Sigm1=ivh(Lvec(:,1))*ivh(Lvec(:,1))';
        Sigm2=ivh(Lvec(:,2))*ivh(Lvec(:,2))';
        Sigm3=ivh(Lvec(:,3))*ivh(Lvec(:,3))';
        Sigm4=ivh(Lvec(:,4))*ivh(Lvec(:,4))';
        qlambda=    @(S) wvec(1).^2/sum(wvec.^2)*(sqrt(det(Sigm1))^(-1)*(2*pi)^(-d/2)*...
                                                   exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,1)),pagemldivide(Sigm1,S-muvec(:,1))) )) + ...
                                 wvec(2).^2/sum(wvec.^2)*(sqrt(det(Sigm2))^(-1)*(2*pi)^(-d/2)*...
                                                   exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,2)),pagemldivide(Sigm2,S-muvec(:,2))) )) + ...
                                 wvec(3).^2/sum(wvec.^2)*(sqrt(det(Sigm3))^(-1)*(2*pi)^(-d/2)*...
                                                   exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,3)),pagemldivide(Sigm3,S-muvec(:,3))) )) + ...
                                 wvec(4).^2/sum(wvec.^2)*(sqrt(det(Sigm4))^(-1)*(2*pi)^(-d/2)*...
                                                   exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,4)),pagemldivide(Sigm4,S-muvec(:,4))) ));
        diffqlambda=@(S) wvec(1).^2/sum(wvec.^2)*(sqrt(det(Sigm1))^(-1)*(2*pi)^(-d/2)*...
                                           pagemldivide(Sigm1,(muvec(:,1)-S)).*...
                                           exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,1)),pagemldivide(Sigm1,S-muvec(:,1))) )) + ...
                         wvec(2).^2/sum(wvec.^2)*(sqrt(det(Sigm2))^(-1)*(2*pi)^(-d/2)*...
                                           pagemldivide(Sigm2,(muvec(:,2)-S)).*...
                                           exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,2)),pagemldivide(Sigm2,S-muvec(:,2))) )) + ...
                         wvec(3).^2/sum(wvec.^2)*(sqrt(det(Sigm3))^(-1)*(2*pi)^(-d/2)*...
                                           pagemldivide(Sigm3,(muvec(:,3)-S)).*...
                                           exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,3)),pagemldivide(Sigm3,S-muvec(:,3))) )) + ...
                         wvec(4).^2/sum(wvec.^2)*(sqrt(det(Sigm4))^(-1)*(2*pi)^(-d/2)*...
                                           pagemldivide(Sigm4,(muvec(:,4)-S)).*...
                                           exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,4)),pagemldivide(Sigm4,S-muvec(:,4))) ));

                

end

switch use_tolerance
    case true
        difflogqlambda=@(S) feval_with_tol(diffqlambda,qlambda,S);
    case false
        difflogqlambda=@(S)  diffqlambda(S)./qlambda(S);
end


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
%%

cc=counts;

dMU_mean_mat= zeros(d,NG);
dL_mean_mat = zeros(d*(d+1)/2,NG);
dW_mean_mat = zeros(NG,1);

dxi_dw=@(i,l) get_dadw(wvec,i-1,l) + get_dadw(wvec,i,l) ;

LBnow=0;
for ii=1:NG
    numU=cc(ii);
    if numU>0
        dxidw_l=arrayfun(@(ind) dxi_dw(ii,ind),1:NG);
    
        switch useGPU
            case true
                epsloc    = randn(d,1,numU,'gpuArray');
            case false
                epsloc    = randn(d,1,numU);
        end
        Sloc      = muvec(:,ii)+pagemtimes(ivh(Lvec(:,ii)),epsloc);
        
%         histogram(Sloc(1,:,:))
%         xlim([-20 20])
%         pause(0.1)
        diffhlambdaSloc=diffhlambda(Sloc);
        dMU       = diffhlambdaSloc; % derivative wrt mu
        dMU_mean  = sum(reshape(dMU,[d,numU]),2)/NS;
        
        Lder1     = pagemtimes(dMU,pagetranspose(epsloc));
        Lder2     = pagefun(@tril,Lder1);
        dL_d2     = reshape(Lder2,[d^2,numU]);
        dL_d2     = sum(dL_d2,2)/NS;
        dL_mean   = vech(reshape(dL_d2,[d d]));
        
        
        dW_mean=sum(pagemtimes(pagetranspose(dxidw_l.*Sloc),diffhlambdaSloc),3)/NS;
        if isnan(dW_mean)
            disp("stop")
        end
        dW_mean_mat = dW_mean_mat+dW_mean;
        dMU_mean_mat(:,ii)= dMU_mean;
        dL_mean_mat(:,ii) = dL_mean;

        LBnow=LBnow+sum(LB(Sloc,tobs,yobs))/NS;
    end
end

vMU= dMU_mean_mat.^2;
vL = dL_mean_mat.^2;
vW = dW_mean_mat.^2;

if itervar==1
    dMUbar=dMU_mean_mat
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

if itervar>0
    muvec=muvec+alphat*dMUbar./sqrt(vMUbar);
    Lvec=Lvec+alphat*dLbar./sqrt(vLbar);
    wvec=wvec+alphat*dWbar./sqrt(vWbar);
else
    muvec=muvec+alphat*dMUbar;
    Lvec=Lvec+alphat*dLbar;
    wvec=wvec+alphat*dWbar;
end

if sum(isnan(muvec))>0
    disp('bad muvec')
end

if sum(isnan(Lvec))>0
    disp('bad Lvec')
end


if sum(isnan(wvec))>0
    disp('bad wvec')
end



end

muvec
Lvec
wvec
%%
clear

function da_idw_l = get_dadw(wvec,i,l)
if i==0 || i==length(wvec)
    da_idw_l=0;
elseif l<=i
    da_idw_l=2*wvec(l)*(1/sum(wvec.^2) - sum(wvec(1:i).^2)/sum(wvec.^2));
else
    da_idw_l=-2*wvec(l)*sum(wvec(1:i).^2)/sum(wvec.^2);
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

function  y = feval_with_tol(fnum,fdenom,S)
tol         = 1e-16;
y           = gpuArray(zeros(size(S,1),size(S,2),size(S,3)));
num         = fnum(S);
denom       = fdenom(S);
flag        = abs(denom)<tol;
y(:,:,~flag)= num(:,:,~flag)./denom(:,:,~flag);
end