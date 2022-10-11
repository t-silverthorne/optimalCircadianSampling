%%
estimateBayesianFIMdet(5,5,1/3)

%% Construct measurements
clear
tic
Nleft=5; % 5 
Nright=5; % 5 
tau=1/3;

Ntimes=Nleft+Nright;
mt1=linspace(0,tau,Nleft+1);
mt1=mt1(1:end-1);
mt2=linspace(tau,1,Nright+1);
mt2=mt2(1:end-1);
mt_nu=[mt1 mt2]; % construct non-uniform grid 
mt_unif=linspace(0,1,Ntimes+1); % construct uniform grid 
mt_unif=mt_unif(1:end-1);


% Construct information matrix
syms x beta1 beta2 beta3 beta4 beta5 beta6
xv=sym('x',[1 Ntimes]);
assume([x beta1 beta2 beta3 beta4 beta5 beta6],'real')
f=beta1.*sin(2*pi*x./beta5)+beta2.*cos(2*pi*x./beta5) + ...
                     beta3.*sin(2*pi*x./beta6)+beta4.*cos(2*pi*x./beta6);
gradf=gradient(f,[beta1 beta2 beta3 beta4 beta5 beta6]);
GF=subs(gradf,x,xv);
GFfun=matlabFunction(GF);
M=matlabFunction(GF*GF');

nreps=30000;
beta=rand(nreps,4);% first four parameters are unif [0,1]
beta(:,1:4)=1+beta(:,1:4);
%beta(:,1)=10;
%beta(:,2)=10;
accepted=NaN(nreps,2);
ind=1;
while ind<nreps
% final two need rejection sampling
samp=rand(floor(nreps*1.1),2);
acceptedloc=samp(abs(samp(:,1)-samp(:,2))>.2,:);
accepted(ind:ind+size(acceptedloc,1)-1,:)=acceptedloc;
ind=ind+size(acceptedloc,1);
end
beta(:,5:6)=accepted(1:nreps,:);


beta=num2cell(beta);
mtc_unif=num2cell(mt_unif);
mtc_nu=num2cell(mt_nu);

nu_vals=arrayfun(@(ind) log(det(M(beta{ind,:},mtc_nu{:}))),1:nreps);
[~,~,nuci]=normfit(nu_vals);
unif_vals=arrayfun(@(ind) log(det(M(beta{ind,:},mtc_unif{:}))),1:nreps);
[~,~,unifci]=normfit(unif_vals);
fprintf('nu     CI: [%2.3f,%2.3f]\n',nuci(1),nuci(2)); % todo: add sample sdev
fprintf('unif   CI: [%2.3f,%2.3f]\n\n',unifci(1),unifci(2));
%%
nu_val=mean(nu_vals);
unif_val=mean(unif_vals);
fprintf('nu:    %d\n',nu_val);
fprintf('unif:  %d\n\n',unif_val);
toc
%det(GF*GF')
%double(subs(subs(GF*GF',[beta1 beta2 beta3 beta4 beta5 beta6], [1 2 3 4 5 6]),xv,rand(1,7)));
% compute determinant numerically


% ind=3
% log(det(M(beta{ind,:},mtc_nu{:})))
% log(woodburyDet(GFfun(beta{ind,:},mtc_nu{:})))
function detM=woodburyDet(A)
d=size(A,1);nconv=size(A,2);% need nconv>d
M=A(:,1:d)*A(:,1:d)';
detM=det(M);
ii=d+1;
while ii<=nconv
    detM = (1+ A(:,ii)'*(M\A(:,ii)))*detM;
    M=M+A(:,ii)*A(:,ii)';
    ii=ii+1;
end
end

function detM=choleskyDet(M)
detM=prod(diag(chol(M)))^2;
end




