%% simulated data
clf
clear all
rng(12)
A0=10;
B0=7;
f0=1.2;
fmax=100;
Lcut=fmax;
T=1000;
Nsamp=1e3;
Nquad=1e3;
Nmeas=100;
muA=A0;
sigA=3;
muB=B0;
muf=0;
sigf=.05;
sigB=3;
sig1=.5;sig2=.5;sig3=.5;

t_obs=linspace(0,1,Nmeas+1);
t_obs=t_obs(1:end-1);
y_obs=A0*cos(2*pi*f0*t_obs)+B0*sin(2*pi*f0*t_obs) + randn(1,numel(t_obs));

%% MCMC
tic
X=randn(Nsamp,3);
X(:,1)=X(:,1)*sigA+muA; X(:,2)=X(:,2)*sigB+muB;
X(:,3)=fmax*rand(1,Nsamp);

pfun=@(par) -sum((y_obs - par(:,1).*cos(2*pi*par(:,3)*t_obs) ...
    - par(:,2).*sin(2*pi*par(:,3)*t_obs)   ).^2/2,2) - ...
    (muA-par(:,1)).^2/2/sigA^2 -  ...
    (muB-par(:,2)).^2/2/sigB^2 - ...
    -(log(abs(par(:,3))) - muf).^2/2/sigf^2 - par(:,3)*sigf;

Y=randn(Nsamp,3,T);
Y(:,1)=Y(:,1)*sig1; Y(:,2)=Y(:,2)*sig2; Y(:,3)=Y(:,3)*sig3;
for t=2:T
    Z=Y(:,:,t)+X;
    diffp=pfun(Z)-pfun(X); 
    logdiff=log(prod(normpdf(X,Z,[sig1 sig2 sig3]),2))- ...
        log(prod(normpdf(Z,X,[sig1 sig2 sig3]),2));
    alphamat=min(diffp+logdiff,0);
    rmat=log(rand(Nsamp,1));
    X=(rmat<alphamat).*(Z(:,3)>0).*Y(:,:,t) + X;
end
toc

tiledlayout(3,1)
ptrue=[A0 B0 f0];
for ind=1:3
    nexttile(ind)
    histogram(X(:,ind),max(120,floor(sqrt(Nsamp))),'normalization','pdf','EdgeColor','none')
    hold on
    xline(ptrue(ind))
end
%%
t=sym('t'); f=sym('f'); A=sym('A'); B=sym('B'); y=sym('y');
assume([t f A B y],'real')
tv=sym('t',[1 Nmeas]);
yv=sym('y',[1 Nmeas]);
t_numc=num2cell(t_obs);
y_numc=num2cell(y_obs);

S1A_int=cos(2*pi*f*t).^2;
S1A_int=matlabFunction(sum(subs(S1A_int,t,tv)));
S1A_int_fun=@(f) S1A_int(f,t_numc{:}) ; % reparameterize from -1 to 1

S2A_int=cos(2*pi*f*t).*sin(2*pi*f*t);
S2A_int=matlabFunction(sum(subs(S2A_int,t,tv)));
S2A_int_fun=@(f) S2A_int(f,t_numc{:}) ; % reparameterize from -1 to 1

S3A_int=cos(2*pi*f*t).*y;
SS=subs(S3A_int,t,tv);
S3A_int=matlabFunction(sum(arrayfun( @(ind) subs(SS(ind),y,yv(ind)),1:Nmeas)));
S3A_int_fun=@(f) S3A_int(f,t_numc{:},y_numc{:});

S1B_int=sin(2*pi*f*t).^2;
S1B_int=matlabFunction(sum(subs(S1B_int,t,tv)));
S1B_int_fun=@(f) S1B_int(f,t_numc{:}) ; % reparameterize from -1 to 1

S3B_int=sin(2*pi*f*t).*y;
SS=subs(S3B_int,t,tv);
S3B_int=matlabFunction(sum(arrayfun( @(ind) subs(SS(ind),y,yv(ind)),1:Nmeas)));
S3B_int_fun=@(f) S3B_int(f,t_numc{:},y_numc{:});

[xf,wf]=gauss(Nquad);

muAINIT=muA;
muBINIT=muB;
sigAINIT=sigA;
sigBINIT=sigB;
%%
qf=@(f) (f>0).*(fmax>f);
for i=1:25
    % get qA terms
    S1A_num=wf*(S1A_int_fun(xf*Lcut/2+Lcut/2).*qf(xf*Lcut/2+Lcut/2)*Lcut/2);
    S2A_num=wf*(S2A_int_fun(xf*Lcut/2+Lcut/2).*qf(xf*Lcut/2+Lcut/2)*Lcut/2);
    S3A_num=wf*(S3A_int_fun(xf*Lcut/2+Lcut/2).*qf(xf*Lcut/2+Lcut/2)*Lcut/2);
    S1A=S1A_num; S2A=S2A_num*muB; S3A=S3A_num; % update with exact part
   
    % get qB terms
    S1B_num=wf*(S1B_int_fun(xf*Lcut/2+Lcut/2).*qf(xf*Lcut/2+Lcut/2)*Lcut/2);
    S2B_num=S2A_num;
    S3B_num=wf*(S3B_int_fun(xf*Lcut/2+Lcut/2).*qf(xf*Lcut/2+Lcut/2)*Lcut/2);
    S1B=S1B_num; S2B=S2B_num*muA; S3B=S3B_num;
    
    % get qf terms
    S1f=matlabFunction((muA^2+sigA^2)*sum(subs(cos(2*pi*f*t).^2,t,t_obs)));
    S2f=matlabFunction((muB^2+sigB^2)*sum(subs(sin(2*pi*f*t).^2,t,t_obs)));
    S3fpre=matlabFunction(subs(sin(2*pi*f*t),t,t_obs)*subs(cos(2*pi*f*t),t,t_obs)');
    S4fpre=matlabFunction(y_obs*subs(cos(2*pi*f*t),t,t_obs)');
    S5fpre=matlabFunction(y_obs*subs(sin(2*pi*f*t),t,t_obs)');
    S3f=@(f) 2*muA*muB*S3fpre(f);
    S4f=@(f) 2*muA*S4fpre(f);
    S5f=@(f) 2*muB*S5fpre(f);

    % update parameters
    muA=(S3A/2 - (S2A*muB)/2 + muAINIT/(2*sigAINIT^2))/(S1A/2 + 1/(2*sigAINIT^2));
    muB=(S3B/2 - (S2B*muA)/2 + muBINIT/(2*sigAINIT^2))/(S1B/2 + 1/(2*sigBINIT^2));
    sigA=1/(S1A + 1/sigAINIT^2)^(1/2);
    sigB=1/(S1B + 1/sigBINIT^2)^(1/2);

    % update functions
    qA=@(A) exp(-(A-muA).^2/2/sigA^2)/sqrt(2*pi)/sigA;
    qB=@(B) exp(-(B-muB).^2/2/sigB^2)/sqrt(2*pi)/sigB;
    
	qf=@(F) exp(-F/muf).*exp(- (S1f(F) + S2f(F) + S3f(F) - S4f(F) - S5f(F))/2)/fmax;
    Zf=wf*qf(xf*Lcut/2+Lcut/2)*Lcut/2; % there might be safer way of doing normalization
    qf=@(f) qf(f)/Zf;
    
    %fprintf('%f\n',norm(pvec-[muA sigA muB sigB],2))
    pvec=[muA sigA muB sigB];

end

for ind=1:3
    t=nexttile(ind);
    xv=t.XLim(1):.1:t.XLim(2);
    switch ind
        case 1
            plot(xv,qA(xv))
        case 2
            plot(xv,qB(xv))
        case 3
            plot(xv,qf(xv))
    end
end

function [x,w] = gauss(N)
	% Golub-Welsh to construct Gaussian quadrature nodes
	beta = .5./sqrt(1-(2*(1:N-1)).^(-2)); T = diag(beta,1) + diag(beta,-1);
	[V,D] = eig(T);
	x = diag(D); [x,i] = sort(x);
	w = 2*V(1,i).^2;
end
