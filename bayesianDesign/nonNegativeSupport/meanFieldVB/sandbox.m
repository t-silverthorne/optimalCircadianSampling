clear
rng(12)
Nmeas=16
t_num=sort(rand(1,Nmeas));
y_num=rand(1,Nmeas);

settings.muA=10;
settings.sigA=0.1;
settings.muB=10;
settings.sigB=0.1;
settings.muf=3;
settings.sigf=1;
settings.Lcut=15;
tic
[qA,qB,qf]=updateVariationalBayes(t_num,y_num,settings);
toc

tiledlayout(3,1)
fvals=0:.01:settings.Lcut;
Avals=-10:.01:10;
Bvals=-10:.01:10;
nexttile(1)
plot(Avals,qA(Avals))
nexttile(2)
plot(Bvals,qB(Bvals))
nexttile(3)
plot(fvals,qf(fvals))
%%
% normalization factor for qf
clear all
syms muf sigf f


function [x,w] = gauss(N)
% Golub-Welsh to construct Gaussian quadrature nodes
beta = .5./sqrt(1-(2*(1:N-1)).^(-2)); T = diag(beta,1) + diag(beta,-1);
[V,D] = eig(T);
x = diag(D); [x,i] = sort(x);
w = 2*V(1,i).^2;
end

function [qA,qB,qf]=updateVariationalBayes(t_num,y_num,settings)
Nmeas=length(t_num);
Nquad=1000;
Lcut=settings.Lcut;

syms t f A B y
assume([t f A B y],'real')
tv=sym('t',[1 Nmeas]);
yv=sym('y',[1 Nmeas]);
t_numc=num2cell(t_num);
y_numc=num2cell(y_num);

% initial parameters
muA=settings.muA;
sigA=settings.sigA;
muB=settings.muB;
sigB=settings.sigB;
muf=settings.muf;
sigf=settings.sigf;

Zf=(2^(1/2)*pi^(1/2)*(erf((2^(1/2)*muf*(1/sigf^2)^(1/2))/2) + 1)*(sigf^2)^(1/2))/2;
% assume([muf sigf f],'real')
% int(exp(-(f-muf)^2/2/sigf^2),f,0,Inf)

% initial distributions
qf=@(f) (f>0).*exp(-(f-muf).^2/2/sigf^2)/Zf; % truncated Gaussian
qA=@(A) exp(-(A-muA).^2/2/sigA^2)/sqrt(2*pi)/sigA; % N(muA,sigA^2)
qB=@(B) exp(-(B-muB).^2/2/sigB^2)/sqrt(2*pi)/sigB; % N(muB,sigB^2)

% use Gaussian quadrature for qf(f) integrals
[xf,wf]=gauss(Nquad);

fprintf('\n')
pvec=[muA sigA muB sigB];
disp(pvec)

% pre compute all terms that are insensitive to shape of distribution
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

for i=1:10

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
    S1f=matlabFunction(sigA^2*sum(subs(cos(2*pi*f*t).^2,t,t_num)));
    S2f=matlabFunction(sigB^2*sum(subs(sin(2*pi*f*t).^2,t,t_num)));
    S3f=matlabFunction(2*muA*muB*subs(sin(2*pi*f*t),t,t_num)*subs(cos(2*pi*f*t),t,t_num)');
    S4f=matlabFunction(2*muA*y_num*subs(cos(2*pi*f*t),t,t_num)');
    S5f=matlabFunction(2*muB*y_num*subs(sin(2*pi*f*t),t,t_num)');
    
    % update parameters of qA, qB, qf density functions
    qf=@(f) (f>0).*exp(-(f-muf).^2/2/sigf^2).*exp(- (S1f(f) + S2f(f) + S3f(f) - S4f(f) - S5f(f))/2);
    Zf=wf*qf(xf*Lcut/2+Lcut/2)*Lcut/2; % might be safer way of doing normalization
    qf=@(f) qf(f)/Zf;
    
    muA_old=muA;
    sigA_old=sigA;
    muB_old=muB;
    sigB_old=sigB;
    
    muA=(S3A/2 - S2A/2 + muA/(2*sigA^2))/(S1A/2 + 1/(2*sigA^2));
    sigA=1/(S1A + 1/sigA^2)^(1/2);
    qA=@(A) exp(-(A-muA).^2/2/sigA^2)/sqrt(2*pi)/sigA;
    
    muB=(S3B/2 - S2B/2 + muB/(2*sigB^2))/(S1B/2 + 1/(2*sigB^2));
    sigB=1/(S1B + 1/sigB^2)^(1/2);
    qB=@(B) exp(-(B-muB).^2/2/sigB^2)/sqrt(2*pi)/sigB; % N(muB,sigB^2)
    
    %fprintf('%f\n',norm(pvec-[muA sigA muB sigB],2))
    pvec=[muA sigA muB sigB];

end
end