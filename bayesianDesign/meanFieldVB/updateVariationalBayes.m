function [qA,qB,qf,pars,qfdom,Zf0,Zf]=updateVariationalBayes(t_num,y_num,settings)
fmax=settings.fmax;
Nmeas=length(t_num);
Nquad=100;
Lcut=settings.fmax*3;
f_upper_bound=false;
t=sym('t'); f=sym('f'); A=sym('A'); B=sym('B'); y=sym('y');
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
%%
%clear all
%syms muf sigf f
%assume([muf sigf f],'positive')
%Zf0=int(exp(-(f-muf).^2/2/sigf^2)/sigf/sqrt(2*pi),f,0,Inf)
Zf0=(2^(1/2)*sigf*pi^(1/2)*(erf((2^(1/2)*muf)/(2*sigf)) + erf((2^(1/2)*(fmax - muf))/(2*sigf))))/2/sigf/sqrt(2*pi);

%TODO: put explicit cut off on the f
% initial distributions
qf=@(f) (fmax>f).*(f>0).*exp(-(f-muf).^2/2/sigf^2)/sigf/sqrt(2*pi)/Zf0; % truncated Gaussian
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
    S1f=matlabFunction((muA^2+sigA^2)*sum(subs(cos(2*pi*f*t).^2,t,t_num)));
    S2f=matlabFunction((muB^2+sigB^2)*sum(subs(sin(2*pi*f*t).^2,t,t_num)));
    S3f=matlabFunction(2*muA*muB*subs(sin(2*pi*f*t),t,t_num)*subs(cos(2*pi*f*t),t,t_num)');
    S4f=matlabFunction(2*muA*y_num*subs(cos(2*pi*f*t),t,t_num)');
    S5f=matlabFunction(2*muB*y_num*subs(sin(2*pi*f*t),t,t_num)');
	

    % update parameters of qA, qB, qf density functions
    muA=(S3A/2 - S2A/2 + muA/(2*sigA^2))/(S1A/2 + 1/(2*sigA^2));
    sigA=1/(S1A + 1/sigA^2)^(1/2);
    qA=@(A) exp(-(A-muA).^2/2/sigA^2)/sqrt(2*pi)/sigA;

    muB=(S3B/2 - S2B/2 + muB/(2*sigB^2))/(S1B/2 + 1/(2*sigB^2));
    sigB=1/(S1B + 1/sigB^2)^(1/2);
    qB=@(B) exp(-(B-muB).^2/2/sigB^2)/sqrt(2*pi)/sigB; % N(muB,sigB^2)
    
	qf=@(f) (fmax>f).*(f>0).*exp(-(f-muf).^2/2/sigf^2 - (S1f(f) + S2f(f) + S3f(f) - S4f(f) - S5f(f))/2);
    Zf=wf*qf(xf*Lcut/2+Lcut/2)*Lcut/2; % might be safer way of doing normalization
    qf=@(f) qf(f)/Zf;
    

    %fprintf('%f\n',norm(pvec-[muA sigA muB sigB],2))
    pvec=[muA sigA muB sigB];

end

pars.muA=muA;
pars.sigA=sigA;
pars.muB=muB;
pars.sigB=sigB;

fvals=0:.01:Lcut;
C=min((S1f(fvals) + S2f(fvals) + S3f(fvals) - S4f(fvals) - S5f(fvals))/2);
qfdom=exp(-C);


function [x,w] = gauss(N)
	% Golub-Welsh to construct Gaussian quadrature nodes
	beta = .5./sqrt(1-(2*(1:N-1)).^(-2)); T = diag(beta,1) + diag(beta,-1);
	[V,D] = eig(T);
	x = diag(D); [x,i] = sort(x);
	w = 2*V(1,i).^2;
end
end
