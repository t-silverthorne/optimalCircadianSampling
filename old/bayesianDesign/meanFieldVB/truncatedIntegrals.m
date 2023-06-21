%% Truncated integrals
syms muf sigf f fmax
assume([muf sigf f fmax],'positive')
Zf0=int(exp(-(f-muf).^2/2/sigf^2),f,0,fmax)
%%
Zf0= (2^(1/2)*sigf*pi^(1/2)*(erf((2^(1/2)*muf)/(2*sigf)) + erf((2^(1/2)*(fmax - muf))/(2*sigf))))/2/sigf/sqrt(2*pi)

%%
%clear
syms x mu sig Z
assume([x mu sig Z],'positive')
int(x^2*exp(-(x-mu).^2/2/sig^2),x,-Inf,Inf)


%%
int(x*exp(-(x-mu).^2/2/sig^2),x,0,Inf)
meanTrunc=@(x,mu,sig,Z) ( sig^2*exp(-mu^2/(2*sig^2)) + ...
    (2^(1/2)*mu*sig*pi^(1/2))/2 + ...
    (2^(1/2)*mu*sig*pi^(1/2)*erf((2^(1/2)*mu)/(2*sig)))/2 )/sig/sqrt(2*pi)/Z;

%int(x^2*exp(-(x-mu).^2/2/sig^2),x,0,Inf)
varTrunc=@(x,mu,sig,Z) (mu*sig^2*exp(-mu^2/(2*sig^2)) + ...
    (2^(1/2)*sig*pi^(1/2)*(mu^2 + sig^2))/2 - ...
    (2^(1/2)*sig*pi^(1/2)*erfi((2^(1/2)*mu*1i)/(2*sig))*(mu^2 + sig^2)*1i)/2)/sig/sqrt(2*pi)/Z;
