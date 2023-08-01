addpath('../utils/')
clear all


Nmeas=7; % good to set small number of measurements
mt=linspace(0,1,Nmeas+1);
mt=mt(1:end-1);
Nres=1e5;

Amp=.2;
freq=2.5;
acro=2;
p.Amp=Amp;
p.freq=freq;
p.acro=acro;
X=constructReducedX(mt,freq);


% verify chisq(lambda)
csq       = cos(2*pi*freq*mt-acro);
lambdapap = Amp^2*csq*csq';

Y=Amp*cos(2*pi*freq*mt-acro)+randn(Nres,Nmeas);

ft=fit_cosinor_model(Y,mt,1/freq);
beta=Amp*[sin(2*pi*freq) cos(2*pi*freq)];

% %% bias lower bound
norm(beta,1)/sqrt(2) 
mean(ft.amplitudes)
%%
norm(beta,2)^2 + trace(inv(X'*X))-norm(beta,1)^2/2
var(ft.amplitudes)