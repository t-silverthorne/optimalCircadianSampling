clear
checkUseGPU

param.NL=4;
param.NR=4;
param.freq_true=1.4; % freq used in regression model
param.Amp=2.1;
param.acro=2;
beta=[0; param.Amp*sin(param.acro); param.Amp*cos(param.acro)];
[~,t]=getSamplingSchedules(param.NL,param.NR,0,0.25);
X=constructX(t,param);
d=length(t);n=d;
L=X*((X'*X)\X');

C1=factorial(n-2)/factorial(n);
C2=factorial(n-1)/factorial(n);

muL=NaN(d,d);
for ii=1:d
    muL(ii,ii)=C2*trace(L);
end

for ii=1:d
    for jj=1:d
        if ~(ii==jj)
            muL(ii,jj)=C1*(sum(sum(L))-trace(L));
        end
    end
end
muL;

I=eye(d);
Nreps=1e5;
mu=0;
for ii=1:Nreps
    pind=randperm(d);
    PI=I(:,pind);
    mu=mu+PI*L*PI';
end
mu/Nreps;

(X*beta)'*L*(X*beta)
