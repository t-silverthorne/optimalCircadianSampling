clear
checkUseGPU
tic
param.NL=4;
param.NR=4;
param.freq_true=1.4; % freq used in regression model
param.Amp=0.05;
param.acro=2;
beta=[0; sin(param.acro); cos(param.acro)];
[t,~]=getSamplingSchedules(param.NL,param.NR,0,0.5);
X=constructX(t,param);
d=length(t);
L=X*((X'*X)\X');
I=eye(d);

Nreps=1e5;
sigma=3;

% estimate first expectation value
T=@(Y) (X*((X'*X)\(X'*Y)) -Y)'*(X*((X'*X)\(X'*Y)) -Y);

Tvals=NaN(1,Nreps);
for ii=1:Nreps
    Y=param.Amp*cos(2*pi*param.freq_true*t'-param.acro)+sigma*randn(d,1);
    Tvals(ii)=T(Y);
end
mean(Tvals)
sigma^2*trace(I-L )

% estimate second
T=@(Y,PI) -norm((L-eye(d))*PI*Y,2).^2;

Tvals=NaN(1,Nreps);
for ii=1:Nreps
    Y=param.Amp*cos(2*pi*param.freq_true*t'-param.acro)+sigma*randn(d,1);
    pind=randperm(d);
    Tvals(ii)=T(Y,I(:,pind));
end
mean(Tvals)
-(sigma^2*trace(I-L) + (X*beta)'*(X*beta)) +   svd(L)'*svd(X*beta*(X*beta)')


