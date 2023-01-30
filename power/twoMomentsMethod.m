clear
checkUseGPU
tic
param.NL=5;
param.NR=3;
param.freq_true=1.4; % freq used in regression model
param.Amp=2.1;
param.acro=2;
beta=[0; param.Amp*sin(param.acro); param.Amp*cos(param.acro)];
[~,t]=getSamplingSchedules(param.NL,param.NR,0,0.25);
X=constructX(t,param);
d=length(t);
L=X*((X'*X)\X');
I=eye(d);

Nreps=1e4;
sigma=2;

% estimate first expectation value
T=@(Y) -norm((L-eye(d))*Y,2)^2;

Tvals=NaN(1,Nreps);
for ii=1:Nreps
    Y=param.Amp*cos(2*pi*param.freq_true*t'-param.acro)+sigma*randn(d,1);
    Tvals(ii)=T(Y);
end
T1=mean(Tvals);
%sigma^2*trace(L-I )

% estimate second
T=@(Y,PI) -norm((L-eye(d))*PI*Y,2)^2;

Tvals=NaN(1,Nreps);
for ii=1:Nreps
    Y=param.Amp*cos(2*pi*param.freq_true*t'-param.acro)+sigma*randn(d,1);
    pind=randperm(d);
    Tvals(ii)=T(Y,I(:,pind));
end
T2=mean(Tvals);


% estimate the term you want to control
tct=NaN(1,Nreps);
for ii=1:Nreps
    Y=param.Amp*cos(2*pi*param.freq_true*t'-param.acro)+sigma*randn(d,1);
    pind=randperm(d);
    PI=I(:,pind);
    tct(ii)=trace(L*PI)+trace((X*beta)*(X*beta)'*PI');
end
T2-T1 - (-(X*beta)'*(X*beta)+mean(tct))

%% calculate terms individually

t1=NaN(1,Nreps);
t2=NaN(1,Nreps);
t3=NaN(1,Nreps);
t4=NaN(1,Nreps);
t4ap=NaN(1,Nreps);
t4ap2=NaN(1,Nreps);

for ii=1:Nreps
    Y=param.Amp*cos(2*pi*param.freq_true*t'-param.acro)+sigma*randn(d,1);
    pind=randperm(d);
    PI=I(:,pind);
    eps=randn(d,1);
    t1(ii)=sigma^2*eps'*PI'*(I-L)*PI*eps;
    t2(ii)=(X*beta)'*PI'*(I-L)*PI*eps;
    t3(ii)=eps'*PI'*(I-L)*PI*(X*beta);
    t4(ii)=(X*beta)'*PI'*(I-L)*PI*(X*beta);
    t4ap(ii)=-(X*beta)'*(X*beta)+(X*beta)'*PI'*L*PI*(X*beta);
    t4ap2(ii)=-(X*beta)'*(X*beta)+abs(trace(PI'*L*PI)*trace((X*beta)*(X*beta)'));
end
t1=mean(t1);t2=mean(t2);t3=mean(t3);t4=mean(t4);t4ap=mean(t4ap);
t4ap2=mean(t4ap2);
%-(t1+t2+t3+t4)
-(t1+t4)
-t1+t4ap
-t1+t4ap2
T2


%%
Nreps=1e5
mu=0;
for ii=1:Nreps
    pind=randperm(d);
    PI=I(:,pind);
    mu=mu+PI*L*PI';
end
mu/Nreps