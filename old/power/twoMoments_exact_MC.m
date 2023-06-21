
clear
addpath('utils/')
simtype='fast';
checkUseGPU
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

T=@(Y) -norm((L-eye(d))*Y,2)^2;
Tp=@(Y,PI) -norm((L-eye(d))*PI*Y,2)^2;

%% Term 1
% term 1 Monte Carlo
Tvals=NaN(1,Nreps);
for ii=1:Nreps
    Y=param.Amp*cos(2*pi*param.freq_true*t'-param.acro)+sigma*randn(d,1);
    Tvals(ii)=T(Y);
end
T1_MC=mean(Tvals);

% term 1 analytical
T1_exact=sigma^2*trace(L-I );


%% Term 2
% term 2 Monte Carlo
Tvals=NaN(1,Nreps);
for ii=1:Nreps
    Y=param.Amp*cos(2*pi*param.freq_true*t'-param.acro)+sigma*randn(d,1);
    pind=randperm(d);
    Tvals(ii)=Tp(Y,I(:,pind));
end
T2=mean(Tvals);

Lhat=NaN(d,d);
for ii=1:d
    for jj=1:d
        if ii==jj
            Lhat(ii,jj)=trace(L)/d;
        else
            Lhat(ii,jj)=sum(sum((1-eye(d)).*L))/d/(d-1);
        end
    end
end
T2_exact=sigma^2*trace(L-I )+ (X*beta)'*(Lhat-I)*X*beta;


%% Term 4 (VarP) 
U=L-I;
T3_exact= sqrt( (norm(diag(U),2)^2*3*sigma^4 + ...
            trace( ((1-I).*U)*U')*sigma^4 + ...
            trace(U*U')*sigma^4 + ...
            sum(sum(((1-I).*(diag(U)*diag(U)'))))*sigma^4) - ...
            T1_exact^2); % actually non-negative















