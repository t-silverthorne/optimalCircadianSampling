%% Amplitude
addpath('../utils')
nreps=1e1;  % how many estimates to run
Nbatch=1e4; % how many samples to use per estimate
ubvec=NaN(1,nreps); % check upper bound on A
lbvec=NaN(1,nreps); % check lower bound on A
varvec=NaN(1,nreps); % check bound on variance
% verify upper bound for amplitude bias
for ii=1:nreps % chec
    N=randsample(4:25,1); % random number of measurements
    mt=rand(1,N); % random measurement schedule
    
    Amp=10*rand; % amplitude
    freq=10*rand; % freq
    phi=2*pi*rand; % acro
    beta1=rand;
    beta2=Amp*cos(phi);
    beta3=Amp*sin(phi);
    
    fit=fit_cosinor_model(beta1+ beta2*cos(2*pi*freq*mt)+beta3*sin(2*pi*freq*mt) +randn(Nbatch,N),mt,1/freq);
    
    beta_vec = [beta1;beta2;beta3];
    p.freq=freq;
    X = constructX(mt,p);
    L=(X'*X)\X';
    P=eye(3);
    P(1,1)=0;
    
    ubvec(ii) = mean(fit.amplitudes)/sqrt(beta_vec'*P*beta_vec + trace(L'*P*L)) <1;
    if ~ubvec(ii)
        "bound violated"
        mean(fit.amplitudes)
        sqrt(beta_vec'*P*beta_vec + trace(L'*P*L))
    end
    lbvec(ii) = norm(P*beta_vec,1)/sqrt(3)/mean(fit.amplitudes) <1;
    varvec(ii) = var(fit.amplitudes)/(beta_vec'*P*beta_vec + trace(L'*P*L) - norm(P*beta_vec,1)/3 ) < 1;
end
fprintf("\nUpper amplitude bound:    %d\n",prod(ubvec))
fprintf("lower amplitude bound:    %d\n", prod(lbvec))
fprintf("variance amplitude bound: %d\n", prod(varvec))

% verify lower bound
%% Acrophase

N=randsample(10:20,1); % random number of measurements
mt=rand(1,N); % random measurement schedule

Amp=randsample(5:15,1)+rand; % amplitude
freq=1+rand; % freq
phi=2*pi*rand; % acro
beta1=rand;
beta2=Amp*cos(phi);
beta3=Amp*sin(phi);

Nbatch=1e5;
Nbvals =logspace(3,5,3);

for ii=1:length(Nbvals)
    Nbatch=Nbvals(ii);
    fit   = fit_cosinor_model(beta1+ beta2*cos(2*pi*freq*mt)+beta3*sin(2*pi*freq*mt) +randn(Nbatch,N),mt,1/freq);
    Xr    = constructReducedX(mt,freq);
    Delta = mean(fit.amplitudes-Amp)/Amp;
    est=mean(abs(exp(1j*fit.acrophases_rad) - exp(1j*phi) ).^2) ;
    exact=(var(fit.betas(1,:))+var(fit.betas(2,:)))/Amp^2;
    (est-exact)/exact
end
%(1-2*Delta)*mean(var(fit.betas(1,:))+var(fit.betas(2,:)))/Amp^2
%trace(inv(Xr'*Xr))
%mean(var(fit.betas(1,:))+var(fit.betas(2,:)))