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
    beta0=rand;
    beta1=Amp*cos(phi);
    beta2=Amp*sin(phi);
    
    fit=fit_cosinor_model(beta0+ beta1*cos(2*pi*freq*mt)+beta2*sin(2*pi*freq*mt) +randn(Nbatch,N),mt,1/freq);
    
    beta_vec = [beta0;beta1;beta2];
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
for ii=1:1e3
    N=randsample(5:10,1); % random number of measurements
    mt=rand(1,N); % random measurement schedule
    Amp=rand; % amplitude
    freq=10*rand; % freq
    phi=2*pi*rand; % acro
    beta0=rand;
    beta1=Amp*cos(phi);
    beta2=Amp*sin(phi);
    
    Nbatch=1e3;
    fit   = fit_cosinor_model(beta0+ beta1*cos(2*pi*freq*mt)+beta2*sin(2*pi*freq*mt) +randn(Nbatch,N),mt,1/freq);
    Xr    = constructReducedX(mt,freq);
    est=mean(abs(exp(1j*fit.acrophases_rad) - exp(1j*phi) ).^2) ;
    exact=sqrt(1+2*abs(beta1+beta2)*sqrt(3/N)/(abs(beta1)+abs(beta2))  + 3 * norm([beta1 beta2],2)^2/norm([beta1 beta2],1)^2 )*...
        sqrt(1+2*(beta1+beta2)/Amp/sqrt(N) + norm([beta1 beta2],2)^2/Amp^2+ trace(inv(Xr'*Xr))/Amp^2) - (beta1+beta2)/sqrt(N)/Amp + abs(beta1+beta2)*sqrt(3/N)/norm([beta1 beta2],1);
    if est/exact>0.5
        disp("HOORAY")
    end
end
%(1-2*Delta)*mean(var(fit.betas(1,:))+var(fit.betas(2,:)))/Amp^2
%trace(inv(Xr'*Xr))
%mean(var(fit.betas(1,:))+var(fit.betas(2,:)))