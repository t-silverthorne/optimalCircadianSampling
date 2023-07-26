close all
clear
rng(111)
%% simulate measured spectrum
% setup
f0=4;        % frequency in data
Amp=2.3;     % amp
sig=0.1;     % noise strength
t=(0:.2:1)'; % measurement times
t=t(1:end-1);
omegas=2*pi*(0:.1:7);  % frequency grid
omegas=omegas(2:end-1); % useful to avoid end frequencies

% data
y=Amp*cos(2*pi*f0*t-0.2)+.05*randn(length(t),1);

% estimate LSP
Pomegas=custom_LSP(t,y,omegas);


%% set coarseness of grid for next round of optimization
tau=(0:.005:1)'; % candidate measurement times
tau=tau(1:end-1);
mu0= ones(length(t),1); % existing times not being varied
mu = ones(length(tau),1)/length(tau); % weights for candidate meas times
t3d=reshape(t,[1,1,length(t)]);


%% construct quadprog
% constant term (useful for comparison but not in optimization)
a=0;
for ii=1:length(omegas)
    for jj=1:length(omegas)
        a=a+sum(Pomegas(ii)*Pomegas(jj)*exp(1j*(omegas(ii)-omegas(jj))*(t-t')),'all');
    end
end
a=real(a);

% linear term
f=arrayfun(@(ind) sum((Pomegas*Pomegas').*exp(1j*(omegas-omegas').*(tau(ind)-t3d)),'all'),(1:length(mu))');
f=f+conj(f);

% quadratic term
Htau=zeros(length(tau),length(tau));
for ii=1:length(omegas)
    for jj=1:length(omegas)
        Htau=Htau+Pomegas(ii)*Pomegas(jj)*exp(1j*(omegas(ii)-omegas(jj))*(tau-tau'));
    end
end
Htau=real(Htau);

%% run quadprog
clf
options = optimoptions('quadprog','Display','none');


tic; [muopt,copt]=quadprog(2*Htau,f,[],[],[],[],zeros(length(mu),1),ones(length(mu),1),mu,options);
time_quadprog=toc;

plot(tau,muopt/sum(muopt),'.k')
hold on

%% try simulated annealing

maxit=1e4;

opts = optimoptions(@simulannealbnd,'Display','none', ...
            'MaxIterations',maxit,'DisplayInterval',1e3,'ReannealInterval',50);

copt

tic;
muopt_sa = simulannealbnd(@(mu) mu'*Htau*mu+f'*mu, ...
           mu,zeros(length(tau),1),ones(length(tau),1),opts);
time_sa_cont=toc;

plot(tau,muopt_sa/sum(muopt_sa),'.b')

%% try discrete simulated annealing
Nmeas=5;
fdisc = @(mu) mu'*Htau*mu+f'*mu;
X=zeros(length(tau),1); % initial guess
X(randsample(1:length(tau),Nmeas,false))=1;

fX=fdisc(X); % initial f value
t=1; % counter
maxiter=1e4; % max numbe rof iterations
Tnow=20; % intial temperature
beta=.99999; % cooling rate (for annealing schedule)
tic
while t<maxiter
    ind1=randsample(find(X~=0),1); % indices to swap
    ind2=randsample(find(X==0),1);

    Xcand=X; % make swap
    Xcand(ind1)=0; 
    Xcand(ind2)=1;
    
    fcand=fdisc(Xcand);
    alpha=min( exp( - (fcand - fX )/Tnow) , 1);
    if rand<=alpha
        X=Xcand;
        fX=fcand;
    end
    Tnow=beta*Tnow;
    %fprintf('%f \n',fX)
    t=t+1;
end
time_sa_discrete=toc;
xline(tau(X==1),'-r')
fprintf('Quadprog:       %f\n',time_quadprog)
fprintf('SA continuous:  %f\n',time_sa_cont)
fprintf('SA discrete:    %f\n',time_sa_discrete)
ylim([0 0.05])
%fprintf('Temperature: %f\n', Tnow)
%fprintf('Solution: %f \n',fX)

% %%
% % compare against naive implementations
% tic
% a+mu'*Htau*mu+f'*mu
% toc
% tic
% Jspec_wL2([t;tau],[mu0;mu],omegas,Pomegas)
% toc
% tic
% Jspec_wL2_slow([t;tau],[mu0;mu],omegas,Pomegas)
% toc
%  

function J = Jspec_wL2(t,mu,omegas,Pomegas)
% weighted L2 cost function for spectral optimization
t=reshape(t,1,1,length(t));   
mu=reshape(mu,1,1,length(mu));
J1=exp(1j*(omegas-omegas').*t);
J1=sum(mu.*J1,3);
J1=J1.*conj(J1);
J2=Pomegas*Pomegas';
J=sum(J1.*J2,'all');
end

function J = Jspec_wL2_slow(t,mu,omegas,Pomegas)
J=0;
t=reshape(t,length(t),1);
mu=reshape(mu,length(mu),1);
for p=1:length(omegas)
    for k=1:length(omegas)
        J=J+Pomegas(p)*Pomegas(k)*abs(exp(1j*(omegas(p)-omegas(k))*t)'*mu)^2;
    end
end
end

