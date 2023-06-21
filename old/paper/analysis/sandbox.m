%% general version of loewr bound for bias false
addpath('../utils_core')
addpath('../utils_cost_fun')
d = 2;
N = 1e4;

A = 100*randn(d,d);
B = 0.1*randn(d,d);
x = rand(d,1);
veps = 5*randn(d,N);

M = A*x + B*veps;

% want LHS >= RHS
LHS = mean(sqrt(sum(M.*M,1)));
RHS = norm(A*x,1)/sqrt(d);
RHS/LHS

%% other version


%% verify bias estimates
p.Nmeas     = 19;
p.Nresidual = 1e4;
p.freq      = 3;
p.Amp       = 2;
p.noise     = 1;
acro        = 2*pi*rand; % choose random acrophase
p.Nbatch    = 2;
[t,~]=getSamplingSchedules(p.Nmeas,0,0,0);


eps=randn(p.Nresidual,p.Nmeas);
Y=p.Amp*cos(2*pi*t*p.freq-acro)+p.noise*eps;
X=constructX(t,p);

betas_obs = (X'*X)\(X'*Y');%pagemldivide(X'*X,pagemtimes(X',pagetranspose(Y))); 
phi_est   = atan2(betas_obs(2,:),betas_obs(3,:));
amp_est   = sqrt(betas_obs(2,:).^2 + betas_obs(3,:).^2);

amp_true  = p.Amp;
beta      = [0; p.Amp*sin(acro); p.Amp*cos(acro)];

bias = mean(amp_est);
std2 = var(amp_est);

P=diag([0,1,1]);

L = (X'*X)\X';


% lower bound is acrophase dependent, upper bound isnt
LHS_bias = norm(P*beta,1)/sqrt(p.Nmeas);
bias;
RHS_bias = sqrt(beta'*P*beta + trace(L'*P*L));
fprintf('\n')
fprintf('Using N=%.0e samples\n',p.Nresidual)
fprintf('  Lower bound: %.4f\n',LHS_bias)
fprintf('  Bias:        %.4f\n',bias)
fprintf('  Upper bound: %.4f\n',RHS_bias)










