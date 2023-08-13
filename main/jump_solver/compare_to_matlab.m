addpath('../utils/')
%%
c1rgb=[1.0, 0.01, 0.2]; % sparse color
c4rgb=[.6 .4 .8];       % dense color

Nmeas    = 8;
Amp      = 2.23;
freq     = 2.8;
[mt,~]   = getSamplingSchedules(Nmeas,0,0,0);


% Method 1: directly optimize power 
f_unif  = getMinPower(Amp,freq,mt);
tic;
[x,fval1]=fmincon(@(t) -getMinPower(Amp,freq,t),mt,[],[],[],[],zeros(1,Nmeas),ones(1,Nmeas));
t1=toc;
fval1


% Method 2: gradient based for lambda (no user derivative)
tic;
[x,~]=fmincon(@(t) -getMinLambda(t,freq),mt,[],[],[],[],zeros(1,Nmeas),ones(1,Nmeas));
fval2=getMinPower(Amp,freq,x);
t2=toc;


% Method 3: Gradient based for lambda (with user derivative)
options=optimoptions('fmincon',...
  'CheckGradients',false,...
  'SpecifyObjectiveGradient',true,'FiniteDifferenceStepSize',1e-10,'FiniteDifferenceType','central');

tic;
[x,~]=fmincon(@(t) get_MINUS_MinLambdaAndDiff(t,freq),mt,[],[],[],[],zeros(1,Nmeas),ones(1,Nmeas),[],options);
fval3=getMinPower(Amp,freq,x);
t3=toc;

% Method 4: mulitstart fmincon
% problem = createOptimProblem('fmincon','objective',...
%     @(t) get_MINUS_MinLambdaAndDiff(t,freq),'x0',mt,'lb',zeros(1,Nmeas),'ub',ones(1,Nmeas),'options',options);
% ms = MultiStart;
% tic
% [x,~] = run(ms,problem,20);
% fval4 = getMinPower(Amp,freq,x);
% t4=toc;


fprintf('\nUniform :                     %2.3f\n',f_unif)
fprintf('fmincon naive:                  %2.3f   %2.4f\n'  ,-fval1,t1)
fprintf('fmincon lambda:                 %2.3f   %2.4f\n'  ,fval2, t2)
fprintf('fmincon diff lambda:            %2.3f   %2.4f\n'  ,fval3, t3)
fprintf('multistart fmincon diff lambda: %2.3f   %2.4f\n'  ,fval4, t4)
