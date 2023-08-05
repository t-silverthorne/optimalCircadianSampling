addpath('../utils/')
%%
Nmeas    = 12;
Amp      = 2.23;
freq     = 3.8;
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
  'SpecifyObjectiveGradient',true);

tic;
[x,~]=fmincon(@(t) get_MINUS_MinLambdaAndDiff(t,freq),mt,[],[],[],[],zeros(1,Nmeas),ones(1,Nmeas),[],options);
fval3=getMinPower(Amp,freq,x);
t3=toc;

%%Check that derivative matches Julia implementation
%[~,grad]=get_MINUS_MinLambdaAndDiff(mt,freq);
%-grad
 

fprintf('\nUniform :            %2.3f\n',f_unif)
fprintf('fmincon naive:       %2.3f   %2.4f\n'  ,-fval1,t1)
fprintf('fmincon lambda:      %2.3f   %2.4f\n'  ,fval2, t2)
fprintf('fmincon diff lambda: %2.3f   %2.4f\n'  ,fval3, t3)
