clear
param.useGPU=false
param.NL=8;
param.NR=8;
param.freq_true=3; % freq used in regression model
Nmeastot=param.NL+param.NR;

% options for patternsearch
swarm_opts=optimoptions("patternsearch",'Algorithm',"classic",... 
            'Display','Iter', ...
            'MeshTolerance',1e-3,'MaxIterations',1000,'UseParallel',true);

%'MeshContractionFactor',0.9,'MeshExpansionFactor',4.0, ...

Aineq=eye(Nmeastot-1,Nmeastot); % inequality constraints
for ii=1:Nmeastot-1
    Aineq(ii,ii+1)=-1;
end
bineq=ones(Nmeastot-1,1);

[t_unif,~]=getSamplingSchedules(param.NL,param.NR,0,0.5); % initial guess for sampling

% function to be optimized
costfun=@(t) get_Aoptimal_cost(t,param);

[topt,fval]=patternsearch(costfun,t_unif,Aineq,bineq,[],[],zeros(Nmeastot,1),ones(Nmeastot,1),[],swarm_opts);


function cost = get_Aoptimal_cost(t,param)
X=constructX(t,param);
cost = trace(inv(X'*X));
end