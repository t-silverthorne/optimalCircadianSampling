function run_patt_search_batch(Numfreq,freqmin,freqmax,Numamp,ampmin,ampmax)
simtype='fast';          % how many Monte Carlo samples to include
solver='pattern_search'; % options simulanneal or pswarm
addpath('utils/')
checkUseGPU              % uses simType to construct param

param.NL=4;              % samples in left interval
param.NR=4;              % samples in right interval
Nmeastot=param.NL+param.NR;
param.noise=1;         % the noise actually used in the simulation

% options for patternsearch
swarm_opts=optimoptions("patternsearch",'Algorithm',"nups",... 
            'Display','None', ...
            'MeshTolerance',1e-3,'MaxIterations',1000,'UseParallel',true);

% inequality constraints
Aineq=eye(Nmeastot-1,Nmeastot);
for ii=1:Nmeastot-1
    Aineq(ii,ii+1)=-1;
end
bineq=ones(Nmeastot-1,1);

[t_unif,~]=getSamplingSchedules(param.NL,param.NR,0,0.5); % initial guess for sampling

freqvals=logspace(log10(freqmin),log10(freqmax),Numfreq);
ampvals =logspace(log10(ampmin),log10(ampmax),Numamp);
nu_mat=NaN(Numfreq,Numamp);
unif_mat=NaN(Numfreq,Numamp);
for ii=1:Numfreq
    fprintf('Running %d\n',ii)
    for jj=1:Numamp
        
        param.freq_true=freqvals(ii);     % freq used in regression model
        param.Amp=ampvals(jj);             % signal amplitude
    
        % function to be optimized
        costfun=@(t) -min(wrapsimulatePWR(param,t));
           
        [topt,~]=patternsearch(costfun,t_unif,Aineq,bineq,[],[],zeros(Nmeastot,1),ones(Nmeastot,1),[],swarm_opts);
           
        simtype='medium';
        checkUseGPU
        [~,p_unif]=simulatePWR(param,t_unif);
        [~,p_nu]=simulatePWR(param,topt);
        
        nu_mat(ii,jj)  = min(p_nu);
        unif_mat(ii,jj)= min(p_unif);
   end
end
filename="results/"+ "Numfreq_" +  string(Numfreq) + ...
                  "_min_" + string(freqmin) + ...
                  "_max_" + string(freqmax) + ...
                  "_Numamp_" +string(Numamp) + ...
                  "_min_"   +string(ampmin) + ...
                  "_max_"   +string(ampmax) + '.mat';
save(filename);%,'nu_mat','unif_mat')
end