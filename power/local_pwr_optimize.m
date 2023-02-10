% optimize power by finding best place to move measurement to at each
% iteration, since only one node moves, the nodes can be checked
% independently

close all
clear
simtype='medium';          % how many Monte Carlo samples to include
addpath('utils/')
checkUseGPU              % uses simType to construct param
param.NL=4;              % samples in left interval
param.NR=4;              % samples in right interval
Nmeastot=param.NL+param.NR;
param.freq_true=1;     % freq used in regression model
param.Amp=1;             % signal amplitude
param.noise=.5;          % the noise actually used in the simulation

param.Nacro=16;

[t_unif,~]=getSamplingSchedules(param.NL,param.NR,0,0.5); % initial guess for sampling

iter=1;
Niter=5;
t_now=t_unif;

Nunif_std_est=6*4;
minpwr_unif_vec=NaN(1,Nunif_std_est);
parfor ii=1:Nunif_std_est
    [~,pwr_unif]=simulatePWR(param,t_unif);
    minpwr_unif_vec(ii)=min(pwr_unif)
end
std(minpwr_unif_vec)
%%

fprintf('Uniform measurement, best power: %1.2f\n',mean(minpwr_unif_vec))
dt=.02;
%t_now=sort(rand(1,param.NL+param.NR));
while iter<Niter
    t_prop=NaN(numel(t_unif),numel(t_unif)); % list of proposed measurement schedules
    pwr_prop=NaN(numel(t_unif),1); % min powers corresponding to measurement schedules
    parfor ii=1:numel(t_unif)-1
        t_loc   = [t_now(1:ii-1) t_now(ii+1:end)]; % remove ii-th measurement time
        t_cands = t_now(ii):dt:t_now(ii+1);        % consider all places between ii-th and ii+1-th times
        t_cands=t_cands(1:end-1); % last point is irrelevant
        pwrloc  = NaN(1,numel(t_cands));           % each gives a minimum power
        for jj=1:numel(t_cands) % simulate each power
            tloc_cand= [t_cands(jj) t_loc];
            [~,pwr_vec]=simulatePWR(param,tloc_cand);
            pwrloc(jj)=min(pwr_vec);
        end
        [pwr_min_max,ind_max]=max(pwrloc); % find best choice
        pwr_prop(ii)=pwr_min_max; 
        t_prop(ii,:)=sort([t_cands(ind_max) t_loc]);
    end    
    [best_pwr,best_move_ind]=max(pwr_prop);
    t_now=t_prop(best_move_ind,:);
    fprintf('Iteration %d, best power: %1.2f\n',iter,best_pwr)
    iter=iter+1;
    clf
    xlim([0,1])
    xline(t_now)
    xlim([0,1])
    drawnow
end

%%
[~,pwr]=simulatePWR(param,t_now)
min(pwr)
% Ncand=2^5;
% t_cands=linspace(0,1,Ncand);

% old iteration, naive, step through all time points and consider all
% places in the interval where they could go
% while iter<Niter
%     t_prop=NaN(numel(t_unif),numel(t_unif));
%     pwr_prop=NaN(numel(t_unif),1);
%     for ii=1:numel(t_unif)
%         disp(ii)
%         t_loc=[t_now(1:ii-1) t_now(ii+1:end)];
%         pwrloc=NaN(1,numel(t_cands));
%         for jj=1:numel(t_cands)
%             tloc_cand= [t_cands(jj) t_loc];
%             [~,pwr_vec]=simulatePWR(param,tloc_cand);
%             pwrloc(jj)=min(pwr_vec);
%         end
%         [pwr_min_max,ind_max]=max(pwrloc);
%         pwr_prop(ii)=pwr_min_max;
%         t_prop(ii,:)=sort([t_cands(ind_max) t_loc]);
%     end    
%     [best_pwr,best_move_ind]=max(pwr_prop);
%     t_now=t_prop(best_move_ind,:);
%     fprintf('Iteration %d, best power: %1.2f\n',iter,best_pwr)
%     iter=iter+1;
%     clf
%     xlim([0,1])
%     xline(t_now)
%     xlim([0,1])
%     drawnow
% end
fprintf('complete\n')