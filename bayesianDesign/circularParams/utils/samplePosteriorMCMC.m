function Sampvec = samplePosteriorMCMC(Nsamp,fnames,t_obs_MAT,Y_obs_MAT,model,method,settings)
% inputs:
%   Nsamp:  desired number of samples
%   fnames: names of fields in parameter struct
%   t_obs_MAT: matrix with rows correspond to distinct measurements, columns
%   measurement times
%   Y_obs_MAT: simulated data corresponding to t_obs_MAT
%   model,method,settings: specify choice of posterior for model being used

% returns: Sampvec, a matrix of parameters sampled from the posterior 
proposal_method='iterative';
switch model
    case 'cosinorOneFreq'
        sd_vec=[5 1 1];
    case 'cosinorTwoFreq'
        sd_vec=[5 1 1 5 1 1];
end
T=100; % length of Markov Chain
theta_len=length(fnames);
Sampvec=NaN(Nsamp,theta_len);
% TODO: decide if log helps
logp=@(Y) evalLogImproperPosterior(t_obs_MAT, ...
                 Y_obs_MAT,...
                 getTheta(Y,fnames),model,method,settings);
if Nsamp>1
    numworkers=inf;
else
    numworkers=0;
end

switch proposal_method
    case 'fixed'
        Xtvec=samplePrior(Nsamp,model,method,settings);
    case 'iterative'
         if size(t_obs_MAT,1)>1
                Xtvec=samplePosteriorMCMC(Nsamp,fnames,t_obs_MAT(2:end,:),Y_obs_MAT(2:end,:),model,method,settings);
            else
                Xtvec=samplePrior(Nsamp,model,method,settings);
         end
end

parfor(N=1:Nsamp,numworkers)
%for N=1:Nsamp
    Xt=Xtvec(N,:);
    for t=1:T
        Y=randn(1,theta_len).*sd_vec+Xt;
        diffp=logp(Y)-logp(Xt);
        alpha=min( diffp+log(normpdfvec(Xt,Y))-log(normpdfvec(Y,Xt)),0);
        if log(rand)<alpha
            if evalPrior(getTheta(Y,fnames),model,method,settings)==0
                disp('oh no')
            end
            Xt=Y;
        end

    end
    Sampvec(N,:)=Xt;
end

       


end

