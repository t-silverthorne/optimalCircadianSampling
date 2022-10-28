function Sampvec = samplePosteriorMCMC(Nsamp,fnames,t_obs_MAT,Y_obs_MAT,model,method,settings)
% inputs:
%   Nsamp:  desired number of samples
%   fnames: names of fields in parameter struct
%   t_obs_MAT: matrix with rows correspond to distinct measurements, columns
%   measurement times
%   Y_obs_MAT: simulated data corresponding to t_obs_MAT
%   model,method,settings: specify choice of posterior for model being used

% returns: Sampvec, a matrix of parameters sampled from the posterior 

T=100; % length of Markov Chain
if ~isfield(settings,'parallel_mode')
    disp('parallel not specified, using default')
    settings.parallel_mode='vectorize';
end
if ~isfield(settings,'proposal_method')
    disp('proposal method not specified, using default')
    settings.proposal_method='iterative';
end
proposal_method=settings.proposal_method;

switch settings.parallel_mode
    case 'parfor'
        switch model
            case 'cosinorOneFreq'
                sd_vec=[5 1 1];
            case 'cosinorTwoFreq'
                sd_vec=[5 1 1 5 1 1];
        end
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
                %Xtvec=randn(Nsamp,theta_len);
            case 'iterative'
                 if size(t_obs_MAT,1)>1
                        Xtvec=samplePosteriorMCMC(Nsamp,fnames,t_obs_MAT(2:end,:),Y_obs_MAT(2:end,:),model,method,settings);
                    else
                        Xtvec=samplePrior(Nsamp,model,method,settings);
                        %Xtvec=randn(Nsamp,theta_len);
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
                    Xt=Y;
                end
        
            end
            Sampvec(N,:)=Xt;
        end


    case 'vectorize'
        switch model
            case 'cosinorOneFreq'
                sd_vec=[5 1 1];
                sdmat=repmat(sd_vec,Nsamp,1);
            case 'cosinorTwoFreq'
                sd_vec=[5 1 1 5 1 1];
                sdmat=repmat(sd_vec,Nsamp,1);
        end
        theta_len=length(fnames);
        Sampvec=NaN(Nsamp,theta_len,T);
        Y=randn([Nsamp,theta_len,T]);

        % TODO: decide if log helps
        logp=@(Y) transpose(arrayfun(@(ind) evalLogImproperPosterior(t_obs_MAT, ...
                         Y_obs_MAT,...
                         getTheta(Y(ind,:),fnames),model,method,settings),1:size(Y,1)));
        if Nsamp>1
            numworkers=inf;
        else
            numworkers=0;
        end
        
        switch proposal_method
             case 'fixed'
                Sampvecloc=samplePrior(Nsamp,model,method,settings);
                %Sampvecloc=randn(Nsamp,theta_len);
            case 'iterative'
                 if size(t_obs_MAT,1)>1
                        Sampvecloc=samplePosteriorMCMC(Nsamp,fnames,t_obs_MAT(2:end,:),Y_obs_MAT(2:end,:),model,method,settings);
                    else
                        Sampvecloc=samplePrior(Nsamp,model,method,settings);
                        %Sampvecloc=randn(Nsamp,theta_len);
                 end
        end
        
        Sampvec(:,:,1)=Sampvecloc;
        for t=2:T
            Z=Y(:,:,t).*sdmat+Sampvec(:,:,t-1);
            diffp=logp(Z)-logp(Sampvec(:,:,t-1));
            logdiff=log(prod(normpdf(Sampvec(:,:,t-1),Z,sdmat),2))- ...
                        log(prod(normpdf(Z,Sampvec(:,:,t-1),sdmat),2));
            alphamat=min( diffp+logdiff,0);
            rmat=log(rand(Nsamp,1));
            Sampvec(:,:,t)=(rmat<alphamat).*Y(:,:,t).*sdmat + Sampvec(:,:,t-1);
        end
        Sampvec=Sampvec(:,:,T);    


end

