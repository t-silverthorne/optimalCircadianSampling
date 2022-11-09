function Sampvec = samplePosteriorMCMC(Nsamp,fnames,t_obs_MAT,Y_obs_MAT,model,method,settings)
speed=settings.speed; 
% inputs:
%   Nsamp:  desired number of samples
%   fnames: names of fields in parameter struct
%   t_obs_MAT: matrix with rows correspond to distinct measurements, columns
%   measurement times
%   Y_obs_MAT: simulated data corresponding to t_obs_MAT
%   model,method,settings: specify choice of posterior for model being used

% returns: Sampvec, a matrix of parameters sampled from the posterior 

T=1000; % length of Markov Chain
if ~isfield(settings,'parallel_mode')
    disp('parallel not specified, using default')
    settings.parallel_mode='vectorize';
end
if ~isfield(settings,'proposal_method')
    disp('proposal method not specified, using default')
    settings.proposal_method='iterative';
end
proposal_method=settings.proposal_method;
switch model
    case 'cosinorOneFreq'
        sd_vec=[.1 .1 .1]/5;
    case 'cosinorTwoFreq'
        sd_vec=[5 1 1 5 1 1];
end
switch settings.parallel_mode
    case 'parfor'
    
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
        theta_len=length(fnames);
        Sampvec=NaN(Nsamp,theta_len);
        Y=randn([Nsamp,theta_len,T]);
                
        if settings.run_gpu
            Sampvec=gpuArray(Sampvec);
            Y=gpuArray(Y);
        end


        % TODO: decide if log helps
        switch speed
            case 'slow'
                logp=@(Y) transpose(arrayfun(@(ind) evalLogImproperPosterior(t_obs_MAT, ...
                                 Y_obs_MAT,...
                                 getTheta(Y(ind,:),fnames),model,method,settings),1:size(Y,1)));
            case 'fast'
                tvec=reshape(t_obs_MAT,1,numel(t_obs_MAT));
                yvec=reshape(Y_obs_MAT,1,numel(Y_obs_MAT));
                
                mutheta1=settings.acro_est;
                sigtheta1=settings.sig;
                muamp1=settings.amp_est;
                sigamp1=settings.sig;
                muT1=settings.per_est;
                sigT1=settings.sig;

                if settings.run_gpu
                    tvec=gpuArray(tvec);
                    yvec=gpuArray(yvec);
                end
                logp=@(pmat) -sum((yvec-pmat(:,1).*cos(2*pi*tvec./pmat(:,3)-pmat(:,2))).^2/2,2) - ...
                             (muamp1-pmat(:,1)).^2/2/sigamp1^2 - ...
                                (mutheta1-pmat(:,2)).^2/2/sigtheta1^2 - ...
                                  (muT1-pmat(:,3)).^2/2/sigT1^2 ;
        end

        % pmat(ii,:) row of parameter matrix
        % yvec, flattened obserations, tvec flattened obs times
        
        if Nsamp>1
            numworkers=inf;
        else
            numworkers=0;
        end
        
        switch proposal_method
            case 'fixed'
                Sampvecloc=samplePrior(Nsamp,model,method,settings);
            case 'iterative'
                 if size(t_obs_MAT,1)>1
                        Sampvecloc=samplePosteriorMCMC(Nsamp,fnames,t_obs_MAT(2:end,:),Y_obs_MAT(2:end,:),model,method,settings);
                    else
                        Sampvecloc=samplePrior(Nsamp,model,method,settings);
                 end
        end

        Sampvec=Sampvecloc;
        if settings.run_gpu
            Sampvec=gpuArray(Sampvec); 
            sd_vec=gpuArray(sd_vec); % todo: check if this is necessary
        end
        for t=2:T
            Z=Y(:,:,t).*sd_vec+Sampvec;
            diffp=logp(Z)-logp(Sampvec);
            logdiff=log(prod(normpdf(Sampvec,Z,sd_vec),2)) - ...
                        log(prod(normpdf(Z,Sampvec,sd_vec),2));
            if settings.run_gpu
                alphamat=min( diffp+logdiff,gpuArray(0)); % might be overkill
                rmat=gpuArray(log(rand(Nsamp,1)));
            else
                alphamat=min( diffp+logdiff,0); % might be overkill
                rmat=log(rand(Nsamp,1));
            end
            Sampvec=(rmat<alphamat).*Y(:,:,t).*sd_vec + Sampvec;
        end

        Sampvec=gather(Sampvec);    
        

end

