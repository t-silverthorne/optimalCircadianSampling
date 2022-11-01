function Sampvec = samplePosteriorMCMC(Nsamp,fnames,t_obs_MAT,Y_obs_MAT,model,method,settings)
% returns: Sampvec, a matrix of parameters sampled from the posterior 

% settings for MCMC
T=300; % length of Markov Chain


switch model
    case 'cosinorOneFreq'
        sd_vec=[.1 .1 .1];
    case 'cosinorTwoFreq'
        sd_vec=[5 1 1 5 1 1];
    case 'gaussian'
        sd_vec=[1];
end

% pre load gaussian noise, set up GPU arrays
switch model
    case 'cosinorOneFreq'
        theta_len=length(fnames);
        Sampvec=NaN(Nsamp,theta_len);
        Y=randn([Nsamp,theta_len,T]);
        tvec=reshape(t_obs_MAT,1,numel(t_obs_MAT));
        yvec=reshape(Y_obs_MAT,1,numel(Y_obs_MAT));
        if settings.run_gpu
            Sampvec=gpuArray(Sampvec);
            Y=gpuArray(Y);
            tvec=gpuArray(tvec);
            yvec=gpuArray(yvec);
        end
    case 'gaussian'
        Sampvec=NaN(Nsamp,1); 
        Y=randn([Nsamp,1,T]);
        yvec=reshape(Y_obs_MAT,1,numel(Y_obs_MAT));
        if settings.run_gpu
            Sampvec=gpuArray(Sampvec);
            Y=gpuArray(Y);
            yvec=gpuArray(yvec);
        end
end

% define target distribution
switch model
    case 'cosinorOneFreq'
        mutheta1=settings.acro_est;
        sigtheta1=settings.sig;
        muamp1=settings.amp_est;
        sigamp1=settings.sig;
        muf1=settings.freq_est;
        sigf1=settings.sig;
        logp=@(pmat) - (sum((yvec-abs(pmat(:,1)).*cos(2*pi*tvec.*( abs(pmat(:,3)) ) - ...
            ( pi+ pi*tanh(pmat(:,2)) ) )).^2/2,2) + ...
             (muamp1-pmat(:,1)).^2/2/sigamp1^2 + ...
                (mutheta1-pmat(:,2)).^2/2/sigtheta1^2 + ...
                  (muf1-pmat(:,3)).^2/2/sigf1^2) ;

    case 'gaussian'
        s0=settings.sigma0;
        sigma=settings.sigma;
        mu0=settings.mu0;
        n=length(yvec);

        logp=@(mu) - sum((yvec-mu).^2/2/sigma^2,2) - (mu-mu0).^2/2/s0^2;
end

switch model
    case 'cosinorOneFreq'
    if ~isfield(settings,'proposal_method')
        disp('proposal method not specified, using default')
        settings.proposal_method='iterative';
    end
    proposal_method=settings.proposal_method;
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
    case 'gaussian'
        Sampvecloc=mu0+randn([Nsamp,1])*s0;
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

