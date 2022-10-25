function Sampvec = samplePosteriorMCMC(Nsamp,fnames,t_obs_MAT,Y_obs_MAT,theta,model,method,settings)
proposal_method='fixed';
T=1000; % length of Markov Chain
theta_len=length(fieldnames(theta));
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
parfor(N=1:Nsamp,numworkers)
    switch proposal_method
        case 'fixed'
            Xt=samplePrior(1,model,method,settings);
        case 'iterative'
             if size(t_obs_MAT,1)>1
                    Xt=samplePosteriorMCMC(1,fnames,t_obs_MAT(2:end,:),Y_obs_MAT(2:end,:),theta,model,method,settings);
                else
                    Xt=samplePrior(1,model,method,settings);
             end
    end
    %Xt=rand(1,theta_len);

    for t=1:T
        Y=randn(1,theta_len)+Xt;
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

