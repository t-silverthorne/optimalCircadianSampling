function [acrovec,pwr_master,est] = simulatePWR_rank4(param,nodes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the power, bias, and variability of amplitude and
% acrophase. From the outputs, you can calculate the cost function used in
% multi-objective optimization
% OUTPUT:
% acrovec    =  vector of phases used in calucation
% pwr_master =  power as a function of acrophase
% est        =  struct which stores various statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjust compute size based on if parallel sessions are running
mem_size=1e7;
if size(gcp('nocreate'),1)>0
    mem_size=mem_size/(1+gcp('nocreate').NumWorkers);
end

% partition into batches
total_samples=param.Nmeas*param.Nperm*param.Nresidual*param.Nacro;

%fprintf('Cutting into batches\n')
Nres_init       = param.Nresidual;
param.Nperm     = min(factorial(param.Nmeas),param.Nperm);
param.Nresidual = min(Nres_init,floor(mem_size/param.Nmeas/param.Nacro/param.Nperm));
pwr_master      = zeros(1,param.Nacro);

amp_sum_master     = zeros(1,param.Nacro);
amp_sum_sq_master  = zeros(1,param.Nacro);
phi_harm_sum       = zeros(1,param.Nacro);
nbatch             = 0;

% run each batch
while nbatch*param.Nresidual < Nres_init
    NL=param.NL;
    NR=param.NR;
    Nresidual=param.Nresidual;
    Nacro=param.Nacro;
    Amp=param.Amp;
    Nmeas=NL+NR;
    
    Nperm=param.Nperm;
    freq_true=param.freq_true;
    
    if isnumeric(nodes) % construct time grid
        t=nodes;
        Nmeas=length(t);
    else
        switch nodes
            case 'uniform'
                [t,~]=getSamplingSchedules(NL,NR,0,0.5);
            case 'non-uniform'
                [~,t]=getSamplingSchedules(NL,NR,0,0.3);
        end
    
    end
    
    if size(t,1)>size(t,2)% make sure t is row vector
        t=t';
    end
    
    
    acrovec=linspace(0,2*pi,Nacro); % get acrophases
    acromat=reshape(acrovec,1,1,1,Nacro);
    
    if param.useGPU
        eps=randn(Nresidual,Nmeas,1,Nacro,'gpuArray');
    else
        eps=randn(Nresidual,Nmeas,1,Nacro);
    end
   
    Y=Amp*cos(2*pi*t*freq_true-acromat)+param.noise*eps;  % simulated data
    %Y=repmat(Y(1,:,1,1),size(Y,1),1,size(Y,3),size(Y,4));
    clear eps
   
    X=constructX(t,param); % construct linear model
    if ~param.useGPU
        X=gather(X);
    end
    
    betas_obs=pagemldivide(X'*X,pagemtimes(X',pagetranspose(Y))); % observed error
    
    % accuracy of amplitude and phase estimates
    est.amp_sum      = sum(sqrt(betas_obs(2,:,:,:).^2 + betas_obs(3,:,:,:).^2),2);
    est.amp_sq_sum   = sum(betas_obs(2,:,:,:).^2 + betas_obs(3,:,:,:).^2,2);
    phi              = atan2(betas_obs(2,:,:,:),betas_obs(3,:,:,:));
    est.phi_harm_sum = sum(exp(1j*phi),2);
    
    
    fits_obs=pagetranspose(pagemtimes(X,betas_obs));
    SSres_obs=sqrt(sum((fits_obs-Y).^2,2));
    clear betas_obs fits_obs
    
    %YI=NaN(size(Y,1),size(Y,2),Nperm,size(Y,4));
    %YI(:,:,1,:)=Y
    YI=get_permuted_Y_rank4(repmat(Y,1,1,Nperm,1),param); 
    
    betas=pagemldivide(X'*X,pagemtimes(X',pagetranspose(YI)));
    
    fits =pagetranspose(pagemtimes(X,betas));
    SSres=sqrt(sum((fits-YI).^2,2));
    pwr=sum(sum(SSres>SSres_obs,3)/Nperm>.95)/Nresidual;
    pwr=reshape(pwr,1,Nacro);        


    pwr_master        = pwr_master+pwr;
    amp_sum_master    = amp_sum_master + reshape(est.amp_sum,1,param.Nacro);
    amp_sum_sq_master = amp_sum_sq_master + reshape(est.amp_sq_sum,1,param.Nacro);
    phi_harm_sum      = phi_harm_sum + reshape(est.phi_harm_sum,1,param.Nacro);
    nbatch=nbatch+1;
end

pwr_master=pwr_master/nbatch;



clear est
est.amp_mu   = amp_sum_master/param.Nresidual/nbatch;   
amp_m2       = amp_sum_sq_master/param.Nresidual/nbatch;
est.amp_st   = sqrt(amp_m2 - est.amp_mu.^2);
est.phi_mu   = angle(phi_harm_sum/param.Nresidual/nbatch);
est.phi_cvar = abs(phi_harm_sum/param.Nresidual/nbatch);
end
