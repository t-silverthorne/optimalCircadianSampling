function [acrovec,pwr_master,est] = wrap_simulatePWR_matperm_fv(param,nodes)
if param.useParallel
    GPUsize=1e8/param.poolobj.NumWorkers;
else
    GPUsize=1e8;
end

% useful to construt perms before forming batches
if ~isfield(param,'Pmat')
    param.Pmat=construct_all_perms(param.Nmeas);
end

% partition into batches and run
if param.Nmeas*param.Nperm*param.Nresidual*param.Nacro > GPUsize
    %fprintf('Cutting into batches\n')
    Nres_init       = param.Nresidual;
    param.Nperm     = min(factorial(param.Nmeas),param.Nperm);
    param.Nresidual = floor(GPUsize/param.Nmeas/param.Nacro/param.Nperm);
    pwr_master      = zeros(1,param.Nacro);

    amp_sum_master     = zeros(1,param.Nacro);
    amp_sum_sq_master  = zeros(1,param.Nacro);
    phi_harm_sum       = zeros(1,param.Nacro);
    nbatch             = 0;
    while nbatch*param.Nresidual < Nres_init
        %fprintf(' batch %d\n',nbatch)
        [acrovec,pwr,est] = simulatePWR_matperm_fv(param,nodes);
        pwr_master        = pwr_master+pwr;
        amp_sum_master    = amp_sum_master + reshape(est.amp_sum,1,param.Nacro);
        amp_sum_sq_master = amp_sum_sq_master + reshape(est.amp_sq_sum,1,param.Nacro);
        phi_harm_sum      = phi_harm_sum + reshape(est.phi_harm_sum,1,param.Nacro);
        nbatch=nbatch+1;
    end


    pwr_master=pwr_master/nbatch;
else
    [acrovec,pwr,est] = simulatePWR_matperm_fv(param,nodes);
    pwr_master        = pwr;
    amp_sum_master    = reshape(est.amp_sum,1,param.Nacro);
    amp_sum_sq_master = reshape(est.amp_sq_sum,1,param.Nacro);
    phi_harm_sum      = reshape(est.phi_harm_sum,1,param.Nacro);
    nbatch=1;
end
clear est
est.amp_mu   = amp_sum_master/param.Nresidual/nbatch;   
amp_m2       = amp_sum_sq_master/param.Nresidual/nbatch;
est.amp_st   = sqrt(amp_m2 - est.amp_mu.^2);
est.phi_mu   = angle(phi_harm_sum/param.Nresidual/nbatch);
est.phi_cvar = abs(phi_harm_sum/param.Nresidual/nbatch);
end

