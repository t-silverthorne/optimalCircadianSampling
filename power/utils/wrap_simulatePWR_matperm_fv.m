function [acrovec,pwr_master] = wrap_simulatePWR_matperm_fv(param,nodes)
GPUsize=1e8;

% useful to construt before forming batches
Nmeas=param.Nmeas;
pall=perms(1:Nmeas);
Pmat=NaN(Nmeas,Nmeas,factorial(Nmeas));
I=eye(Nmeas);
for ii=1:factorial(Nmeas)
    Pmat(:,:,ii)= I(:,pall(ii,:)); % might be possible to vectorize
end
param.Pmat=Pmat;

% partition into batches and run
if param.Nmeas*param.Nperm*param.Nresidual*param.Nacro > GPUsize
    fprintf('Cutting into batches\n')
    Nres_init       = param.Nresidual;
    param.Nresidual = floor(GPUsize/param.Nmeas/param.Nacro/param.Nperm);
    pwr_master      = zeros(1,param.Nacro);
    nbatch          = 0;
    while nbatch*param.Nresidual < Nres_init
        fprintf(' batch %d\n',nbatch)
        [acrovec,pwr]=simulatePWR_matperm_fv(param,nodes);
        pwr_master=pwr_master+pwr;
        nbatch=nbatch+1;
    end
    pwr_master=pwr_master/nbatch;
else
    pwr_master      = simulatePWR_matperm_fv(param,nodes)
end
end

