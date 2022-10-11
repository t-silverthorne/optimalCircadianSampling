function [unif_vals,nu_vals] = estimateBayesianFIMdet(NL,NR,tau,cutoff_option,nreps)
if nargin < 4
    cutoff_option='nreps'; % options: nreps, sdev
end
if nargin < 5
    nreps=30000;
end
nreps_loc=3000;

% construct sample grids and information matrix
Ntimes=NL+NR;
mt1=linspace(0,tau,NL+1);
mt1=mt1(1:end-1);
mt2=linspace(tau,1,NR+1);
mt2=mt2(1:end-1);
mt_nu=[mt1 mt2]; % construct non-uniform grid 
mt_unif=linspace(0,1,Ntimes+1); % construct uniform grid 
mt_unif=mt_unif(1:end-1);
M=getBayesianFIM(Ntimes);

switch cutoff_option
    case 'nreps'
        [unif_vals,nu_vals] = getSimulatedData(M,nreps,mt_unif,mt_nu);
    case 'sdev'
        sdev_now = Inf;
        unif_vals=[];
        while sdev_now>1e-2
            fprintf('%d\n',sdev_now)
            [valsloc,~] = getSimulatedData(M,nreps_loc,mt_unif,mt_nu,true,false);
            unif_vals = [unif_vals valsloc];
            sdev_now = std(unif_vals)/mean(unif_vals);
        end
        sdev_now = Inf;
        nu_vals=[];

        while sdev_now>1e-2
            [~,valsloc] = getSimulatedData(M,nreps_loc,mt_unif,mt_nu,false,true);
            nu_vals=[nu_vals valsloc]
            sdev_now = std(nu_vals)/mean(nu_vals);
        end
end
    %case 'sdev'
end

