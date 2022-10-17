function [unif_vals,sd_unif,nu_vals,sd_nu] = estimateBayesianFIMdet(NL,NR,tau,cutoff_option,nreps)
sdev_cut=1e-1;
if nargin < 4
    cutoff_option='nreps'; % options: nreps, sdev
end
if nargin < 5
    nreps=30000;
end
nreps_loc=500;

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
        sd_unif=std(unif_vals)/sqrt(numel(unif_vals));
        sd_nu=std(nu_vals)/sqrt(numel(nu_vals));
    case 'sdev'
        sdev_now = Inf;
        unif_vals=[];
        while sdev_now>sdev_cut
            %fprintf('%d\n',sdev_now)
            [valsloc,~] = getSimulatedData(M,nreps_loc,mt_unif,mt_nu,true,false);
            unif_vals = [unif_vals valsloc];
            sdev_now = std(unif_vals)/sqrt(numel(unif_vals));
        end
        sd_unif=sdev_now;
        sdev_now = Inf;
        nu_vals=[];

        while sdev_now>sdev_cut
            %fprintf('%d\n',sdev_now)
            [~,valsloc] = getSimulatedData(M,nreps_loc,mt_unif,mt_nu,false,true);
            nu_vals=[nu_vals valsloc];
            sdev_now = std(nu_vals)/sqrt(numel(nu_vals));
        end
        sd_nu=sdev_now;
end
    %case 'sdev'
end

