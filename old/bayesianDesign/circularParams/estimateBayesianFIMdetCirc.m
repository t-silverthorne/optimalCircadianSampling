function [unif_vals,sd_unif,nu_vals,sd_nu] = estimateBayesianFIMdetCirc(NL,NR,tauA,tauB,cutoff_option,nreps)
%ESTIMATEBAYESIANFIMDETCIRC Summary of this function goes here
%   Detailed explanation goes here
sdev_cut=1e-1;
if nargin < 5
    cutoff_option='nreps'; % options: nreps, sdev
end
if nargin < 6
    nreps=30000;
end
nreps_loc=500;

% construct sample grids and information matrix
[mt_unif,mt_nu]=getSamplingSchedules(NL,NR,tauA,tauB);
Ntimes=NL+NR;
M=getBayesianFIMcirc(Ntimes,2);

switch cutoff_option
    case 'nreps'
        [unif_vals,nu_vals] = getSimulatedDataCirc(M,nreps,mt_unif,mt_nu);
        sd_unif=std(unif_vals)/sqrt(numel(unif_vals));
        sd_nu=std(nu_vals)/sqrt(numel(nu_vals));
    case 'sdev'
        sdev_now = Inf;
        unif_vals=[];
        while sdev_now>sdev_cut
            %fprintf('%d\n',sdev_now)
            [valsloc,~] = getSimulatedDataCirc(M,nreps_loc,mt_unif,mt_nu,true,false);
            unif_vals = [unif_vals valsloc];
            sdev_now = std(unif_vals)/sqrt(numel(unif_vals));
        end
        sd_unif=sdev_now;
        sdev_now = Inf;
        nu_vals=[];

        while sdev_now>sdev_cut
            %fprintf('%d\n',sdev_now)
            [~,valsloc] = getSimulatedDataCirc(M,nreps_loc,mt_unif,mt_nu,false,true);
            nu_vals=[nu_vals valsloc];
            sdev_now = std(nu_vals)/sqrt(numel(nu_vals));
        end
        sd_nu=sdev_now;
end
end