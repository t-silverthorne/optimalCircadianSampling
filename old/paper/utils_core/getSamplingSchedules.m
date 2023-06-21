function [mt_unif,mt_nu] = getSamplingSchedules(NL,NR,tauA,tauB)
%%%%%%%%%%%%%%
% GETSAMPLINGSCHEDULES returns uniform and non-uniform sampling schedules.
% INPUT: 
%   NL      number of measurements in left interval
%   NR      number of measurements in right interval
%   tauA    start of left interval
%   tauB    start of right interval
% OUTPUT:
%   mt_unif uniform measurement grid with NL+NR points
%   mt_nu   non_uniform measurement grid
%%%%%%%%%%%%%%
Ntimes=NL+NR;
mt1=linspace(0,tauB-tauA,NL+1);
mt1=mt1(1:end-1);
mt2=linspace(tauB-tauA,1,NR+1);
mt2=mt2(1:end-1);
mt_nu=[mt1 mt2]; % construct non-uniform grid 
mt_nu=mod(mt_nu+tauA,1);
mt_unif=linspace(0,1,Ntimes+1); % construct uniform grid 
mt_unif=mt_unif(1:end-1);
end

