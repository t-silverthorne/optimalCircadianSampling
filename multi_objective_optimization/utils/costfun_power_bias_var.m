function Jvec = costfun_power_bias_var(param,nodes)
% assuming we want to minimize J1 = bias + var - min power
% J2 is a weaker cost function, only containing terms that you can
% calculate without knowing the (in principal unknown) parameter values
switch param.method
    case 'matperm_fv'
        [acrovec,pwr,est]=wrap_simulatePWR_matperm_fv(param,nodes);
    case '4tensor'
        [acrovec,pwr,est]=simulatePWR_rank4(param,nodes);
end
 
Jvec(1)=max(abs(param.Amp-est.amp_mu));                 % amplitude error
Jvec(2)=max(min(mod([acrovec-est.phi_mu; est.phi_mu-acrovec],2*pi),[],1));   % phase error (angle of order parameter)
Jvec(3)=max(est.amp_st);                                % amp stdev
Jvec(4)=-min(est.phi_cvar);                             % phase variance (magnitude of order parameter)
Jvec(5)=-min(pwr);                                      % power

end

