function Jvec = costfun_power_bias_var(param,nodes)
% assuming we want to minimize J1 = bias + var - min power
% J2 is a weaker cost function, only containing terms that you can
% calculate without knowing the (in principal unknown) parameter values
[acrovec,pwr,est]=wrap_simulatePWR_matperm_fv(param,nodes);


% J1=max(abs(param.Amp-est.amp_mu)) + max(abs(exp(1j*acrovec)-exp(1j*est.phi_mu))) + ...
%     + max(est.amp_st) - min(est.phi_cvar) -...
%       min(pwr);
% J2 = max(est.amp_st) - min(est.phi_cvar) - min(pwr);
 
Jvec(1)=max(abs(param.Amp-est.amp_mu));                 % amplitude error
Jvec(2)=max(abs(exp(1j*acrovec)-exp(1j*est.phi_mu)));   % phase error (angle of order parameter)
Jvec(3)=max(est.amp_st);                                % amp stdev
Jvec(4)=-min(est.phi_cvar);                             % phase variance (magnitude of order parameter)
Jvec(5)=-min(pwr);                                      % power

end

