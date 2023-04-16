function Y = getSimulatedData(t,p)
%%%%%%%%%%%%%%
%GETSIMULATEDDATA simulate data from cosinor model with Gaussian noise
% INPUT:   everything is stored in struct p (the parameters)
%   t               measurement times
% size of data
%   p.Nmeas         number of measurements in each experiment
%   p.Nacro         number of acrophases to simulate
%   p.Nresidual     number of simulated experiment (for each acrophase)
% osc features
%   p.Amp           amplitude of cosinor model
%   p.freq          frequency of cosinor model
%   p.noise         noise strength in cosinor model
% OUTPUT:
%   Y               Nresidual x Nmeas x 1 x Nacro array of simulated data
%%%%%%%%%%%%%%
eps=randn(p.Nresidual,p.Nmeas,1,p.Nacro);
acrovec=linspace(0,2*pi,p.Nacro+1);
acrovec=acrovec(1:end-1); % get acrophases
acromat=reshape(acrovec,1,1,1,p.Nacro);

Y=p.Amp*cos(2*pi*t*p.freq-acromat)+p.noise*eps;
end

