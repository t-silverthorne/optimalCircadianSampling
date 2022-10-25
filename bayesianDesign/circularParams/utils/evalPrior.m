function p=evalPrior(theta,model,method,settings)
% evaluate prior distribution, used in MCMC sampling of posterior
% method:
%   pseudo-uniform: uniform distribution with requirement that periods are
%   not too similar
%   test-spt: exp distribution on period and amplitudes, uniform for phases
if nargin<4
    settings.dT=.2;
end
p=NaN;
switch model
    case 'cosinorTwoFreq'
        switch method 
            case 'pseudo-uniform'
                p=0;
                if (abs(theta.T1-theta.T2)>settings.dT) && ...
                    (0<theta.phi1 && 2*pi>theta.phi1) && ...
                    (0<theta.phi2 && 2*pi>theta.phi2) && ...
                    (0<theta.A1   && 10>theta.A1) && ...
                    (0<theta.A2   && 10>theta.A2) && ...
                    (0<theta.T1 && 1>theta.T1) && ...
                    (0<theta.T2 && 1>theta.T2)
                    p=1;%(1/10)^2*(1/2/pi)^2*(2/(1-sqrt(2)*settings.dT)^2);
                end
            case 'test-spt'
                T1mean=2/3;
                T2mean=1/4;
                A1mean=10;
                A2mean=5;

                p=exppdf(theta.T1,T1mean)*exppdf(theta.T2,T2mean)*...
                    exppdf(theta.A1,A1mean)*exppdf(theta.A2,A2mean);
        end
    case 'cosinorOneFreq'
        switch method
            case 'test-spt'
                T1mean=2/3;
                A1mean=10;
                p=exppdf(theta.T1,T1mean)*exppdf(theta.A1,A1mean);
        end
end
end

