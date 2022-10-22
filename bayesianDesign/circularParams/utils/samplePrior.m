function [thetaMat,fieldnames] = samplePrior(numParamSets,model,prior,settings)
% thetaMat: matrix whose rows correspond to parameter sets for model
switch model
    case {'cosinorOneFreq','cosinorTwoFreq'}
        fieldnames={'A1','phi1','T1','A2','phi2','T2'};
        switch prior
            case 'pseudo-uniform'
                if nargin<4
                    settings.A1max=100;
                    settings.A1min=1;
                    settings.A2max=settings.A1max;
                    settings.A2min=settings.A1min;
                    settings.T1max=1;
                    settings.T1min=0.1;
                    settings.T2max=settings.T1max;
                    settings.T2min=settings.T1min;
                end
                thetaMat=NaN(numParamSets,6);
                thetaMat(:,[1 4])=10*rand(numParamSets,2); % get amplitudes
                thetaMat(:,[2 5])= 2*pi*rand(numParamSets,2); % get acrophases
                accepted=NaN(numParamSets,2); 
                % final two need rejection sampling
                ind=1;
                while ind<numParamSets
                    samp=rand(numParamSets,2);
                    acceptedloc=samp(abs(samp(:,1)-samp(:,2))>.2,:);
                    accepted(ind:ind+size(acceptedloc,1)-1,:)=acceptedloc;
                    ind=ind+size(acceptedloc,1);
                end
                thetaMat(:,[3 6])=accepted(1:numParamSets,:); % get periods

        end
        if strcmp(model,'cosinorOneFreq')
            thetaMat=thetaMat(:,1:3);
            fieldnames={'A1','phi1','T1'};
        end

end



end
