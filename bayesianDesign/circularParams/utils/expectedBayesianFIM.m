function Mexp = expectedBayesianFIM(M,fnames,t_obs_MAT,Y_obs_MAT,model,method,settings)
N=10000;
thetasamp=samplePosteriorMCMC(N,fnames,t_obs_MAT,Y_obs_MAT,model,method,settings);

Mexpvec=NaN(1,N);
tmeasc=num2cell(t_obs_MAT(1,:));
parfor ii=1:N
    thetac=num2cell(thetasamp(ii,:));    
    Mexpvec(ii)=log(det(M(thetac{:},tmeasc{:})));
end
Mexp=mean(Mexpvec);
    %function Mt=evaluateM(M,theta,tmeas)
%
 %   end
end

