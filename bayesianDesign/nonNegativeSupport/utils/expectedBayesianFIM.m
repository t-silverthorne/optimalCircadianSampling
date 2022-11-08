function Mexp = expectedBayesianFIM(M,tmeas_prop,t_obs_MAT,Y_obs_MAT,settings)
tvec=reshape(t_obs_MAT,1,numel(t_obs_MAT));
yvec=reshape(Y_obs_MAT,1,numel(Y_obs_MAT));

switch settings.FIM_expectation_method
    case 'fixed'
        low_enough_variance=true;
        N=1000;
    case 'variance'
        low_enough_variance=false;
        var_cut=settings.var_cut;
        N=settings.batch_size;
end

Mexpvec=[];
while ~low_enough_variance
    if settings.verbose
        disp('running')
    end
    thetasamp=samplePosteriorMCMC(N,yvec,tvec,settings);
    
    Mexpvecloc=NaN(1,N);
    tmeasc=num2cell(tmeas_prop);
    
    for ii=1:N
        thetac=num2cell(thetasamp(ii,:));    
        Mexpvecloc(ii)=real(log(det(M(thetac{:},tmeasc{:})))); % TODO make less ad hoc
    end
    Mexp=mean(Mexpvecloc);
    Mexpvec=horzcat(Mexpvec,Mexpvecloc);
    low_enough_variance= ( std(Mexpvec)/sqrt(numel(Mexpvec))<var_cut );
end

end

