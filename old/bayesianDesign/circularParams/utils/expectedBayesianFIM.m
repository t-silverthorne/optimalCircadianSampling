function Mexp = expectedBayesianFIM(M,fnames,tmeas_prop,t_obs_MAT,Y_obs_MAT,model,method,settings)
switch settings.FIM_expectation_method
    case 'fixed'
        low_enough_variance=true;
        N=1000;
    case 'variance'
        low_enough_variance=false;
<<<<<<< HEAD
        var_cut=1e-1; % todo put back
        N=1000;
=======
        var_cut=settings.var_cut;
        N=settings.batch_size;
>>>>>>> 8bf398575ce343a1f7308b1b3de471dcd097206f
end
Mexpvec=[];
while ~low_enough_variance
    if settings.verbose
        disp('running')
    end
    thetasamp=samplePosteriorMCMC(N,fnames,t_obs_MAT,Y_obs_MAT,model,method,settings);
    
    Mexpvecloc=NaN(1,N);
    tmeasc=num2cell(tmeas_prop);
    switch settings.parallel_mode
        case 'parfor'
            Numworkers=inf;
        case 'vectorize'
            Numworkers=0;
    end
    
    parfor (ii=1:N,Numworkers)
        thetac=num2cell(thetasamp(ii,:));    
        Mexpvecloc(ii)=log(det(M(thetac{:},tmeasc{:})));
    end
    Mexp=mean(Mexpvecloc);
    Mexpvec=horzcat(Mexpvec,Mexpvecloc);
    low_enough_variance= ( std(Mexpvec)/sqrt(numel(Mexpvec))<var_cut );
end

end

