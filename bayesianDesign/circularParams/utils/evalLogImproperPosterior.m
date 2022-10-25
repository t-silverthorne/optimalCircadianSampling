function p = evalLogImproperPosterior(t_obs_MAT,Y_obs_MAT,theta,model,method,settings)
if size(t_obs_MAT,1)>1
    p=evalLogLikelihood(t_obs_MAT(1,:),Y_obs_MAT(1,:),theta,model) + ...
       evalLogImproperPosterior(t_obs_MAT(2:end,:),Y_obs_MAT(2:end,:),theta,model,method,settings);
else
    p=evalLogLikelihood(t_obs_MAT,Y_obs_MAT,theta,model)+log(evalPrior(theta,model,method,settings));
end
end

