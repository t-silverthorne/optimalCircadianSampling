function p = evalLogLikelihood(t_obs,Y_obs,theta,model)
p=NaN;
switch model
    case 'cosinorTwoFreq'
        p= -sum((Y_obs - cosinorTwoFreq(t_obs,theta)).^2/2);
    case 'cosinorOneFreq'
        p= -sum((Y_obs - cosinorOneFreq(t_obs,theta)).^2/2);
end
end

