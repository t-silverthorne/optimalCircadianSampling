function [amp_est,acro_est,per_est]=convertToCircularParams(c,model)
    switch model
        case 'cosinorOneFreq'
            per_est = c(3);
            acro_est=sincos2acro(c(1),c(2));
            amp_est =sincos2amp(c(1),c(2));
    end
    
    function acro=sincos2acro(sin,cos)
    
    acro=mod(atan2(sin,cos),2*pi);
    end
    
    function amp=sincos2amp(sin,cos)
    amp=sqrt(sin^2+cos^2);
    end
end
