function Jvec = getCostFun(p,pwr_mat,amp_mat,acro_mat,active_inds)
acrovec=linspace(0,2*pi,p.Nacro+1);
acrovec=acrovec(1:end-1);
[acro_bias,acro_var] = getAcroStats(acro_mat,acrovec);
[amp_bias,amp_var]   = getAmpStats(amp_mat,p.Amp);
pwr                  = mean(pwr_mat,1);

Jvec(1)=max(abs(amp_bias)); 
Jvec(2)=max(abs(acro_bias));
Jvec(3)=max(amp_var);  
Jvec(4)=-min(acro_var);
Jvec(5)=-min(pwr);                                      
Jvec   =Jvec(active_inds);
end

