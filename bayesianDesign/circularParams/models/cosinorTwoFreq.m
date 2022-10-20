function Y = cosinorTwoFreq(t,theta)
A1=theta.A1;A2=theta.A2; % extract params
phi1=theta.phi1;phi2=theta.phi2;
T1=theta.T1;T2=theta.T2;
Y=A1.*cos(2*pi*t./T1-phi1) + A2.*cos(2*pi*t./T2-phi2); 
end

