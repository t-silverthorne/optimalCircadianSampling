function Y = cosinorOneFreq(t,theta)
A1=theta.A1; % extract params
phi1=theta.phi1;
T1=theta.T1;
Y=A1.*cos(2*pi*t./T1-phi1);
end