function X = constructX(t,p)
% build design matrix X
x0=ones(1,length(t));
x1=sin(2*pi*p.freq*t);
x2=cos(2*pi*p.freq*t);
X= [x0' x1' x2'];
end