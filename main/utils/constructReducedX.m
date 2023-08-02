
function X=constructReducedX(t,freq)
% build reduced design matrix X
x1=cos(2*pi*freq*t);
x2=sin(2*pi*freq*t);
X= [x1' x2'];
end
