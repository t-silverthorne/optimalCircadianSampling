syms x mu sig S1 S2 S3
c=fliplr(coeffs(collect((x-mu)^2/2/sig^2 + (S1*x^2 +2*x*(S2-S3))/2,x),x))

m=-c(2)/2/c(1); n=c(3) - (c(2)^2/4/c(1));
%%
simplify(c(1)*(x-m)^2+n +(  -(x-mu)^2/2/sig^2 - (S1*x^2 +2*x*(S2-S3))/2 ))
%%
sigNew=1/sqrt(2*c(1))
%%
simplify((x-m)^2/2/sigNew^2+n +(  -(x-mu)^2/2/sig^2 - (S1*x^2 +2*x*(S2-S3))/2 ))
