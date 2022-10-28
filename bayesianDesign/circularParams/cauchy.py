import numpy as np
from scipy.stats import norm
N=100 # desired number of samples
T=10 # length of Markov chain
k=3
sig=2

def p(x):
 return	np.prod(1./(1+pow(x/sig,2)),1) # column-wise product


X=np.random.rand(N,k,T)+3
Y=np.random.rand(N,k,T)
for t in range(T):	
	Z=Y[:,:,t]+X[:,:,t-1]
	diffp=np.log(p(Z))-np.log(p(X[:,:,t-1]))
	logdiff=np.prod(norm.pdf(X[:,:,t-1],Z),1) - np.prod(norm.pdf(Z,X[:,:,t-1]),1)
	transprob=diffp+logdiff
	alpha=np.multiply(transprob<0,transprob)
	rmat=np.log(np.random.rand(N))
	X[:,:,t]= np.multiply((rmat<alpha),Y[:,:,t]) + X[:,:,t-1]

print(X[:,:,T-1])
