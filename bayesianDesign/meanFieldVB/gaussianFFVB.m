% assume posterior is sum of mutlivariate normal distributions
vec =@(A) reshape(A,size(A,1)*size(A,2),1); % d^2 column vector


function Avech=vech(A)
d=size(A,1);
Avech=[];
for ii=1:size(A,2)
    for jj=ii:d
        Avech(end+1)=A(jj,ii);
    end
end
end
function S=sample_from_q(Nsamp,muCell,cholMatCell,wvec)
% input:
% muCell    - list of means
% sigmaCell - list of covariance matrices
% wvec      - weights for each gaussian 
% output: samples from the summed distribution
% generate random samples from qlambda(theta)
mu_inds = cumsum(mnrnd(Nsamp,wvec));
S=randn([Nsamp,3]);

for ii=1:length(mu_inds)
    R=cholMatCell{ii};
    if ii>1
        Sstart=mu_inds(ii-1)+1;
    else
        Sstart=1;
    end
    Sstop=mu_inds(ii);
    S(Sstart:Sstop)=muCell{ii} + S(Sstart:Sstop)*R;
end
end