%% (Nsamp,p)
p=@(x1,x2,x3)
%%
Nsamp=10;
randn_array=randn(Nsamp*T,4);
T=1e2; % length of Markov chain
Sampvec=NaN(Nsamp,8);
for N=1:Nsamp
    Xt=rand(1,4);
    for t=1:T
        Y=randn_array(N*(T-1)+t,:)+Xt;
        alpha=min( p(Y)*normpdfvec(Xt,Y)/p(Xt)/normpdfvec(Y,Xt),1 );
        if rand<alpha
            Xt=Y;
        end
    end
    Sampvec(N,:)=Xt;
end
%%
normpdfvec([1 1 1 1],[1 1 1 1])
%%
function p = normpdfvec(X,Y)
p=prod(arrayfun(@(ind) normpdf(X(ind),Y(ind)),1:numel(X)));
end