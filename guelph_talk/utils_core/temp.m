% d=8;
% n=1e3;
% P=repmat(1:d,n,1);
% for ii=1:d
%     inds=randsample((ii:d)',n,true);
%     ptemp=P(:,ii);
%     P(:,ii)=P(sub2ind(size(P),(1:n)',inds))
%     P(sub2ind(size(P),(1:n)',inds))=ptemp;
% end
% sum(P,2)
% %%
%%
d=8
tic
for ii=1:10
    randsample((1:d)',1e3,true);
end
toc

tic
randsample((1:d)',1e4,true);
toc

%%
d=8;
n=1e3;
Y=repmat(1:d,n,1);
for kk=1:N4;
    for ii=1:d
        inds=randsample((ii:d)',n,true);
        ptemp=Y(:,ii);
        Y(:,ii)=Y(sub2ind(size(Y),(1:n)',inds))
        Y(sub2ind(size(Y),(1:n)',inds))=ptemp;
    end
end
sum(Y,2)
%%


%%
p.Nresidual=1e3;
p.Nperm=1e2;
p.Nacro=16;
p.Nmeas=12;
Y_in=randn(p.Nresidual,p.Nmeas,1,p.Nacro);
p.permMethod='FY_randsample2';
tic
getPermDataAllMethods(p,Y_in);
toc

p.permMethod='FY_double_for_loop';
p.permActionMethod='index';
tic
getPermDataAllMethods(p,Y_in);
toc
