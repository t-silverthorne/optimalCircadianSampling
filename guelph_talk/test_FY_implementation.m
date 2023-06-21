addpath('utils_core/')
close all

p.Nperm=1e5;
p.Nmeas=4;
p.Nresidual=1;
p.Nacro=1;
avec=1:p.Nmeas;
A=repmat(avec,[p.Nresidual 1 1 p.Nacro]);

p.permMethod       = 'FY_double_for_loop'; 
p.permActionMethod = 'index';
Ap=transpose(reshape(pagetranspose(getPermDataAllMethods(p,A)),[p.Nmeas p.Nresidual*p.Nacro*p.Nperm]));

B=cell(size(A,1),1);
for ii=1:size(Ap,1)
    B{ii}=convertStringsToChars(strjoin(string(Ap(ii,:)),''));
end

B=categorical(B)
hist(B)
%%
%histogram(B)