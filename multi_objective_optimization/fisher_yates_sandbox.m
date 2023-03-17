clear
addpath('utils/')
Y=rand(1e3,8,2);
Nperm=1e3;
tic
get_permuted_pages_recursive(repmat(Y,1,1,Nperm),Nperm);
toc

%%
randperm_recursive(NaN(5,5),1,5,5)

function Ypages = get_permuted_pages_recursive(Ypages,ind)
N=size(Ypages,1);
d=size(Ypages,2);
Nperm=size(Ypages,3);
Ypages=NaN(N,d,Nperm);
R=randperm_recursive(NaN(N,d),1,d,N);
Ypages(:,:,ind)=Ypages(sub2ind(size(R),repmat((1:N)',1,d),R));
if ind<Nperm
    Ypages=get_permuted_pages_recursive(Ypages,ind+1);
end
end

function R=randperm_recursive(R,ind,d,N)
R(ind,:)=randperm(d,d);
if ind<N
    R=randperm_recursive(R,ind+1,d,N);
end 
end

function Ypages = get_permuted_pages_fast(Y,Nperm)
N=size(Y,1);
d=size(Y,2);

% for constructing FY sequences
modf=@(a,b) (a==0).*(b+a)+a;
Ypages=NaN(N,d,Nperm);
for pp=1:Nperm
    A=Y;

    % construct Fisher Yates sequences
    r=randi(d,[N,d-1]);
    for ii=1:(d-2)
        r(:,ii+1)=modf(mod(r(:,ii+1),d-ii),d-ii); 
    end
    r(:,end+1)=1;

    source =sub2ind(size(A),repmat((1:N)',1,d),repmat(1:d,N,1));
    target =sub2ind(size(A),repmat((1:N)',1,d),repmat(d+1,N,1)-r(:,1:d));
    Atmp=A(target);
    A(target)=A(source);
    A(source)=Atmp;
    
    Ypages(:,:,pp)=A;
end
end

function Ypages = get_permuted_pages(Y,Nperm)
N=size(Y,1);
d=size(Y,2);

% for constructing FY sequences
modf=@(a,b) (a==0).*(b+a)+a;
Ypages=NaN(N,d,Nperm);
for pp=1:Nperm
    A=Y;

    % construct Fisher Yates sequences
    r=randi(d,[N,d-1]);
    for ii=1:(d-2)
        r(:,ii+1)=modf(mod(r(:,ii+1),d-ii),d-ii); 
    end

    for ii=1:(d-1)
        source =sub2ind(size(A),(1:N)',repmat(ii,N,1));
        target =sub2ind(size(A),(1:N)',repmat(d+1,N,1)-r(:,ii));
        Atmp=A(target);
        A(target)=A(source);
        A(source)=Atmp;
    end
    Ypages(:,:,pp)=A;
end
end

% %%
% 
% for ii=1:(d-1)
%     Atmp=A(r(:,ii))
%     %%
%     A(r(1,ii))=A(ii);
%     A(ii)=Atmp;
%     A
% end
