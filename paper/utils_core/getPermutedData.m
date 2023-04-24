function YI = getPermutedData(Y,R,I3,I4)
%%%%%%%%%%%%%%
%GETPERMUTEDDATA act on simulated data Y with permutations stored in R
% INPUT:   
%   Y          simulated data (size Nresidual x Nmeas x 1 x Nacro)
%   R          permutations   (size Nresidual x Nmeas x Nperm x Nacro)
%   I3         bookkeeping matrix for sub2ind, generated in constructUtilMats
%   I4         bookkeeping matrix for sub2ind, generated in constructUtilMats
% OUTPUT:
%   YI         permuted data (size Nresidual x Nmeas x Nperm x Nacro)
%%%%%%%%%%%%%%

% method 1
N1=size(R,1);N2=size(R,2);N3=size(R,3);N4=size(R,4);
Y=repmat(Y,1,1,N3,1);

method='m1';
switch method
    case 'm1'
        YI=Y(sub2ind(size(Y),repmat((1:N1)',1,N2,N3,N4),R,repmat(I3,1,1,1,N4),I4));
    case 'm2'
        d=N2;
        I=eye(d);
        I=repmat(I,[1 1 N3 N4]);
        ind3=repmat(ones(d,d),[1 1 N3 N4]);
        for ii=1:N3
            ind3(:,:,ii,:)=ind3(:,:,ii,:)*ii;
        end
        ind4=repmat(ones(d,d),[1 1 N3 N4]);
        for ii=1:N4
            ind4(:,:,:,ii)=ind4(:,:,:,ii)*ii;
        end
        Yloc=Y;
        YI=NaN(N1,N2,N3,N4); % construct in transpose form
        for ii=1:N1
            Rloc=R(ii,:,:,:);
            YI(ii,:,:,:)=pagemtimes(Yloc(ii,:,:,:),I(sub2ind(size(I),repmat((1:d)',[1 d N3 N4]),repmat(Rloc,[d 1 1 1]),ind3,ind4)));
        end
end

end
