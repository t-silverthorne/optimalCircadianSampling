function Y = getPermDataAllMethods(p,Y_in)
% options: naive, naive_make_perms_first, naive_reuse_perms,
%permMethod       = 'naive';
%permMethod       = 'naive_reuse_perms'; 
%permMethod       = 'exact_in_batches';
%permMethod       = 'FY_double_for_loop'; 

switch p.permMethod
    case 'naive'
    % loop over all data and generate permutations one at a time (slowest method)
        Y = repmat(Y_in,[1 1 p.Nperm 1]);
        for kk=1:p.Nacro
            for jj=1:p.Nresidual
                for ii=1:p.Nperm
                    Y(jj,:,ii,kk)=Y(jj,randperm(p.Nmeas,p.Nmeas),ii,kk);
                end
            end
        end

    case 'exact_in_batches'
    % use entire permutation group, only works for small Nmeas
        pind1 = p.perm_start_ind;
        pind2 = p.perm_stop_ind;
        Y = repmat(Y_in,[1 1 p.Nperm 1]);
        all_perms=perms(1:p.Nmeas);
        for ii=pind1:pind2
            rp=all_perms(ii,:);
            I1=repmat((1:p.Nresidual)',[1 p.Nmeas,1,p.Nacro]);
            I2=repmat(rp,[p.Nresidual 1 1 p.Nacro]);
            I3=ii*ones([p.Nresidual,p.Nmeas,1,p.Nacro]);
            I4=reshape(1:p.Nacro,[1 1 1 p.Nacro]).*ones([p.Nresidual,p.Nmeas,1,p.Nacro]);
            Y(:,:,ii,:) = Y(sub2ind(size(Y),I1,I2,I3,I4));
        end

    case 'naive_reuse_perms'
    % generate one batch of permutations and reuse it on all residuals and acrophases
    % likely that this will introduce bias since samples are shared
        Y = repmat(Y_in,[1 1 p.Nperm 1]);
        for ii=1:p.Nperm
            rp=randperm(p.Nmeas,p.Nmeas);
            I1=repmat((1:p.Nresidual)',[1 p.Nmeas,1,p.Nacro]);
            I2=repmat(rp,[p.Nresidual 1 1 p.Nacro]);
            I3=ii*ones([p.Nresidual,p.Nmeas,1,p.Nacro]);
            I4=reshape(1:p.Nacro,[1 1 1 p.Nacro]).*ones([p.Nresidual,p.Nmeas,1,p.Nacro]);
            Y(:,:,ii,:) = Y(sub2ind(size(Y),I1,I2,I3,I4));
        end

    case 'FY_double_for_loop'
    % implementation of the Fisher-Yates shuffle, no reusing of
    % permutations
        Y = repmat(Y_in,[1 1 p.Nperm 1]);
        N1=p.Nresidual;d=p.Nmeas;N2=d;N3=p.Nperm;N4=p.Nacro;
        n=N1*N3*N4;
        I3=NaN(N1,N2,N3); % make indexing matrices
        I4=NaN(N1,N2,N3,N4);
        for ii=1:N3
            I3(:,:,ii)=ii*ones(N1,N2);
        end
        for ii=1:N4
            I4(:,:,:,ii)=ii*ones(N1,N2,N3);
        end
        sqr
        P=repmat(1:d,n,1); % make permutation matrix using FY shuffle
        P=P';
        for ii=1:d-1
            avec=ii-1+randi(d+1-ii,[n 1]);
            for jj=1:n
                Ytemp=P(ii,jj);
                P(ii,jj)=P(avec(jj),jj);
                P(avec(jj),jj)=Ytemp;
            end
        end
        P=P';
        P=pagetranspose(reshape(P',d,N1,N3,N4));
        % reindex to act with permutation
        switch p.permActionMethod
            case 'index'
                % reindex matrix using P the FY permutations
                Y=Y(sub2ind(size(Y),repmat((1:N1)',1,N2,N3,N4),P,repmat(I3,1,1,1,N4),I4));
            case 'matrix'
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
                    Ploc=P(ii,:,:,:);
                    YI(ii,:,:,:)=pagemtimes(Yloc(ii,:,:,:),I(sub2ind(size(I),repmat((1:d)',[1 d N3 N4]),repmat(Ploc,[d 1 1 1]),ind3,ind4)));
                end
        end

    
end
end

