function [expand_ind,bmass] = cme_update_ss_qtt_find_indices(v,sv,t,t_step,marginal_old,marginal_new,marg_diff,is,acc_error,bmass_old)
    global FSP_qtt_size no_species tol_fsp final_time
%   Find boundary mass
    N2          = length(FSP_qtt_size);
    bmass       = zeros(1,N2);
    w2          = reshape(v,2.^FSP_qtt_size);
    fspsize2    = 2.^FSP_qtt_size;
    C2          = cell(N2,1);
    for i2=1:N2
        C2{i2}  = core(w2,i2);
    end
    D2          = cell(N2,1);
    for i2=1:N2
        F2          = C2{i2};
        [r1,~,r2]   = size(F2);
        E2          = zeros(r1,r2);
        for j2=1:fspsize2(i2)
            G2      = reshape(F2(:,j2,:),[r1 r2]);
            E2      = E2+G2;
        end
        D2{i2}      = E2;
    end
%   Contract to find the N marginal distributions
    marginal2   = cell(N2,1);
    for i2=1:N2
        p3      = zeros(fspsize2(i2),1);
        E2      = 1;
        for j2=1:i2-1
            E2  = E2*D2{j2};
        end
        if (i2<N2)
            F2  = D2{i2+1};
        else
            F2  = 1;
        end
        for j2=i2+2:N2
            F2  = F2*D2{j2};
        end
        G2  = C2{i2};
        [r1,~,r2]   = size(G2);
        for j2=1:fspsize2(i2)
            p3(j2)  = E2*reshape(G2(:,j2,:),[r1 r2])*F2;
        end
        numax2          = max(sv(i2,:));
        bmass(i2)       = sum(p3(fspsize2(i2)-numax2+1:fspsize2(i2)));
        marginal2{i2}   = p3;
        for j2=1:no_species
            for k2=1:size(marginal2{j2},1)
                marginal_new(j2,k2) =  marginal2{j2}(k2);
            end
        end
        if (size(marginal_new,1)==size(marginal_old,1)&&(size(marginal_new,2)==size(marginal_old,2)))
            marg_diff   = 0;
            for k2=1:size(marginal_new,1)
                for h2=1:size(marginal_new,2)
                    marg_diff   = marg_diff+abs(marginal_old(k2,h2)-marginal_new(k2,h2));
                end
            end
        end
        marginal_old    = marginal_new;
        difftracker(is) = marg_diff;
    end
%   Find edges with a lot of mass to expand
    expand_ind = (no_species*(acc_error+(bmass-bmass_old))>tol_fsp*((t+t_step)/final_time));
end
