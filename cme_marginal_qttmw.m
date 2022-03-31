function [vec_ss_minsize_new,vec_ss_maxsize_new,marginal_reduced] = cme_marginal_qttmw(vec_ss_minsize_old,vec_ss_maxsize_old,p_current)
    global tol_reduce_ss_tt no_species
% %==============================================================Version 1
%     size_old = vec_ss_maxsize_old - vec_ss_minsize_old + 1;
%     vec_ss_minsize_new = vec_ss_minsize_old;
%     vec_ss_maxsize_new = vec_ss_maxsize_old;
%
%     marginal_reduced = cell(no_species,1);
%
%     left_mode = ones(1,no_species+1);
%     right_mode = ones(1,no_species+1);
%     right_mode(1) = prod(size_old)/size_old(1);
%     for i=2:no_species
%         left_mode(i) = left_mode(i-1)*size_old(i-1);
%         right_mode(i) = right_mode(i-1)/size_old(i);
%     end
%
%     zero_p_current = -1; % Suppose that p_current is not identically zero
%
%     for i=1:no_species
%         p_reshape = tt_reshape(p_current,[left_mode(i) size_old(i) right_mode(i)]);
%         count_up = 0;
%         while (tt_max_abs(p_reshape(:,1+count_up,:)) < tol_reduce_ss_tt)
%             if 1+count_up < size_old(i)
%                 count_up = count_up + 1;
%             else
%                 zero_p_current = 1; % p_current is identically zero
%                 break;
%             end
%         end
%         vec_ss_minsize_new(i) = vec_ss_minsize_old(i) + count_up;
%         count_down = 0;
%         while (tt_max_abs(p_reshape(:,size_old(i)-count_down,:)) < tol_reduce_ss_tt)
%             if size_old(i)-count_down > 1
%                 count_down = count_down + 1;
%             else
%                 zero_p_current = 1; % p_current is identically zero
%                 break;
%             end
%         end
%         vec_ss_maxsize_new(i) = vec_ss_maxsize_old(i) - count_down;
%         size_new(i) = vec_ss_maxsize_new(i) - vec_ss_minsize_new(i) + 1;
%
%         marginal = zeros(size_new(i),1);
%         for j=1:size_new(i)
%             marginal(j) = marginal(j) + sum(p_reshape(:,j+count_up,:));
%         end
%         marginal_reduced{i} = marginal;
%     end
%     if zero_p_current == 1
%         fprintf('p current is identically zeros.\n');  % p_current is identically zero
%     end
%==============================================================Version 2
%   Generate the marginal distributions

    size_old = vec_ss_maxsize_old - vec_ss_minsize_old + 1; % It is 2.^FSP_qtt_size of Huy
    fspsize = vec_ss_maxsize_old - vec_ss_minsize_old + 1; % fspsize of Huy
    % Do not need. Bin remove it later
    vec_ss_minsize_new = vec_ss_minsize_old;
    vec_ss_maxsize_new = vec_ss_maxsize_old;
    % vec_ss_minsize_new = vec_ss_minsize_old; % Do not need now
    % vec_ss_maxsize_new = vec_ss_maxsize_old; % Do not need now
    FSP_qtt_size = log2(fspsize); % Our l2size_now, we should transfer it to our program
    % marginal_reduced = cell(no_species,1); % It will declare later

    % w       = reshape(v,2.^FSP_qtt_size);
    w       = reshape(p_current,fspsize);
    % fspsize = 2.^FSP_qtt_size; % Do not need
    C       = cell(no_species,1);
    for i=1:no_species
        C{i}    = core(w,i);
    end
    D           = cell(no_species,1);
    for i=1:no_species
        F           = C{i};
        [r1,n1,r2]  = size(F);
        E           = zeros(r1,r2);
        for j=1:fspsize(i)
            G   = reshape(F(:,j,:),[r1 r2]);
            E   = E+G;
        end
        D{i}    = E;
    end
    % marginal    = cell(no_species,1);
    marginal_reduced    = cell(no_species,1);
    for i=1:no_species
        p   = zeros(FSP_qtt_size(i),1);
        E   = 1;
        for j=1:i-1
            E   = E*D{j};
        end
        if (i<no_species)
            F   = D{i+1};
        else
            F   = 1;
        end
        for j=i+2:no_species
            F   = F*D{j};
        end
        G           = C{i};
        [r1,n1,r2]  = size(G);
        for j=1:fspsize(i)
            p(j)    = E*reshape(G(:,j,:),[r1 r2])*F;
        end
        % p



        count_up    = 0;
        while (sum(p(1:1+count_up))<tol_reduce_ss_tt)
            if (1+count_up)<size_old(i)
                count_up    = count_up+1;
            else
                zero_p_current = 1; % p_current is identically zero
                break;
            end
        end
        count_down  = 0;
        while (sum(p(end-count_down:end))<tol_reduce_ss_tt)
            if length(p)-count_down>1
                count_down      = count_down+1;
            else
                zero_p_current = 1; % p_current is identically zero
                break;
            end
        end
        vec_ss_minsize_new(i)   = vec_ss_minsize_new(i)+count_up;
        vec_ss_maxsize_new(i)   = vec_ss_maxsize_new(i)-count_down;
        p                       = p(1+count_up:end-count_down);


        % marginal{i} = p;
        marginal_reduced{i} = p;
    end
end
