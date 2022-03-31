function [vec_ss_minsize_new,vec_ss_maxsize_new] = cme_update_ss_qttmw_reduce(vec_ss_minsize_old,vec_ss_maxsize_old,p_current)
    global tol_reduce_ss_tt no_species
%---------------------------------------------------------------Method 1
    % no_states           = prod(vec_ss_maxsize_old-vec_ss_minsize_old+1);
    %
    % vec_ss_minsize_new  = inf(1,no_species);
    % vec_ss_maxsize_new  = zeros(1,no_species);
    % for i=1:no_states
    %     state   = index_to_state(i,vec_ss_minsize_old,vec_ss_maxsize_old);
    %     pos     = state-vec_ss_minsize_old+1;
    %     prop    = p_current(pos);
    %     if prop>tol_reduce_ss_tt
    %         for species=1:no_species
    %             if vec_ss_maxsize_new(species)<state(species)
    %                 vec_ss_maxsize_new(species) = state(species);
    %             end
    %             if vec_ss_minsize_new(species)>state(species)
    %                 vec_ss_minsize_new(species) = state(species);
    %             end
    %         end
    %     end
    % end
    % fprintf('State space after 1st kind of reduction:\n');
    % disp(vec_ss_minsize_new)
    % disp(vec_ss_maxsize_new)
%---------------------------------------------------------------Method 2
%     vec_ss_minsize_new  = vec_ss_minsize_old;
%     vec_ss_maxsize_new  = vec_ss_maxsize_old;
%     for species=1:no_species
%         vec_minsize_temp            = vec_ss_minsize_new;
%         vec_maxsize_temp            = vec_ss_maxsize_new;
% %       Record the current bounds of this species
%         species_minsize             = vec_minsize_temp(species);
%         species_maxsize             = vec_maxsize_temp(species);
% %       Record the current bounds of other species
%         vec_minsize_temp(species)   = [];
%         vec_maxsize_temp(species)   = [];
% %       Number of states to consider when reducing the bounds for this species
%         no_states_temp              = prod(vec_maxsize_temp-vec_minsize_temp+1);
% %       Increasing minsize for this species...
%         logic_increase_min          = 1;
%         while logic_increase_min==1
% %           For every possible state to consider...
%             for i=1:no_states_temp
%                 vec_others      = index_to_state(i,vec_minsize_temp,vec_maxsize_temp);
%                 vec_whole       = zeros(1,no_species);
%                 vec_whole(1:species-1)  = vec_others(1:species-1);
%                 vec_whole(species)      = species_minsize;
%                 vec_whole(species+1:no_species) = vec_others(species:no_species-1);
% %               Find if the minsize state exceeds the tolerance...
%                 pos                     = vec_whole-vec_ss_minsize_old+1;
%                 prop                    = p_current(pos);
%                 if prop>tol_reduce_ss_tt
%                     logic_increase_min  = -1;
%                     break;
%                 end
%             end
%             if logic_increase_min==1
%                 species_minsize = species_minsize+1;
%             end
%         end
%         vec_ss_minsize_new(species)   = species_minsize;
% %       Decreasing maxsize for this species...
%         logic_decrease_max          = 1;
%         while logic_decrease_max==1
% %           For every possible state to consider...
%             for i=1:no_states_temp
%                 vec_others      = index_to_state(i,vec_minsize_temp,vec_maxsize_temp);
%                 vec_whole       = zeros(1,no_species);
%                 vec_whole(1:species-1)  = vec_others(1:species-1);
%                 vec_whole(species)      = species_maxsize;
%                 vec_whole(species+1:no_species) = vec_others(species:no_species-1);
% %               Find if the maxsize state exceeds the tolerance...
%                 pos                     = vec_whole-vec_ss_minsize_old+1;
%                 prop                    = p_current(pos);
%                 if prop>tol_reduce_ss_tt
%                     logic_decrease_max  = -1;
%                     break;
%                 end
%             end
%             if logic_decrease_max==1
%                 species_maxsize = species_maxsize-1;
%             end
%         end
%         vec_ss_maxsize_new(species)   = species_maxsize;
%     end
%     fprintf('State space after method 2:\n');
%     disp(vec_ss_minsize_new)
%     disp(vec_ss_maxsize_new)
%---------------------------------------------------------------Method 3
%     vec_ss_minsize_new  = zeros(1,no_species);
%     vec_ss_maxsize_new  = zeros(1,no_species);
% %   Find the marginal distributions
%     fspsize = vec_ss_maxsize_old-vec_ss_minsize_old+1;
%     C       = cell(no_species,1);
%     for i=1:no_species
%         C{i}    = core(p_current,i);
%     end
%     D           = cell(no_species,1);
%     for i=1:no_species
%         F           = C{i};
%         [r1,n1,r2]  = size(F);
%         E           = zeros(r1,r2);
%         for j=1:fspsize(i)
%             G   = reshape(F(:,j,:),[r1 r2]);
%             E   = E+G;
%         end
%         D{i}    = E;
%     end
%     marginal    = cell(no_species,1);
%     for i=1:no_species
%         p   = zeros(fspsize(i),1);
%         E   = 1;
%         for j=1:i-1
%             E   = E*D{j};
%         end
%         if (i<no_species)
%             F   = D{i+1};
%         else
%             F   = 1;
%         end
%         for j=i+2:no_species
%             F   = F*D{j};
%         end
%         G           = C{i};
%         [r1,n1,r2]  = size(G);
%         for j=1:fspsize(i)
%             p(j)    = E*reshape(G(:,j,:),[r1 r2])*F;
%         end
%         marginal{i} = p;
%     end
% %   Find the new bounds for each species based on its marginal dist.
%     for species=1:no_species
% %       Old bounds for this species
%         minsize_old = vec_ss_minsize_old(species);
%         maxsize_old = vec_ss_maxsize_old(species);
% %       Extract the marginal distribution for this species
%         dist_marg   = marginal{species};
% %       Find the bounds that satisfy the tolerance
%         pos         = find(dist_marg>tol_reduce_ss_tt,1,'first');
%         vec_ss_minsize_new(species) = minsize_old+pos-1;
%         pos         = find(dist_marg>tol_reduce_ss_tt,1,'last');
%         vec_ss_maxsize_new(species) = minsize_old+pos-1;
%     end
% %---------------------------------------------------------------Method 4
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
%     end
%     if zero_p_current == 1
%         fprintf('p current is identically zeros.\n');  % p_current is identically zero
%     end
%---------------------------------------------------------------Method 5
    [vec_ss_minsize_new,vec_ss_maxsize_new,~] = cme_marginal_qttmw(vec_ss_minsize_old,vec_ss_maxsize_old,p_current);
end

function state = index_to_state(index,vec_minsize,vec_maxsize)
%---Check if the index is within the bounds
    no_total_states = prod(vec_maxsize-vec_minsize+1);
    if (index<1)|(index>no_total_states)
        state   = -ones(1,length(vec_minsize));
        return
    end
%---Conversion from index to state
    state   = zeros(1,length(vec_minsize));
    K       = index;
    for i=length(vec_minsize):-1:1
        C   = 1;
        for j=1:i-1
            C   = C*(vec_maxsize(j)-vec_minsize(j)+1);
        end
        % state(i)    = idivide(int64(K),int64((C+vec_minsize(i))));
        state(i)    = idivide(int64(K),int64((C)));
        % if (mod(K,(C+vec_minsize(i)))==0)
        if (mod(K,(C))==0)
            state(i)    = state(i)-1;
        end
        % K           = K-(state(i)-vec_minsize(i))*C;
        K           = K-(state(i))*C;
    end
    state   = state+vec_minsize;
end
