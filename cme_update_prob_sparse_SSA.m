function p_new = cme_update_prob_sparse_SSA(p_old,state_space_list_old,state_space_list_new)

    p_new   = zeros(size(state_space_list_new,1),1);
    for i=1:length(p_new)
        state   = state_space_list_new(i,:);
        for j=1:length(p_old)
            state_old   = state_space_list_old(j,:);
            if isequal(state,state_old)
                p_new(i)    = p_old(j);
                break
            end
        end
    end
end
