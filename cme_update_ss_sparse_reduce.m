function [vec_ss_minsize_SSA_new,vec_ss_maxsize_SSA_new] = cme_update_ss_sparse_reduce(vec_ss_minsize_SSA_old,vec_ss_maxsize_SSA_old,p_current,state_space_list)
    global tol_reduce_ss_sparse no_species

    vec_ss_minsize_SSA_new  = inf(1,no_species);
    vec_ss_maxsize_SSA_new  = zeros(1,no_species);
    for i=1:size(state_space_list,1)
        state   = state_space_list(i,:);
        prop    = p_current(i);
        if prop>tol_reduce_ss_sparse
            for species=1:no_species
                if vec_ss_maxsize_SSA_new(species)<state(species)
                    vec_ss_maxsize_SSA_new(species) = state(species);
                end
                if vec_ss_minsize_SSA_new(species)>state(species)
                    vec_ss_minsize_SSA_new(species) = state(species);
                end
            end
        end
    end
end
