function [vec_ss_minsize_SSA_new,vec_ss_maxsize_SSA_new] = cme_update_ss_qttmw_SSA(vec_ss_minsize_SSA_old,vec_ss_maxsize_SSA_old,t_step,propen_func);
    global no_SSA_ss_sparse no_species

    vec_ss_minsize_SSA_new  = inf(1,no_species);
    vec_ss_maxsize_SSA_new  = zeros(1,no_species);
    for i_bound=1:2
        if i_bound==1
            state_start = vec_ss_minsize_SSA_old;
        else
            state_start = vec_ss_maxsize_SSA_old;
        end
        for i=1:no_SSA_ss_sparse
            [~,~,state_max,state_min] = ssa(state_start,t_step,propen_func);
            for species=1:no_species
                if vec_ss_maxsize_SSA_new(species)<state_max(species)
                    vec_ss_maxsize_SSA_new(species) = state_max(species);
                end
                if vec_ss_minsize_SSA_new(species)>state_min(species)
                    vec_ss_minsize_SSA_new(species) = state_min(species);
                end
            end
        end
    end
end
