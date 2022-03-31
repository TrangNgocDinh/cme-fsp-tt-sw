function marginal_dist = cme_marginal_traditional(p,state_space_list)
    global no_species

    marginal_dist   = cell(no_species,1);
    for species=1:no_species
        max_species             = max(state_space_list(:,species));
        vec_marginal            = zeros(max_species+1,1);
        for row=1:size(state_space_list,1)
            pos                 = state_space_list(row,species)+1;
            vec_marginal(pos)   = vec_marginal(pos)+p(row);
        end
        marginal_dist{species}  = vec_marginal;
    end
end
