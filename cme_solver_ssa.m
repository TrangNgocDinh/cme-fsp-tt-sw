function [marginal_dist,lb,ub] = cme_solver_ssa(propen_func)
    global no_SSA_traj no_species initial_state final_time

    marginal_dist_mat   = zeros(1,no_species);
    lb                  = zeros(1,no_species);
    ub                  = zeros(1,no_species);

    for i=1:no_SSA_traj
        fprintf('SSA trajectory no. %d/%d...\n',i,no_SSA_traj)
        [~,~,~,~,final_state] = ssa(initial_state,final_time,propen_func);
        for species=1:no_species
            pos = final_state(species)+1;
            if pos>size(marginal_dist_mat,1)
                marginal_dist_mat(pos,species)  = 1;
            else
                marginal_dist_mat(pos,species)  = marginal_dist_mat(pos,species)+1;
            end
        end
    end
    marginal_dist_mat   = marginal_dist_mat/no_SSA_traj;
    marginal_dist   = cell(no_species,1);
    for species=1:no_species
        vec_marginal            = marginal_dist_mat(:,species);
        lb(species)             = find(vec_marginal>0,1,'first')-1;
        ub(species)             = find(vec_marginal>0,1,'last')-1;
        marginal_dist{species}  = vec_marginal((lb(species)+1):(ub(species)+1));
        vec_marginal((lb(species)+1):(ub(species)+1));
    end
end
