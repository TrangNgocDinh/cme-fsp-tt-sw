function [w,state_space_list,vec_ss_minsize,vec_ss_maxsize] = cme_solver_SSA_expokit(propen_func)
    global vec_ss_minsize_SSA_initial vec_ss_maxsize_SSA_initial
    global initial_state final_time t_step_initial_SSA
    global tol_expokit_SSA deg_expokit_SSA tol_local_fsp_SSA
%   Build the CME matrix according to the initial state space
    vec_ss_minsize          = vec_ss_minsize_SSA_initial;
    vec_ss_maxsize          = vec_ss_maxsize_SSA_initial;
    [A,state_space_list]    = cme_generator_SSA_expokit(vec_ss_minsize,vec_ss_maxsize,propen_func);
%   Build the initial probability vector according to the state space
    p_0                     = cme_initial_traditional_sparse(state_space_list,initial_state);
    v                       = p_0;

    t                       = 0;
    t_step                  = t_step_initial_SSA;

    iter                    = 0;

    while (t<final_time)
        iter        = iter+1;
%       Find the tolerance on the FSP error
        tol_EXPV    = tol_local_fsp_SSA*(t+t_step)/final_time;
        tol_real    = min(tol_EXPV,tol_expokit_SSA)
%       Find the next probability vector by using EXPOKIT
        p_out       = cme_solver_expokit(t_step,A,v,tol_real,deg_expokit_SSA);
%       Find the FSP error
        error_fsp   = 1-sum(p_out);
%       Check if the FSP error satisfies the local tolerance...
        if (error_fsp>tol_local_fsp_SSA*(t+t_step)/final_time)
%           If it doesn't, then expand the state space...
            error_fsp

            [vec_ss_minsize,vec_ss_maxsize] = cme_update_ss_sparse_SSA(vec_ss_minsize,vec_ss_maxsize,t_step,propen_func);
            [A,state_space_list_new]        = cme_generator_SSA_expokit(vec_ss_minsize,vec_ss_maxsize,propen_func);
            fprintf('At time %f: expanding the state space from %d states to %d states...\n',t,size(state_space_list,1),size(state_space_list_new,1));
            fprintf('New state space:\n');
            disp(vec_ss_minsize)
            disp(vec_ss_maxsize)



%           and change the probability vector to the new state space...
            v   = cme_update_prob_sparse_SSA(v,state_space_list,state_space_list_new);
            state_space_list                = state_space_list_new;
        else
%           If it does, then advance to the next step...
            v   = p_out;
            t   = t+t_step;
%           and reduce the state space...
            [vec_ss_minsize,vec_ss_maxsize] = cme_update_ss_sparse_reduce(vec_ss_minsize,vec_ss_maxsize,v,state_space_list);
            [A,state_space_list_new]        = cme_generator_SSA_expokit(vec_ss_minsize,vec_ss_maxsize,propen_func);
            fprintf('At time %f: reducing the state space from %d states to %d states...\n',t,size(state_space_list,1),size(state_space_list_new,1));
            v   = cme_update_prob_sparse_SSA(v,state_space_list,state_space_list_new);
            state_space_list                = state_space_list_new;
        end
    end
    w   = v;
end
