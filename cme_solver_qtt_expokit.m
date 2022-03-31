function [w,marginal] = cme_solver_qtt_expokit(p_0,propen_func_partial)
    global tol_expokit_tt deg_expokit_tt final_time stoich_mat_total
    global no_species FSP_qtt_size t_step_initial_qtt_expokit

    is              = 1;
    iter            = 1;
    marginal_old    = zeros(1,1,1,1);
    marginal_new    = marginal_old;
    marg_diff       = 10000;
    sv              = stoich_mat_total';
    bmass_old       = zeros(1,no_species);
    acc_error       = zeros(1,no_species);
%   Initial probability vector in QTT format
    v               = p_0;
%   CME matrix in QTT format
    A               = cme_generator_qtt(propen_func_partial);

    t_step          = t_step_initial_qtt_expokit;

    t               = 0;
    while (t<final_time)
        fprintf('At time %f...\n',t);
        iter                = iter+1;
        t_step              = min(t_step,final_time-t);
%       Find the indices to expand the state space if necessary
        [expand_ind,bmass]  = cme_update_ss_qtt_find_indices(v,sv,t,t_step,marginal_old,marginal_new,marg_diff,is,acc_error,bmass_old);
%       Expand the state space if necessary
        if (sum(expand_ind)>0)
%           FSP criteria failed at least at one dimension, expand and redo the step
            v                   = cme_update_ss_qtt_expand(v,FSP_qtt_size,expand_ind);
            FSP_qtt_size        = FSP_qtt_size + expand_ind;
%           Re-generate the matrix and the one-step linear system
            fprintf('Expand A-qtt...\n');
            A                   = cme_generator_qtt(propen_func_partial);
            disp('Generate A complete.');

            bmass_old(expand_ind) = 0;

            disp('FSP criteria fails, expand');
            disp(FSP_qtt_size);
        else
%           Find the next probability tensor by using TT-Expokit
            [w,err]     = tt_expv(t_step,A,v,tol_expokit_tt,deg_expokit_tt);
%           Advance to the next step
            v           = w;
            t           = t + t_step;

            t_step      = min(t_step,final_time-t);

        end
    end
%   Generate the marginal distributions
    w       = reshape(v,2.^FSP_qtt_size);
    fspsize = 2.^FSP_qtt_size;
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
    marginal    = cell(no_species,1);
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
        marginal{i} = p;
    end




end
