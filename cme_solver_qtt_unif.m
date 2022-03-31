function [w,marginal,vec_t,vec_stepsize,vec_ss_size,T_total_elapsed_qtt,t_amen_total_qtt, iter_qtt] = cme_solver_qtt_unif(p_0,propen_func_partial)
    global tol_fsp tol_unif tol_amen sc_factor final_time no_species
    global FSP_qtt_size stoich_mat_total
    global tol_fsp_probsum
%-------------------------Internal statistics---------------------------
    T_total = tic;
    % t_expand_elapsed_qtt         = 0;
    t_amen_total_qtt         = 0;
    % t_elapsed = 0;
%-----------------------------------------------------------------------
    is              = 1;
    marginal_old    = zeros(1,1,1,1);
    marginal_new    = marginal_old;
    iter            = 1;
    difftracker     = zeros(10000,1);
    difftrack       = zeros(10000,1);
    timetrack       = zeros(10000,1);
    marg_diff       = 10000;
    difftracker(is) = marg_diff;
    difftol         = .01;
    difftrack(is)   = marg_diff;
    sv              = stoich_mat_total';
%   Initial probability vector in QTT format
    v               = p_0;
%   CME matrix in QTT format
    A               = cme_generator_qtt(propen_func_partial);
%   Alpha = bound on the sum of propensities
    alpha           = tt_max(-diag(A));
    nstep           = ceil(alpha*final_time/sc_factor);
    t_step          = final_time/nstep;
%   Performance statistics
    memp            = zeros(10000,1);
    memW            = zeros(10000,1);
    cput            = zeros(10000,1);
    fsperr          = zeros(10000,1);
    fsize           = zeros(10000,no_species);
    W_max_rank      = zeros(10000,1);
    p_max_rank      = zeros(10000,1);
    deg_taylor      = zeros(10000,1);
    step_size       = zeros(10000,1);
    maxdA           = zeros(10000,1);
    ncore           = zeros(10000,1);

    t               = 0;
    save_step       = 1;
    bmass_old       = zeros(1,no_species);
    acc_error       = zeros(1,no_species);
%   Prepare vectors of step-size and state space size
    count           = 0;
    vec_t           = zeros(1,1);
    vec_stepsize    = zeros(1,1);
    vec_ss_size     = zeros(1,1);
    tic
    while (t<final_time)
        fprintf('==================================================\n');
        fprintf('For iteration %d\n', iter);
        iter                = iter+1;
%       Find the indices to expand the state space if necessary
        [expand_ind,bmass]  = cme_update_ss_qtt_find_indices(v,sv,t,t_step,marginal_old,marginal_new,marg_diff,is,acc_error,bmass_old);
%       Expand the state space if necessary
        if (sum(expand_ind)>0)
%           FSP criteria failed at least at one dimension, expand and redo the step
            v               = cme_update_ss_qtt_expand(v,FSP_qtt_size,expand_ind);
            FSP_qtt_size    = FSP_qtt_size + expand_ind
%           Re-generate the matrix and the one-step linear system
            fprintf('Expand A-qtt...\n');
            A               = cme_generator_qtt(propen_func_partial);
            disp('Generate A complete.');

            alpha           = tt_max(-diag(A));
            nstep           = ceil(alpha*(final_time-t)/sc_factor);
            t_step          = (final_time-t)/nstep;

            bmass_old(expand_ind) = 0;

            disp('FSP criteria fails, expand');
            disp(FSP_qtt_size);

            fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
            fprintf('At time %f; final time %f\n',t,final_time);
            fprintf('Expansion occurs!\n');
            fprintf('   Time step = %f\n',t_step);
            fprintf('   Theta = %f\n', alpha);
            fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
        else
%           Find the next probability tensor by using Uniformization & AMEn
            [w,W,d,m, t_amen_elapsed]   = unif(A,v,t_step,nstep);
            fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
            fprintf('At time %f; final time %f\n',t,final_time);
            fprintf('No expansion!\n');
            fprintf('   Time step = %f\n',t_step);
            fprintf('   Theta = %f ; m = %d; time for AMEn is %f\n', alpha, m, t_amen_elapsed);
            t_amen_total_qtt    = t_amen_total_qtt   + t_amen_elapsed;
            fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');


% %           Sanity test: if the probability sum is too small, then
% %           reduce the time step and redo
%             prob_sum    = sum(w);
%             tol_probsum = 1-tol_fsp_probsum*((t+t_step)/final_time);
%             if prob_sum<tol_probsum
%                 disp('Probability sum is too low, redo this step with reduced time step.');
%                 fprintf('Probability sum = %f, tolerance = %f, time step = %f\n',prob_sum,tol_probsum,t_step);
%                 t_step  = min(0.5*t_step,final_time-t);
%                 continue;
%             end


%       Record the step-size and state space size
            count               = count+1;
            vec_t(count)        = t;
            vec_stepsize(count) = t_step;
            vec_ss_size(count)  = 2^(sum(FSP_qtt_size));
%           Advance to the next step
            v           = w;
            t           = t + t_step;
%           Initial guess for the next AMEn solve
            W0          = W;
            rW          = rank(W);
            if (rW(d+1)==1)
%               If the rank of the 'degree' core drops to 1, we have happy
%               breakdown...
                disp('Happy breakdown');
                break
            end
            acc_error       = acc_error+(bmass-bmass_old);
            bmass_old       = bmass;
%           Record some statistics
            timetrack(is)   = toc;
            difftrack(is)   = difftracker(is);
            memp(is)        = mem(v);
            memW(is)        = mem(W);
            fsperr(is)      = max(bmass);
            W_max_rank(is)  = max(rank(W));
            p_max_rank(is)  = max(rank(v));
            deg_taylor(is)  = m;
            step_size(is)   = t_step;
            maxdA(is)       = alpha;
            ncore(is)       = length(size(v));
            cput(is)        = toc;
            fsize(is,:)     = FSP_qtt_size;
            is              = is+1;


% %           Return to the normal time step, in case the current time
% %           step is too small because of probability sum in the past
%             alpha           = tt_max(-diag(A));
%             nstep           = ceil(alpha*(final_time-t)/sc_factor);
%             t_step          = (final_time-t)/nstep;
        end

    end
    cpu_time    = toc;
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
%-------------------------Report Internal time--------------------------
      iter_qtt = iter;
      T_total_elapsed_qtt = toc(T_total);
      % fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
      fprintf('   Total iterations: %d\n',iter);
      fprintf('   Total time for cme_solver_qtt_unif: %f\n',T_total_elapsed_qtt);
      fprintf('   Total time for AMEn: %f\n',t_amen_total_qtt);


      fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');

      fprintf('   Percentage of runtime for AMEn is %d%%\n',round(100*t_amen_total_qtt/T_total_elapsed_qtt));
      fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
%-----------------------------------------------------------------------
