function [all_w,count,w,lb_w,ub_w,marginal,lb_now,ub_now,vec_t,vec_stepsize,vec_ss_size, vec_windows_lb, vec_windows_l2size, T_total_elapsed_qttmw, t_reduce_elapsed_qttmw, t_expand_elapsed_qttmw, t_update_elapsed_qttmw,t_generator_elapsed_qttmw, t_transform_elapsed_qttmw,t_amen_total_qttmw,iter_qttmw] = cme_solver_qttmw_unif(propen_func,propen_func_partial)

    global tol_fsp tol_unif tol_fsp_probsum tol_amen sc_factor final_time no_species
    global stoich_mat_total initial_state

    all_w = {};
%-------------------------Internal statistics---------------------------
      % global T_total_elapsed_qttmw t_reduce_elapsed_qttmw
      % global t_expand_elapsed_qttmw t_update_elapsed_qttmw
      % global t_generator_elapsed_qttmw t_transform_elapsed_qttmw
      % global t_amen_total_qttmw

      T_total = tic;
      t_reduce_elapsed_qttmw        = 0;
      t_expand_elapsed_qttmw         = 0;
      t_update_elapsed_qttmw         = 0;
      t_generator_elapsed_qttmw     = 0;
      t_transform_elapsed_qttmw      = 0;
      % t_unif_elapsed           = 0;
      t_amen_total_qttmw             = 0;
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

%   Initial window using lower_bound, upper_bound
    % [~,~,lb_now,ub_now] = ssa(initial_state,final_time);
    % l2size_now = ceil(log2(ub_now-lb_now+1));
    % ub_now = lb_now + 2.^l2size_now;
    lb_now          = initial_state;
    l2size_now      = 2*ones(1,no_species);
    ub_now          = lb_now+2.^l2size_now-1;
%   Initial probability vector in QTT format
    v               = cme_initial_qttmw(lb_now,l2size_now);
%   CME matrix in QTT format
    t_generator     = tic;

    A               = cme_generator_qttmw(lb_now,l2size_now,propen_func_partial);

    t_generator_elapsed_qttmw      = t_generator_elapsed_qttmw+ toc(t_generator);
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
%   Prepare vectors of step-size and state space size
    count           = 1; %new
    vec_t           = zeros(1,1);
    vec_stepsize    = zeros(1,1);
    vec_ss_size     = zeros(1,1);

%New: 06/09/20
%    all_w:
    vec_t(count)            = t;
    vec_stepsize(count)     = t_step;
    vec_ss_size(count)      = 2^(sum(l2size_now));
    all_w{count}            = tt_reshape(v,2.^l2size_now);
%       Record the windows
    vec_windows_lb{count}       = lb_now;
    vec_windows_l2size{count}   = l2size_now;

%==============end new

%    bmass_old       = zeros(1,no_species);
    acc_error       = zeros(1,no_species);

    logic_redo      = 0;

    tic
    while (t<final_time)
        fprintf('==================================================\n');
        % fprintf('At time %f; final time %f\n',t,final_time);
        fprintf('For iteration %d\n', iter);
        iter                = iter+1;
%       Keep parameters of old window, v and w
        fprintf('Old state space...\n');
        l2size_old          = ceil(log2(ub_now-lb_now+1));
        lb_old              = lb_now;
        ub_old              = lb_old + 2.^l2size_old-1;
        disp(lb_now)
        disp(ub_now)
%       Reduce the state space
        if logic_redo==0
            fprintf('Reducing the state space...\n');

            t_reduce            = tic;

            mode_now            = ub_now - lb_now + 1;
            v_physic            = tt_reshape(v,mode_now);

            [lb_tmp,ub_tmp]     = cme_update_ss_qttmw_reduce(lb_now,ub_now,v_physic);

            t_reduce_elapsed_qttmw    = t_reduce_elapsed_qttmw + toc(t_reduce);

            disp(lb_tmp)
            disp(ub_tmp)
        else
            lb_tmp              = lb_old;
            ub_tmp              = ub_old;
        end

%       Expand the state space using SSA
        fprintf('Expanding the state space...\n');

        t_expand            = tic;

        [lb_new,ub_new]     = cme_update_ss_qttmw_SSA(lb_tmp,ub_tmp,t_step,propen_func);

        t_expand_elapsed_qttmw    = t_expand_elapsed_qttmw + toc(t_expand);

        l2size_now          = ceil(log2(ub_new-lb_new+1));
        l2size_now          = max(l2size_now,2);
        disp(lb_new)
        disp(ub_new)
%       Transform the bounds into QTT format
        fprintf('Transforming the state space into QTT format...\n');

        t_transform         = tic;

        for species=1:no_species
%           Calculate how many states will be redundant because of QTT storage
            excess  = 2.^l2size_now(species)-ub_new(species)+lb_new(species)-1;
%           Distribute the excess states evenly on both sides of the bounds
            lb_now(species) = max(0,lb_new(species)-floor(excess/2));
            ub_now(species) = lb_now(species)+2^l2size_now(species)-1;
        end

        t_transform_elapsed_qttmw  = t_transform_elapsed_qttmw + toc(t_transform);

        disp(lb_now)
        disp(ub_now)
%       Update the probability vector
        fprintf('Updating the probability vector...\n')

        t_update            = tic;

        v                   = cme_update_prob_qttmw(v,lb_old,l2size_old,lb_now,l2size_now);

        t_update_elapsed_qttmw    = t_update_elapsed_qttmw + toc(t_update);
%       Re-generate the matrix and the one-step linear system
        fprintf('Expanding A-qtt...\n');

        t_generator     = tic;

        A               = cme_generator_qttmw(lb_now,l2size_now,propen_func_partial);

        t_generator_elapsed_qttmw      = t_generator_elapsed_qttmw+ toc(t_generator);
%       Decide the time step
        alpha           = tt_max(-diag(A));
        nstep           = ceil(alpha*(final_time-t)/sc_factor);
        t_step          = (final_time-t)/nstep;
        % t_step          = min(t_step,(final_time-t)/nstep);
        % fprintf('Time step = %f...\n',t_step);
%       Find the next probability tensor by using Uniformization & AMEn
        fprintf('Solving for the next probability tensor...\n');

        % t_unif          = tic;

        [w,W,d,m,t_amen_elapsed]           = unif(A,v,t_step,nstep);


        % t_unif_elapsed  = t_unif_elapsed + toc(t_unif);
        fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
        fprintf('At time %f; final time %f\n',t,final_time);
        fprintf('   Time step = %f\n',t_step);
        fprintf('   Theta = %f ; m = %d; time for AMEn is %f\n', alpha, m, t_amen_elapsed);
        t_amen_total_qttmw    = t_amen_total_qttmw   + t_amen_elapsed;
        fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
%       Sanity test: if the probability sum is too small, then
%       reduce the time step and redo
        prob_sum    = sum(w);
        fprintf('Sum of next probability tensor = %f...\n',prob_sum);
        tol_probsum = 1-tol_fsp_probsum*((t+t_step)/final_time);
        if prob_sum<tol_probsum
            disp('Probability sum is too low, redo this step with larger SSA expansion...');
            logic_redo  = 1;
            continue;
        else
            logic_redo  = 0;
        end
%       Record the step-size and state space size
        count               = count+1;
        vec_t(count)        = t;
        vec_stepsize(count) = t_step;
        vec_ss_size(count)  = 2^(sum(l2size_now));
%       Record the windows
        vec_windows_lb{count}       = lb_new;
        vec_windows_l2size{count}   = l2size_now;

        % all_w{count}                = reshape(w,2.^l2size_now);
        all_w{count}                = reshape(w,2.^l2size_now);
%       Advance to the next step
        v       = w;
        t       = t + t_step;
%       Initial guess for the next AMEn solve
        W0      = W;
        rW      = rank(W);
        if (rW(d+1)==1)
%           If the rank of the 'degree' core drops to 1, we have happy
%           breakdown...
            disp('Happy breakdown');
            break
        end
%       Record some statistics
        timetrack(is)   = toc;
        difftrack(is)   = difftracker(is);
        memp(is)        = mem(v);
        memW(is)        = mem(W);
        % fsperr(is)      = max(bmass);
        W_max_rank(is)  = max(rank(W));
        p_max_rank(is)  = max(rank(v));
        deg_taylor(is)  = m;
        step_size(is)   = t_step;
        maxdA(is)       = alpha;
        ncore(is)       = length(size(v));
        cput(is)        = toc;
        is              = is+1;
    end
    cpu_time    = toc;
%   Generate the marginal distributions
    w                           = reshape(v,2.^l2size_now);
    lb_w                        = lb_now;
    ub_w                        = ub_now;
% %   Record probability
% %================ t = 1
%       t1          = vec_t_plot(1);
%       t2          = vec_t_plot(2);
%       t3          = vec_t_plot(3);
%       pos1         = find(vec_t>=t1,1);
%       pos2         = find(vec_t>=t2,1);
%       pos3         = find(vec_t>=t3,1);
%       if ~isempty(pos1)
%           vec_w{1}      = w;
%       end
%
%       if ~isempty(pos2)
%           vec_w{2}      = w;
%       end
%
%       if ~isempty(pos3)
%           vec_w{3}      = w;
%       end


%   Record marginal probability

    [lb_now,ub_now,marginal]    = cme_marginal_qttmw(lb_now,ub_now,w);

%-------------------------Report Internal time--------------------------
      iter_qttmw = iter;
      T_total_elapsed_qttmw = toc(T_total);
      % fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
      fprintf('   Total iterations: %d\n',iter);
      fprintf('   Total time for cme_solver_qttmw_unif: %f\n',T_total_elapsed_qttmw);
      fprintf('   Total time for generating CME matrix: %f\n',t_generator_elapsed_qttmw);
      fprintf('   Total time for reducing state space: %f\n',t_reduce_elapsed_qttmw);
      fprintf('   Total time for expanding state space: %f\n',t_expand_elapsed_qttmw);
      fprintf('   Total time for transforming to QTT-state space: %f\n',t_transform_elapsed_qttmw);
      fprintf('   Total time for updating probability: %f\n',t_update_elapsed_qttmw);

      % fprintf('   Total time for Uniformization: %f\n',t_unif_elapsed);
      fprintf('   Total time for AMEn: %f\n',t_amen_total_qttmw);


      fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
      fprintf('   Percentage of runtime for generating CME matrix is %d%%\n',round(100*t_generator_elapsed_qttmw/T_total_elapsed_qttmw));
      fprintf('   Percentage of runtime for reducing state space is %d%%\n',round(100*t_reduce_elapsed_qttmw/T_total_elapsed_qttmw));
      fprintf('   Percentage of runtime for expanding state space is %d%%\n',round(100*t_expand_elapsed_qttmw/T_total_elapsed_qttmw));
      fprintf('   Percentage of runtime for transforming to QTT-state space is %d%%\n',round(100*t_transform_elapsed_qttmw/T_total_elapsed_qttmw));
      fprintf('   Percentage of runtime for updating probability is %d%%\n',round(100*t_update_elapsed_qttmw/T_total_elapsed_qttmw));

      % fprintf('   Percentage of runtime for Uniformization is %d%%\n',round(100*t_unif_elapsed/T_total_elapsed_qttmw));
      fprintf('   Percentage of runtime for AMEn is %d%%\n',round(100*t_amen_total_qttmw/T_total_elapsed_qttmw));
      fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
%-----------------------------------------------------------------------
