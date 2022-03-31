function [all_w,all_lb_w_now,all_ub_w_now,all_marginal,all_lb_now,all_ub_now,vec_t,vec_stepsize,vec_ss_size,vec_windows_lb,vec_windows_l2size] = cme_solver_qttpsw_unif(propen_func,propen_func_partial)
    global tol_fsp tol_unif tol_fsp_probsum tol_amen sc_factor final_time no_species
    global stoich_mat_total initial_state

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
%   Initial window
    lb_now          = initial_state;
    l2size_now      = 2*ones(1,no_species);
    ub_now          = lb_now+2.^l2size_now-1;
%   Initial probability vector in QTT format
    v               = cme_initial_qttmw(lb_now,l2size_now);



%   Data for the initial window
    no_windows      = 1;
    all_lb_now      = [lb_now];
    all_l2size_now  = [l2size_now];
    all_v           = cell(1,1);
    all_v{1}        = v;



%   CME matrix in QTT format
    A               = cme_generator_qttmw(lb_now,l2size_now,propen_func_partial);
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
    count           = 0;
    vec_t           = zeros(1,1);
    vec_stepsize    = zeros(1,1);
    vec_ss_size     = zeros(1,1);
    acc_error       = zeros(1,no_species);

    logic_redo      = 0;

    tic
    while (t<final_time)
        fprintf('==================================================\n');
        fprintf('At time %f; final time %f\n',t,final_time);
        iter                = iter+1;
%       Update the windows and probability vectors
        [no_windows_new,all_v_new,all_lb_new,all_l2size_new] = cme_update_ss_qttpsw(no_windows,all_v,all_lb_now,all_l2size_now,logic_redo,t_step,propen_func);
%       Re-generate the CME matrices
        fprintf('Expanding A-qtt...\n');
        all_A           = cell(no_windows_new,1);
        for i=1:no_windows_new
            lb_now      = all_lb_new(i,:);
            l2size_now  = all_l2size_new(i,:);
            A           = cme_generator_qttmw(lb_now,l2size_now,propen_func_partial);
            all_A{i}    = A;
        end
%       Decide the time step
        t_step          = Inf;
        all_nstep       = zeros(no_windows_new,1);
        for i=1:no_windows_new
            A           = all_A{i};
            alpha       = tt_max(-diag(A));
            nstep       = ceil(alpha*(final_time-t)/sc_factor);
            all_nstep(i)= nstep;
            t_step      = min(t_step,(final_time-t)/nstep);
        end
        fprintf('Time step = %f...\n',t_step);
%       Find the next probability tensor by using Uniformization & AMEn
        fprintf('Solving for the next probability tensor...\n');
        all_w           = cell(no_windows_new,1);
        all_W           = cell(no_windows_new,1);
        for i=1:no_windows_new
            nstep       = all_nstep(i);
            v           = all_v_new{i};
            A           = all_A{i};
            [w,W,d,m]   = unif(A,v,t_step,nstep);
            all_w{i}    = w;
            all_W{i}    = W;
        end
%       Sanity test: if the probability sum is too small, then
%       reduce the time step and redo
        prob_sum    = 0;
        for i=1:no_windows_new
            w           = all_w{i};
            prob_sum    = prob_sum+sum(w);
        end
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
        vec_ss_size(count)  = 0;
        for i=1:no_windows_new
            l2size_new          = all_l2size_new(i,:);
            vec_ss_size(count)  = vec_ss_size(count)+2^(sum(l2size_now));
        end
%       Record the windows
        vec_windows_lb{count}       = all_lb_new;
        vec_windows_l2size{count}   = all_l2size_new;
%       Advance to the next step
        all_v           = all_w;
        no_windows      = no_windows_new;
        all_lb_now      = all_lb_new;
        all_l2size_now  = all_l2size_new;
        t               = t + t_step;
%       Check for happy breakdown of all windows
        logic_happybreakdown    = 1;
        for i=1:no_windows_new
            W       = all_W{i};
            rW      = rank(W);
            if (rW(d+1)==1)
%               If the rank of the 'degree' core drops to 1, we have
%               happy breakdown for this window
            else
                logic_happybreakdown    = 0;
            end
        end
        if logic_happybreakdown==1
            disp('Happy breakdown');
            break
        end
    end
    cpu_time    = toc;
%   Generate the marginal distributions
    for i=1:no_windows_new
        v                           = all_v{i};
        lb_now                      = all_lb_now(i,:);
        l2size_now                  = all_l2size_now(i,:);
        ub_now                      = lb_now+2.^l2size_now-1;
        w                           = reshape(v,2.^l2size_now);

        all_lb_w_now(i,:)           = lb_now;
        all_ub_w_now(i,:)           = ub_now;

        [lb_now,ub_now,marginal]    = cme_marginal_qttmw(lb_now,ub_now,w);
        all_w{i}                    = w;
        all_marginal{i}             = marginal;
        all_lb_now(i,:)             = lb_now;
        all_ub_now(i,:)             = ub_now;
    end
