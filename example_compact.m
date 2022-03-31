    global FSP_qtt_size_initial FSP_qtt_size

    no_methods              = 4;

% % ==========================================================SSA FREQUENCY
% %   Uncomment  below if run the model the 1st time
%   % And Comment them the below if..end after the 1st run
%     tic
%     [marginal_dist_1,lb,ub] = cme_solver_ssa(propen_func);
%     all_marginal_dist_1     = cell(1);
%     all_marginal_dist_1{1}  = marginal_dist_1;
%     all_lb_now_1            = lb;
%     all_ub_now_1            = ub;
%     t_ssa                   = toc;
%
%     if (strcmp(model_name,'P53'))
%           save data_SSA_p53.mat;
%     elseif (strcmp(model_name,'Gene_toggle'))
%           save data_SSA_gene_toggle.mat;
%     elseif (strcmp(model_name,'Michaelis-Menten'))
%           save data_SSA_MM.mat;
%     else
%           save data_SSA_goutsias.mat;
%     end
% % %=======================================================================
% % %=====================CME IN QTT FORMAT, SOLVED WITH UNIFORMIZATION+AMEN
%     FSP_qtt_size            = FSP_qtt_size_initial;
%     tic
% %   Build the CME matrix
%     A                       = cme_generator_qtt(propen_func_partial);
% %   Build the initial probability vector
%     p_0                     = cme_initial_qtt;
% %   Solve the CME with Uniformization and AMEn
%     [w_2,marginal_dist_2,vec_t_2,vec_stepsize_2,vec_ss_size_2,T_total_elapsed_qtt,t_amen_total_qtt,iter_qtt]    = cme_solver_qtt_unif(p_0,propen_func_partial);
%     all_w_2                 = cell(1);
%     all_w_2{1}              = w_2;
%     all_marginal_dist_2{1}  = marginal_dist_2;
%     all_lb_now_2            = zeros(1,no_species);
%     all_ub_now_2            = 2.^FSP_qtt_size-1;
%     t_qtt                   = toc;
% %===================================CME IN QTT FORMAT WITH MOVING WINDOW
% %========================================SOLVED WITH UNIFORMIZATION+AMEN
    tic
    FSP_qtt_size            = FSP_qtt_size_initial;
    [all_w_3,count_3,w_3,lb_w_3,ub_w_3,marginal_dist_3,lb_now,ub_now,vec_t_3,vec_stepsize_3,vec_ss_size_3, vec_windows_lb_3,vec_windows_l2size_3, T_total_elapsed_qttmw, t_reduce_elapsed_qttmw, t_expand_elapsed_qttmw, t_update_elapsed_qttmw,t_generator_elapsed_qttmw, t_transform_elapsed_qttmw,t_amen_total_qttmw,iter_qttmw] = cme_solver_qttmw_unif(propen_func,propen_func_partial);
    all_marginal_dist_3     = cell(1);
    all_marginal_dist_3{1}  = marginal_dist_3;

    % all_w_3                 = cell(1);
    % all_w_3{1}              = w_3;

    all_lb_now_3            = lb_now;
    all_ub_now_3            = ub_now;

    all_lb_now_3_joint      = lb_w_3;
    all_ub_now_3_joint      = ub_w_3;

    t_qtt_sw                = toc;
% %===================================================ADAPTIVE SSA-FSP CME
%     tic
% %   Solve the CME
%     [p_out_4,state_space_list_4,vec_ss_minsize,vec_ss_maxsize]  = cme_solver_SSA_expokit(propen_func);
%     marginal_dist_4         = cme_marginal_traditional(p_out_4,state_space_list_4);
%     all_marginal_dist_4     = cell(1);
%     all_marginal_dist_4{1}  = marginal_dist_4;
%     all_lb_now_4            = zeros(1,no_species);
%     all_ub_now_4            = vec_ss_maxsize;
%     t_ssa_fsp               = toc;

%====================================
if (strcmp(model_name,'P53'))
      save data_p53.mat;
elseif (strcmp(model_name,'Gene_toggle'))
      save data_gene_toggle.mat;
elseif (strcmp(model_name,'Michaelis-Menten'))
      save data_MM.mat;
else
      save data_goutsias.mat;
end
