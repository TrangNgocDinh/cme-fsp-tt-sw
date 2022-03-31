    global FSP_qtt_size_initial FSP_qtt_size

    no_methods              = 100;

% %==========================================================SSA FREQUENCY
%     tic
%     [marginal_dist_1,lb,ub] = cme_solver_ssa(propen_func);
%     all_marginal_dist_1     = cell(1);
%     all_marginal_dist_1{1}  = marginal_dist_1;
%     all_lb_now_1            = lb;
%     all_ub_now_1            = ub;
%     t_ssa                   = toc;
% %=========================================ONE-STEP DENSE TRADITIONAL CME
%     tic
% %   Build the CME matrix
%     [A,state_space_list_2]  = cme_generator_traditional_dense(propen_func);
% %   Build the initial probability vector
%     p_0                     = cme_initial_traditional_dense(state_space_list_2,initial_state);
% %   Solve the CME
%     p_out_2                 = cme_solver_expm(final_time,A,p_0);
%     marginal_dist_2         = cme_marginal_traditional(p_out_2,state_space_list_2);
%     all_marginal_dist_2     = cell(1);
%     all_marginal_dist_2{1}  = marginal_dist_2;
%     all_lb_now_2            = vec_ss_minsize;
%     all_ub_now_2            = vec_ss_maxsize;
%     t_one_step_dense        = toc;
% %========================================ONE-STEP SPARSE TRADITIONAL CME
%     tic
% %   Build the CME matrix
%     [A,state_space_list_3]  = cme_generator_traditional_sparse(propen_func);
% %   Build the initial probability vector
%     p_0                     = cme_initial_traditional_sparse(state_space_list_3,initial_state);
% %   Solve the CME
%     p_out_3                 = cme_solver_expokit(final_time,A,p_0,tol_expokit,deg_expokit);
%     marginal_dist_3         = cme_marginal_traditional(p_out_3,state_space_list_3);
%     all_marginal_dist_3     = cell(1);
%     all_marginal_dist_3{1}  = marginal_dist_2;
%     all_lb_now_3            = vec_ss_minsize;
%     all_ub_now_3            = vec_ss_maxsize;
%     t_one_step_sparse       = toc;
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
% %=====================CME IN QTT FORMAT, SOLVED WITH UNIFORMIZATION+AMEN
%     FSP_qtt_size            = FSP_qtt_size_initial;
%     tic
% %   Build the CME matrix
%     A                       = cme_generator_qtt(propen_func_partial);
% %   Build the initial probability vector
%     p_0                     = cme_initial_qtt;
% %   Solve the CME with Uniformization and AMEn
%     [~,marginal_dist_5,vec_t_5,vec_stepsize_5,vec_ss_size_5]    = cme_solver_qtt_unif(p_0,propen_func_partial);
%     all_marginal_dist_5{1}  = marginal_dist_5;
%     all_lb_now_5            = zeros(1,no_species);
%     all_ub_now_5            = 2.^FSP_qtt_size-1;
%     t_qtt                   = toc;
% %=================================CME IN QTT FORMAT, SOLVED WITH EXPOKIT
%     FSP_qtt_size            = FSP_qtt_size_initial;
%     tic
% %   Build the initial probability vector
%     p_0                     = cme_initial_qtt;
% %   Solve the CME with EXPOKIT
%     [~,marginal_dist_6]     = cme_solver_qtt_expokit(p_0,propen_func_partial);
%     all_marginal_dist_6{1}  = marginal_dist_6;
%     all_lb_now_6            = zeros(1,no_species);
%     all_ub_now_6            = 2.^FSP_qtt_size-1;
%     t_qtt_expokit           = toc;
%===================================CME IN QTT FORMAT WITH MOVING WINDOW
%========================================SOLVED WITH UNIFORMIZATION+AMEN
    % tic
    lb              = initial_state;
    l2size          = 4*ones(1,no_species);
    A               = cme_generator_qttmw(lb,l2size,propen_func_partial)
    % FSP_qtt_size            = FSP_qtt_size_initial;
    % [w_7,lb_w_7,ub_w_7,marginal_dist_7,lb_now,ub_now,vec_t_7,vec_stepsize_7,vec_ss_size_7] = cme_solver_qttmw_unif(propen_func,propen_func_partial);
    % all_marginal_dist_7     = cell(1);
    % all_marginal_dist_7{1}  = marginal_dist_7;
    % all_w_7                 = cell(1);
    % all_w_7{1}              = w_7;
    % all_lb_now_7            = lb_now;
    % all_ub_now_7            = ub_now;
    % t_qtt_sw                = toc;
%==========================CME IN QTT FORMAT WITH PARALLEL MOVING WINDOW
%========================================SOLVED WITH UNIFORMIZATION+AMEN
    % tic
    % [all_w_8,all_lb_w_now_8,all_ub_w_now_8,all_marginal_dist_8,all_lb_now_8,all_ub_now_8,vec_t_8,vec_stepsize_8,vec_ss_size_8,vec_windows_lb_8,vec_windows_l2size_8] = cme_solver_qttpsw_unif(propen_func,propen_func_partial);
    % t_qtt_psw               = toc;



% %================================================
% %================================================
% %================================================
% %================================================
% %================================================
%     filename    = ['Data_' model_name '_Thanh.mat'];
%     % save(filename);
%     load(filename)
% %================================================
% %================================================
% %================================================
% %================================================
% %================================================



% %================================================ANALYSIS OF THE RESULTS
%     for species=1:no_species
%         fprintf('======================================================\n');
%         fprintf('Marginal distributions for species %d:\n\n',species);
%         for method=1:no_methods
%             if ~isempty(marginal_dist{method,species})
%                 if method==1
%                     fprintf('SSA trajectory-based frequency: ----------------------\n');
%                 elseif method==2
%                     fprintf('One-step dense traditional CME: ----------------------\n');
%                 elseif method==3
%                     fprintf('One-step sparse traditional CME: ---------------------\n');
%                 elseif method==4
%                     fprintf('Adaptive FSP-SSA CME:---------------------------------\n');
%                 elseif method==5
%                     fprintf('CME in QTT format using Uniformization + AMEn: -------\n');
%                 elseif method==6
%                     fprintf('CME in QTT format using TT-Expokit: ------------------\n');
%                 elseif method==7
%                     fprintf('CME in TT format using Uniformization + AMEn: --------\n');
%                 end
%                 disp(marginal_dist{method,species});
%             end
%         end
%     end
%========================================PLOT THE MARGINAL DISTRIBUTIONS
    % figure(1);clf
    % fprintf('======================================================\n');
    % for method=1:no_methods
    %     if (method==100)&&exist('t_ssa','var')
    %         fprintf('RUN TIME FOR THE SSA:\n');
    %         disp(t_ssa);
    %         name_plot               = 'SSA frequency';
    %         style_plot              = 'bar';
    %         cme_plot_marginal_distribution(all_marginal_dist_1,all_lb_now_1,all_ub_now_1,style_plot,name_plot,'g');
    %     elseif (method==2)&&exist('t_one_step_dense','var')
    %         fprintf('RUN TIME FOR THE ONE-STEP DENSE CME:\n');
    %         disp(t_one_step_dense);
    %         name_plot               = 'One-step dense matrix';
    %         style_plot              = 'plot';
    %         cme_plot_marginal_distribution(all_marginal_dist_2,all_lb_now_2,all_ub_now_2,style_plot,name_plot,'-k',4);
    %     elseif (method==3)&& exist('t_one_step_sparse','var')
    %         fprintf('RUN TIME FOR THE ONE-STEP SPARSE CME:\n');
    %         disp(t_one_step_sparse)
    %         name_plot               = 'One-step sparse matrix';
    %         style_plot              = 'scatter';
    %         cme_plot_marginal_distribution(all_marginal_dist_3,all_lb_now_3,all_ub_now_3,style_plot,name_plot,800,'s','b','k',3);
    %     elseif (method==4)&& exist('t_ssa_fsp','var')
    %         fprintf('RUN TIME FOR THE "ADAPTIVE" FSP-SSA CME:\n');
    %         disp(t_ssa_fsp)
    %         name_plot               = 'Adaptive SSA sparse matrix';
    %         style_plot              = 'scatter';
    %         cme_plot_marginal_distribution(all_marginal_dist_4,all_lb_now_4,all_ub_now_4,style_plot,name_plot,300,'o','m','m',2);
    %     elseif (method==5)&& exist('t_qtt','var')
    %         fprintf('RUN TIME FOR CME IN QTT FORMAT:\n');
    %         disp(t_qtt)
    %         name_plot               = 'FSP-QTT';
    %         style_plot              = 'plot';
    %         cme_plot_marginal_distribution(all_marginal_dist_5,all_lb_now_5,all_ub_now_5,style_plot,name_plot,'-r',4)
    %     elseif (method==6)&& exist('t_qtt_expokit','var')
    %         fprintf('RUN TIME FOR CME IN QTT FORMAT USING EXPOKIT:\n');
    %         disp(t_qtt_expokit)
    %         name_plot               = 'QTT-Expokit';
    %         style_plot              = 'scatter';
    %         cme_plot_marginal_distribution(all_marginal_dist_6,all_lb_now_6,all_ub_now_6,style_plot,name_plot,120,'s','k','k',2)
    %     elseif (method==7)&& exist('t_qtt_sw','var')
    %         fprintf('RUN TIME FOR CME IN QTT FORMAT WITH SLIDING WINDOW:\n');
    %         disp(t_qtt_sw)
    %         name_plot               = 'FSP-QTT with sliding windows';
    %         style_plot              = 'semilogy';
    %         cme_plot_marginal_distribution(all_marginal_dist_7,all_lb_now_7,all_ub_now_7,style_plot,name_plot,'-k',4)
    %     elseif (method==8)&& exist('t_qtt_psw','var')
    %         fprintf('RUN TIME FOR CME IN QTT FORMAT WITH PARALLEL SLIDING WINDOWS:\n');
    %         disp(t_qtt_psw)
    %         name_plot               = 'FSP-QTT with parallel sliding windows';
    %         style_plot              = 'semilogy';
    %         cme_plot_marginal_distribution(all_marginal_dist_8,all_lb_now_8,all_ub_now_8,style_plot,name_plot,'--r',4)
    %     end
    % end
    % cme_plot_marginal_distribution_settings(all_lb_now_7,all_ub_now_7);
% %================================PLOT THE STATE SPACE SIZE AND TIME STEP
%     figure(2);clf
%     for method=1:no_methods
%         if (method==5)&& exist('t_qtt','var')
%             name_plot               = 'FSP-QTT';
%             subplot(1,2,1);
%             plot(vec_t_5,vec_ss_size_5,'-k','LineWidth',3,'DisplayName',name_plot);hold on
%             set(gca,'FontSize',20);
%             l = legend;
%             l.FontSize  = 20;
%
%             subplot(1,2,2);
%             plot(vec_t_5,vec_stepsize_5,'-k','LineWidth',3);hold on
%             set(gca,'FontSize',20);
%         elseif (method==7)&&exist('t_qtt_sw','var')
%             name_plot               = 'FSP-QTT with sliding windows';
%             subplot(1,2,1);
%             plot(vec_t_7,vec_ss_size_7,'-b','LineWidth',3,'DisplayName',name_plot);hold on
%             xlabel('State space size over time');
%             set(gca,'FontSize',20);
%             l = legend;
%             l.FontSize  = 20;
%
%             subplot(1,2,2);
%             plot(vec_t_7,vec_stepsize_7,'-b','LineWidth',3);hold on
%             xlabel('Step-size over time');
%             set(gca,'FontSize',20);
%         elseif (method==8)&& exist('t_qtt_psw','var')
%             name_plot               = 'FSP-QTT with parallel sliding windows';
%             subplot(1,2,1);
%             plot(vec_t_8,vec_ss_size_8,'--r','LineWidth',3,'DisplayName',name_plot);hold on
%             xlabel('State space size over time');
%             set(gca,'FontSize',20);
%             l = legend;
%             l.FontSize  = 20;
%
%             subplot(1,2,2);
%             plot(vec_t_8,vec_stepsize_8,'--r','LineWidth',3);hold on
%             xlabel('Step-size over time');
%             set(gca,'FontSize',20);
%         end
%     end
% %=======================PLOT THE HEATMAP FOR 2-DIMENSIONAL DISTRIBUTIONS
%     figure(3);clf
%
%     for method=1:no_methods
%         if (method==7)&& exist('t_qtt_sw','var')
%             if exist('t_qtt_sw','var')&&exist('t_qtt_psw','var')
%                 subplot(1,2,1);
%             end
%             data_type               = 'multiple-window';
%             name_plot               = 'FSP-QTT with sliding windows';
%             cme_plot_joint_distribution(all_w_7,all_lb_now_7,all_ub_now_7,data_type,name_plot);
%         elseif (method==8)&& exist('t_qtt_psw','var')
%             if exist('t_qtt_sw','var')&&exist('t_qtt_psw','var')
%                 subplot(1,2,2);
%             end
%             data_type               = 'multiple-windows';
%             name_plot               = 'FSP-QTT with parallel sliding windows';
%             cme_plot_joint_distribution(all_w_8,all_lb_now_8,all_ub_now_8,data_type,name_plot);
%         end
%     end
% %====================================PLOT THE TRAJECTORY OF STATE SPACES
%     figure(4);clf
%
%     if (method==8)&& exist('t_qtt_psw','var')
%         cme_plot_ss_trajectory(vec_t_8,vec_windows_lb_8,vec_windows_l2size_8)
%     end
