% %-------------------------Report Internal time for ---------------------
% %===================================CME IN QTT FORMAT WITH MOVING WINDOW
% %========================================SOLVED WITH UNIFORMIZATION+AMEN
%       % fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
%       fprintf('   Total iterations: %d\n',iter_qttmw);
%       fprintf('   Total time for cme_solver_qttmw_unif: %f\n',T_total_elapsed_qttmw);
%       fprintf('   Total time for generating CME matrix: %f\n',t_generator_elapsed_qttmw);
%       fprintf('   Total time for reducing state space: %f\n',t_reduce_elapsed_qttmw);
%       fprintf('   Total time for expanding state space: %f\n',t_expand_elapsed_qttmw);
%       fprintf('   Total time for transforming to QTT-state space: %f\n',t_transform_elapsed_qttmw);
%       fprintf('   Total time for updating probability: %f\n',t_update_elapsed_qttmw);
%
%       % fprintf('   Total time for Uniformization: %f\n',t_unif_elapsed);
%       fprintf('   Total time for AMEn: %f\n',t_amen_total_qttmw);
%
%
%       fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
%       fprintf('   Percentage of runtime for generating CME matrix is %d%%\n',round(100*t_generator_elapsed_qttmw/T_total_elapsed_qttmw));
%       fprintf('   Percentage of runtime for reducing state space is %d%%\n',round(100*t_reduce_elapsed_qttmw/T_total_elapsed_qttmw));
%       fprintf('   Percentage of runtime for expanding state space is %d%%\n',round(100*t_expand_elapsed_qttmw/T_total_elapsed_qttmw));
%       fprintf('   Percentage of runtime for transforming to QTT-state space is %d%%\n',round(100*t_transform_elapsed_qttmw/T_total_elapsed_qttmw));
%       fprintf('   Percentage of runtime for updating probability is %d%%\n',round(100*t_update_elapsed_qttmw/T_total_elapsed_qttmw));
%
%       % fprintf('   Percentage of runtime for Uniformization is %d%%\n',round(100*t_unif_elapsed/T_total_elapsed_qttmw));
%       fprintf('   Percentage of runtime for AMEn is %d%%\n',round(100*t_amen_total_qttmw/T_total_elapsed_qttmw));
%       fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
% %=====================CME IN QTT FORMAT, SOLVED WITH UNIFORMIZATION+AMEN
% %-------------------------Report Internal time--------------------------
%       % fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
%       fprintf('   Total iterations: %d\n',iter_qtt);
%       fprintf('   Total time for cme_solver_qtt_unif: %f\n',T_total_elapsed_qtt);
%       fprintf('   Total time for AMEn: %f\n',t_amen_total_qtt);
%
%
%       fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
%       fprintf('   Percentage of runtime for AMEn is %d%%\n',round(100*t_amen_total_qtt/T_total_elapsed_qtt));
%       fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
% %-----------------------------------------------------------------------
%================================================ANALYSIS OF THE RESULTS
    % for species=1:no_species
    %     fprintf('======================================================\n');
    %     fprintf('Marginal distributions for species %d:\n\n',species);
    %     for method=1:no_methods
    %         if ~isempty(marginal_dist{method,species})
    %             if method==1
    %                 fprintf('SSA trajectory-based frequency: ----------------------\n');
    %             elseif method==2
    %                 fprintf('CME in QTT format using Uniformization + AMEn: -------\n');
    %             elseif method==3
    %                 fprintf('CME in QTT SLIDING format using Uniformization + AMEn: --------\n');
    %             elseif method==4
    %                 fprintf('Adaptive FSP-SSA CME:---------------------------------\n');
    %             disp(marginal_dist{method,species});
    %              end
    %         end
    %      end
    %   end
%========================================PLOT THE MARGINAL DISTRIBUTIONS
    figure(1);clf
    fprintf('======================================================\n');
    for method=1:no_methods
        if (method==1)&&exist('t_ssa','var')
            fprintf('RUN TIME FOR THE SSA:\n');
            disp(t_ssa);
            name_plot               = 'SSA frequency';
            style_plot              = 'bar';
            cme_plot_marginal_distribution(all_marginal_dist_1,all_lb_now_1,all_ub_now_1,style_plot,name_plot,'g');
        elseif (method==2)&& exist('t_qtt','var')
            fprintf('RUN TIME FOR CME IN QTT FORMAT:\n');
            disp(t_qtt)
            name_plot               = 'FSP-QTT';
            style_plot              = 'plot';
            cme_plot_marginal_distribution(all_marginal_dist_2,all_lb_now_2,all_ub_now_2,style_plot,name_plot,'-r',4)
        elseif (method==3)&& exist('t_qtt_sw','var')
            fprintf('RUN TIME FOR CME IN QTT FORMAT WITH SLIDING WINDOW:\n');
            disp(t_qtt_sw)
            name_plot               = 'FSP-QTT with sliding windows';
            style_plot              = 'semilogy';
            cme_plot_marginal_distribution(all_marginal_dist_3,all_lb_now_3,all_ub_now_3,style_plot,name_plot,'--k',4)
      elseif (method==4)&& exist('t_ssa_fsp','var')
          fprintf('RUN TIME FOR THE "ADAPTIVE" FSP-SSA CME:\n');
          disp(t_ssa_fsp)
          name_plot               = 'Adaptive SSA sparse matrix';
          style_plot              = 'scatter';
          cme_plot_marginal_distribution(all_marginal_dist_4,all_lb_now_4,all_ub_now_4,style_plot,name_plot,300,'o','m','m',2);

        end
    end
    cme_plot_marginal_distribution_settings(all_lb_now_3,all_ub_now_3);
% %================================PLOT THE STATE SPACE SIZE AND TIME STEP
%     figure(2);clf
%     for method=1:no_methods
%         if (method==2)&& exist('t_qtt','var')
%             subplot(1,2,1);
%             plot(vec_t_2,vec_ss_size_2,'-k','LineWidth',3);hold on
%             set(gca,'FontSize',20);
%             % l = legend;
%             % l.FontSize  = 20;
%
%             subplot(1,2,2);
%             name_plot               = 'FSP-QTT';
%             plot(vec_t_2,vec_stepsize_2,'-k','LineWidth',3,'DisplayName',name_plot); hold on
%             set(gca,'FontSize',20);
%             l = legend;
%             l.FontSize  = 20;
%         elseif (method==3)&&exist('t_qtt_sw','var')
%             subplot(1,2,1);
%             plot(vec_t_3,vec_ss_size_3,'--r','LineWidth',3);hold on
%             % xlabel('State space size over time');
%             xlabel('Time (s)');
%             ylabel('State space size');
%             set(gca,'FontSize',20);
%             % l = legend;
%             % l.FontSize  = 20;
%
%             subplot(1,2,2);
%             name_plot               = 'FSP-QTT with sliding windows';
%             plot(vec_t_3,vec_stepsize_3,'--r','LineWidth',3,'DisplayName',name_plot);
%             % xlabel('Step-size over time');
%             xlabel('Time (s)');
%             ylabel('Step size');
%             set(gca,'FontSize',20);
%             l = legend;
%             l.FontSize  = 20;
%
%       elseif (method==4)&& exist('t_qtt_psw','var')
%             subplot(1,2,1);
%             plot(vec_t_4,vec_ss_size_4,'--r','LineWidth',3);hold on
%             xlabel('State space size over time');
%             set(gca,'FontSize',20);
%             l = legend;
%             l.FontSize  = 20;
%
%             subplot(1,2,2);
%             name_plot               = 'FSP-QTT with parallel sliding windows';
%             plot(vec_t_4,vec_stepsize_4,'--r','LineWidth',3);hold on
%             xlabel('Step-size over time');
%             set(gca,'FontSize',20);
%         end
%     end
%     %=======================PLOT THE HEATMAP FOR 2-DIMENSIONAL DISTRIBUTIONS
%       if (strcmp(model_name,'Gene_toggle'))
          % figure(3);clf
%
%           for method=1:no_methods
%               if (method==2)&& exist('t_qtt','var')
%                   if exist('t_qtt','var')&&exist('t_qtt','var')
%                       subplot(1,2,1);
%                   end
%                   data_type               = 'multiple-window';
%                   name_plot               = 'FSP-QTT';
%                   cme_plot_joint_distribution(all_w_2,all_lb_now_2,all_ub_now_2,data_type,name_plot); hold on
%             elseif (method==3)&& exist('t_qtt_sw','var')
%                   if exist('t_qtt_sw','var')&&exist('t_qtt_sw','var')
%                       subplot(1,2,2);
%                   end
%                   data_type               = 'multiple-windows';
%                   name_plot               = 'FSP-QTT with sliding windows';
                  % cme_plot_joint_distribution(all_w_3,all_lb_now_3_joint,all_ub_now_3_joint,data_type,name_plot);
%               end
%           end
%       end
% =======================PLOT THE HEATMAP FOR 2-DIMENSIONAL DISTRIBUTIONS
      % figure(4);clf
      %
      % for method=1:no_methods
      %     if (method==3)&& exist('t_qtt_sw','var')
      %         if exist('t_qtt_sw','var')&&exist('t_qtt_psw','var')
      %             subplot(1,2,1);
      %         end
      %         data_type               = 'multiple-window';
      %         name_plot               = 'FSP-QTT with sliding windows';
      %         cme_plot_joint_distribution(all_w_3,all_lb_now_3,all_ub_now_3,data_type,name_plot);
      %    elseif (method==4)&& exist('t_qtt_psw','var')
      %         if exist('t_qtt_sw','var')&&exist('t_qtt_psw','var')
      %             subplot(1,2,2);
      %         end
      %         data_type               = 'multiple-windows';
      %         name_plot               = 'FSP-QTT with parallel sliding windows';
      %         cme_plot_joint_distribution(all_w_4,all_lb_now_4,all_ub_now_4,data_type,name_plot);
      %     end
      % end
%====================================PLOT THE TRAJECTORY OF STATE SPACES
      % figure(5);clf
      %
      % if (method==4)&& exist('t_qtt_psw','var')
      %     cme_plot_ss_trajectory(vec_t_4,vec_windows_lb_4,vec_windows_l2size_4)
      % end
% ====================================PLOT THE TRAJECTORY OF STATE SPACES
%       figure(6);clf
%
% %       if (method==3)&& exist('t_qtt_sw','var')
          cme_plot_ss_trajectory(vec_t_3,vec_windows_lb_3,vec_windows_l2size_3)
          % cme_plot_mini_ss_trajectory(all_w_3,vec_t_3,vec_windows_lb_3,vec_windows_l2size_3)
%       end
% [w_3,lb_w_3,ub_w_3,marginal_dist_3,lb_now,ub_now,vec_t_3,vec_stepsize_3,vec_ss_size_3] = cme_solver_qttmw_unif(propen_func,propen_func_partial);
