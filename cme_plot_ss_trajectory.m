function cme_plot_ss_trajectory(vec_t,vec_windows_lb,vec_windows_l2size)
    global species_name
    % vec_t_plot  = [0 1 2 5 10 20];
    vec_t_plot  = [0 1 2 5];
% %=========================== Old versions that work
% for i=1:length(vec_t_plot)
%     t           = vec_t_plot(i);
%     pos         = find(vec_t>=t,1);
%     if ~isempty(pos)
%         all_lb      = vec_windows_lb{pos};
%         all_l2size  = vec_windows_l2size{pos};
%         all_ub      = all_lb+2.^all_l2size-1;
%         no_windows  = size(all_lb,1);
%         for window=1:no_windows
%             lb  = all_lb(window,:);
%             ub  = all_ub(window,:);
%             p1  = plot([lb(1) lb(1)],[lb(2) ub(2)],'LineStyle','-','LineWidth',3);hold on
%             if window>=2
%                 p1.Color    = col;
%             else
%                 col = p1.Color;
%             end
%             plot([ub(1) ub(1)],[lb(2) ub(2)],'LineStyle','-','Color',col,'LineWidth',3);hold on
%             plot([lb(1) ub(1)],[lb(2) lb(2)],'LineStyle','-','Color',col,'LineWidth',3);hold on
%             plot([lb(1) ub(1)],[ub(2) ub(2)],'LineStyle','-','Color',col,'LineWidth',3);hold on
%             if (i == 3)
%               text(ub(1)+2,ub(2)+2,['t = ' num2str(t)],'Color',col,'FontSize',20);
%             else
%               text(ub(1)-2,ub(2)+35,['t = ' num2str(t)],'Color',col,'FontSize',20);
%             end
%
%         end
%     end
% end

%============================= MODIFICATION STARTS HERE!
%================ t = 0
    for i=1:1
        t           = vec_t_plot(i);
        pos         = find(vec_t>=t,1);
        if ~isempty(pos)
            all_lb      = vec_windows_lb{pos};
            all_l2size  = vec_windows_l2size{pos};
            all_ub      = all_lb+2.^all_l2size-1;
            no_windows  = size(all_lb,1);
            for window=1:no_windows
                lb  = all_lb(window,:);
                ub  = all_ub(window,:);
                color = [0, 0.4470, 0.7410];
                p1 = rectangle('Position',[lb(1) lb(2) ub(1)-lb(1) ub(2)-lb(2)],'EdgeColor',color,'LineWidth',3);
                if (i == 3)
                  text(ub(1)+2,ub(2)+2,['t = ' num2str(t)],'Color',color,'FontSize',20);
                else
                  text(ub(1),ub(2)+20,['t = ' num2str(t)],'Color',color,'FontSize',20);
                end

            end
        end
    end
%============================ t = 1
for i=2:2
    t           = vec_t_plot(i);
    pos         = find(vec_t>=t,1);
    if ~isempty(pos)
        all_lb      = vec_windows_lb{pos};
        all_l2size  = vec_windows_l2size{pos};
        all_ub      = all_lb+2.^all_l2size-1;
        no_windows  = size(all_lb,1);
        for window=1:no_windows
           lb  = all_lb(window,:);
           ub  = all_ub(window,:);
           color = [0, 0.5, 0];
           p1 = rectangle('Position',[lb(1) lb(2) ub(1)-lb(1) ub(2)-lb(2)],'LineStyle','--','EdgeColor',color,'LineWidth',3);
           if (i == 3)
              text(ub(1)+2,ub(2)+2,['t = ' num2str(t)],'Color',color,'FontSize',20);
           else
              text(ub(1),ub(2)+20,['t = ' num2str(t)],'Color',color,'FontSize',20);
           end

        end
    end
end
%============================ t = 2
for i=3:3
    t           = vec_t_plot(i);
    pos         = find(vec_t>=t,1);
    if ~isempty(pos)
        all_lb      = vec_windows_lb{pos};
        all_l2size  = vec_windows_l2size{pos};
        all_ub      = all_lb+2.^all_l2size-1;
        no_windows  = size(all_lb,1);
        for window=1:no_windows
           lb  = all_lb(window,:);
           ub  = all_ub(window,:);
           color = [0.8500, 0.3250, 0.0980];
           p1 = rectangle('Position',[lb(1) lb(2) ub(1)-lb(1) ub(2)-lb(2)],'LineStyle','-.','EdgeColor',color,'LineWidth',2);
           if (i == 3)
              text(ub(1)+2,ub(2)+2,['t = ' num2str(t)],'Color',color,'FontSize',20);
           else
              text(ub(1),ub(2)+20,['t = ' num2str(t)],'Color',color,'FontSize',20);
           end

        end
    end
end
%============================ t = 5
for i=4:4
    t           = vec_t_plot(i);
    pos         = find(vec_t>=t,1);
    if ~isempty(pos)
        all_lb      = vec_windows_lb{pos};
        all_l2size  = vec_windows_l2size{pos};
        all_ub      = all_lb+2.^all_l2size-1;
        no_windows  = size(all_lb,1);
        for window=1:no_windows
           lb  = all_lb(window,:);
           ub  = all_ub(window,:);
           color = [0.4940, 0.1840, 0.5560];
           p1 = rectangle('Position',[lb(1) lb(2) ub(1)-lb(1) ub(2)-lb(2)],'LineStyle',':','EdgeColor',color,'LineWidth',5);
           if (i == 3)
              text(ub(1)+2,ub(2)+2,['t = ' num2str(t)],'Color',color,'FontSize',20);
           else
              text(ub(1),ub(2)+20,['t = ' num2str(t)],'Color',color,'FontSize',20);
           end

        end
    end
end
%============================ t = 10
% for i=5:5
%     t           = vec_t_plot(i);
%     pos         = find(vec_t>=t,1);
%     if ~isempty(pos)
%         all_lb      = vec_windows_lb{pos};
%         all_l2size  = vec_windows_l2size{pos};
%         all_ub      = all_lb+2.^all_l2size-1;
%         no_windows  = size(all_lb,1);
%         for window=1:no_windows
%            lb  = all_lb(window,:);
%            ub  = all_ub(window,:);
%            color = [1, 0, 0];
%            p1 = rectangle('Position',[lb(1) lb(2) ub(1)-lb(1) ub(2)-lb(2)],'LineStyle',':','EdgeColor',color,'LineWidth',7);
%            if (i == 3)
%               text(ub(1)+2,ub(2)+2,['t = ' num2str(t)],'Color',color,'FontSize',20);
%            else
%               text(ub(1),ub(2)+20,['t = ' num2str(t)],'Color',color,'FontSize',20);
%            end
%
%         end
%     end
% end
% %==================
    xlim([0 600]);
    ylim([0 600]);
     xlabel(['Species ' num2str(1) ' (' species_name{1} ')']);
     ylabel(['Species ' num2str(2) ' (' species_name{2} ')']);
    % xlabel(['Species ' 1 ' (' species_name{1} ')']);
    % ylabel(['Species ' 2 ' (' species_name{2} ')']);
    % title('State space trajectory of FSP-QTT with parallel sliding windows');
    % title('State space trajectory of FSP-QTT with sliding windows');

    set(gca,'FontSize',20);
end
