function cme_cartoon_plot_ss_trajectory(vec_t,vec_windows_lb,vec_windows_l2size)
    vec_t_plot  = [0 1 7 10];

    for i=1:length(vec_t_plot)
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
                p1  = plot([lb(1) lb(1)],[lb(2) ub(2)],'LineWidth',3);hold on
                if window>=2
                    p1.Color    = col;
                else
                    col = p1.Color;
                end
                plot([ub(1) ub(1)],[lb(2) ub(2)],'Color',col,'LineWidth',3);hold on
                plot([lb(1) ub(1)],[lb(2) lb(2)],'Color',col,'LineWidth',3);hold on
                plot([lb(1) ub(1)],[ub(2) ub(2)],'Color',col,'LineWidth',3);hold on
                text(ub(1)+2,ub(2)+30,['t = ' num2str(t)],'Color',col,'FontSize',20);
            end
        end
    end

    xlim([0 600]);
    ylim([0 600]);

    title('State space trajectory of FSP-QTT with parallel sliding windows');

    set(gca,'FontSize',20);
end
