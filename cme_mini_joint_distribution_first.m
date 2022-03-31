% function cme_plot_joint_distribution(joint_dist_tt,lower_bound,upper_bound,data_type,name_plot)
function cme_mini_joint_distribution_first(joint_dist_tt,lower_bound,upper_bound)
    % if strcmp(data_type,'one-window')
    %     joint_dist_tt_window    = joint_dist_tt{1};
    %     plot_one_window(joint_dist_tt_window,lower_bound,upper_bound,name_plot);
    % elseif strcmp(data_type,'multiple-windows')
        for window=1:size(lower_bound,1)
            joint_dist_tt_window    = joint_dist_tt{window};
            lower_bound_window      = lower_bound(window,:);
            upper_bound_window      = upper_bound(window,:);
            % plot_one_window(joint_dist_tt_window,lower_bound_window,upper_bound_window,name_plot); hold on
            plot_one_window(joint_dist_tt_window,lower_bound_window,upper_bound_window); hold on
        end
    % end
end

% function plot_one_window_middle(joint_dist_tt,lower_bound,upper_bound,name_plot)
function plot_one_window(joint_dist_tt,lower_bound,upper_bound)
    global species_name final_time

    joint_dist  = zeros(upper_bound(1)-lower_bound(1)+1,upper_bound(2)-lower_bound(2)+1);

    for row=1:size(joint_dist,1)
        for col=1:size(joint_dist,2)
            joint_dist(row,col) = joint_dist_tt(row,col);
        end
    end

    [m,n]       = size(joint_dist);

    nonzeroInd  = find(joint_dist>=10^-5);
    [x,y]       = ind2sub([m n],nonzeroInd);

    x           = x+lower_bound(1)-1;
    y           = y+lower_bound(2)-1;

    grid on;
    P           = joint_dist(nonzeroInd);
    P           = log10(P);

    hp = patch(x,y,P,'Marker','s','MarkerFaceColor','flat','MarkerSize',4,'EdgeColor','none','FaceColor','none');
    % set(gca,'PlotBoxAspectRatio',[700,700,1],'FontSize',20);
    set(gca, 'XLim', [0, 600], 'YLim', [0, 600], ...%'YDir', 'reverse', ...
    'PlotBoxAspectRatio', [700, 700, 1],'FontSize',20);
    % set(gca, 'XLim', [0, 600], 'YLim', [0, 600]);
    colorbar;
    hold on;

    % plot([lower_bound(1) lower_bound(1)],[lower_bound(2) upper_bound(2)],'-r','LineWidth',3);hold on
    % plot([upper_bound(1) upper_bound(1)],[lower_bound(2) upper_bound(2)],'-r','LineWidth',3);hold on
    % plot([lower_bound(1) upper_bound(1)],[lower_bound(2) lower_bound(2)],'-r','LineWidth',3);hold on
    % plot([lower_bound(1) upper_bound(1)],[upper_bound(2) upper_bound(2)],'-r','LineWidth',3);hold on

    % title(name_plot)
    % title({name_plot;['time = ', num2str(final_time)]})
    % title({['time = ', num2str(final_time)]})
    % xlim([0 600]);
    % ylim([0 600]);
    % xlabel(['Species ' num2str(1) ' (' species_name{1} ')']);
    % ylabel(['Species ' num2str(2) ' (' species_name{2} ')']);
   % xlabel([' species_name{1} ']);
   % ylabel([' species_name{2} ']);
   xlabel([ species_name{1} ]);
   ylabel([species_name{2} ]);
end
