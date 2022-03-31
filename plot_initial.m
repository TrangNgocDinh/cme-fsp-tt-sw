function plot_initial(initial_state,lower_bound,upper_bound)
    global species_name final_time

    x           =  lower_bound(1)-1;
    y           =  lower_bound(2)-1;

    grid on;


    scatter(initial_state(1), initial_state(2),'MarkerFaceColor',[1 0 0]);

    set(gca, 'XLim', [0, 600], 'YLim', [0, 600], ...%'YDir', 'reverse', ...
    'PlotBoxAspectRatio', [700, 700, 1],'FontSize',20);
    % set(gca, 'XLim', [0, 600], 'YLim', [0, 600]);

    hold on;

    % xlabel(['Species ' num2str(1) ' (' species_name{1} ')']);
    % ylabel(['Species ' num2str(2) ' (' species_name{2} ')']);
    xlabel([ species_name{1} ]);
    ylabel([species_name{2} ]);
end
