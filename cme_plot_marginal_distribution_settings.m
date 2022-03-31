function cme_plot_marginal_distribution_settings(lower_bound,upper_bound)
    no_species  = length(lower_bound);
    if no_species<=3
        no_fig_row  = 1;
        no_fig_col  = no_species;
    elseif no_species<=10
        no_fig_row  = 2;
        no_fig_col  = ceil(no_species/2);
    end
    for species=1:no_species
        subplot(no_fig_row,no_fig_col,species);
        xlim([lower_bound(species) upper_bound(species)]);

        no_ticks    = upper_bound(species)-lower_bound(species)+1

        if no_ticks>100
            % tick_space  = ceil(no_ticks/50)*5;
            tick_space  = ceil(no_ticks/150)*15;
            lb_ticks    = ceil(lower_bound(species)/tick_space)*tick_space;
            ub_ticks    = floor(upper_bound(species)/tick_space)*tick_space;
      elseif no_ticks>18
            tick_space  = ceil(no_ticks/20)*2;
            lb_ticks    = ceil(lower_bound(species)/tick_space)*tick_space;
            ub_ticks    = floor(upper_bound(species)/tick_space)*tick_space;
        else
            tick_space  = 1;
            lb_ticks    = lower_bound(species);
            ub_ticks    = upper_bound(species);
        end
        xticks([lb_ticks:tick_space:ub_ticks]);
    end
end
