function cme_plot_marginal_distribution(all_marg_dist,all_lower_bound,all_upper_bound,style_plot,name_plot,varargin)
    global species_name
%   Number of windows
    no_dist     = length(all_marg_dist);
%   Number of species
    no_species  = size(all_lower_bound,2);
    if no_species<=3
        no_fig_row  = 1;
        no_fig_col  = no_species;
    elseif no_species<=10
        no_fig_row  = 2;
        no_fig_col  = ceil(no_species/2);
    end

    for window=1:no_dist
        marg_dist   = all_marg_dist{window};
        lower_bound = all_lower_bound(window,:);
        upper_bound = all_upper_bound(window,:);

        for species=1:no_species
            subplot(no_fig_row,no_fig_col,species);
            if strcmp(style_plot,'plot')
                vec_x   = [lower_bound(species):upper_bound(species)];
                vec_y   = marg_dist{species};
                if window==1
                    p = plot(vec_x,vec_y,varargin{1},'LineWidth',varargin{2},'DisplayName',name_plot);hold on
                else
                    p = plot(vec_x,vec_y,varargin{1},'LineWidth',varargin{2},'HandleVisibility','off');hold on
                end
            elseif strcmp(style_plot,'semilogy')
                vec_x   = [lower_bound(species):upper_bound(species)];
                vec_y   = marg_dist{species};
                if window==1
                    p = semilogy(vec_x,vec_y,varargin{1},'LineWidth',varargin{2},'DisplayName',name_plot);hold on
                else
                    p = semilogy(vec_x,vec_y,varargin{1},'LineWidth',varargin{2},'HandleVisibility','off');hold on
                end
            elseif strcmp(style_plot,'bar')
                vec_x   = [lower_bound(species):upper_bound(species)];
                vec_y   = marg_dist{species};
                if window==1
                    p = bar(vec_x,vec_y,varargin{1},'DisplayName',name_plot);hold on
                else
                    p = bar(vec_x,vec_y,varargin{1},'HandleVisibility','off');hold on
                end
            elseif strcmp(style_plot,'scatter')
                vec_x   = [lower_bound(species):upper_bound(species)];
                vec_y   = marg_dist{species};
                if window==1
                    if nargin>10
                        p = scatter(vec_x(1:varargin{6}:end),vec_y(1:varargin{6}:end),varargin{1},varargin{2},'MarkerFaceColor',varargin{3},'MarkerEdgeColor',varargin{4},'LineWidth',varargin{5},'DisplayName',name_plot);hold on
                    else
                        p = scatter(vec_x,vec_y,varargin{1},varargin{2},'MarkerFaceColor',varargin{3},'MarkerEdgeColor',varargin{4},'LineWidth',varargin{5},'DisplayName',name_plot);hold on
                    end
                else
                    if nargin>10
                        p = scatter(vec_x(1:varargin{6}:end),vec_y(1:varargin{6}:end),varargin{1},varargin{2},'MarkerFaceColor',varargin{3},'MarkerEdgeColor',varargin{4},'LineWidth',varargin{5},'DisplayName',name_plot);hold on
                    else
                        p = scatter(vec_x,vec_y,varargin{1},varargin{2},'MarkerFaceColor',varargin{3},'MarkerEdgeColor',varargin{4},'LineWidth',varargin{5},'HandleVisibility','off');hold on
                    end
                end
            end
            xticks(vec_x);
            ylim([0 inf]);
            xlabel(['Species ' num2str(species) ' (' species_name{species} ')']);
            set(gca,'FontSize',20);
        end
    end
    subplot(no_fig_row,no_fig_col,1);
    l = legend;
    l.FontSize  = 20;
end
