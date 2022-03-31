function [no_windows_new,all_v_new,all_lb_new,all_l2size_new] = cme_update_ss_qttpsw(no_windows,all_v,all_lb_now,all_l2size_now,logic_redo,t_step,propen_func)
    global no_species
%-----------------Get the marginal distributions for all current windows
    all_marginal    = cell(no_windows,1);
    for window=1:no_windows
        lb_now      = all_lb_now(window,:);
        l2size_now  = all_l2size_now(window,:);
        ub_now      = lb_now+2.^l2size_now-1;
        v           = all_v{window};
        [lb,ub,marginal]        = cme_marginal_qttmw(lb_now,ub_now,v);
        all_marginal{window}    = marginal;
        all_marginal_lb{window} = lb;
        all_marginal_ub{window} = ub;
    end
%----------------------------Reduce and check each window for separation
    no_windows_new      = 0;
    all_lb_next         = zeros(1,no_species);
    all_ub_next         = zeros(1,no_species);
    all_l2size_next     = zeros(1,no_species);
    window_new_origin   = zeros(1,1);
    for window=1:no_windows
%       Reduce the window
        lb_now          = all_lb_now(window,:);
        l2size_now      = all_l2size_now(window,:);
        ub_now          = lb_now+2.^l2size_now-1;
        marginal        = all_marginal{window};
        marginal_lb     = all_marginal_lb{window};
        marginal_ub     = all_marginal_ub{window};
        if logic_redo==0
            mode_now        = ub_now - lb_now + 1;
            v               = all_v{window};
            v_physic        = tt_reshape(v,mode_now);
            [lb_tmp,ub_tmp] = cme_update_ss_qttmw_reduce(lb_now,ub_now,v_physic);
        else
            no_windows_new  = no_windows_new+1;
            lb_tmp          = lb_now;
            ub_tmp          = ub_now;
            all_lb_next(no_windows_new,:)    = lb_now;
            all_ub_next(no_windows_new,:)    = ub_now;
            window_new_origin(no_windows_new)   = window;
            % continue
        end
%       Check if the window can be separated
        logic_separated_window  = 0;
        all_lb_next_species     = zeros(1,no_species);
        all_ub_next_species     = zeros(1,no_species);
        for species=1:no_species
            if logic_separated_window==0
                marginal_species    = marginal{species};
%               Find all supports in the window
                pos_tmp             = find(marginal_species>=10^-4);

                lb_species_tmp      = pos_tmp(1);
                ub_species_tmp      = pos_tmp(1);

                no_miniwindows      = 0;

                for i=2:length(pos_tmp)
                    if (pos_tmp(i)==ub_species_tmp+1)&&(i<length(pos_tmp))
                        ub_species_tmp                  = ub_species_tmp+1;
                    else
                        no_miniwindows                  = no_miniwindows+1;
                        lb_next                         = lb_tmp;
                        ub_next                         = ub_tmp;
                        lb_next(species)                = lb_species_tmp+marginal_lb(species)-1;
                        ub_next(species)                = ub_species_tmp+marginal_lb(species)-1;
                        all_lb_next_species(no_miniwindows,:)   = lb_next;
                        all_ub_next_species(no_miniwindows,:)   = ub_next;
                        if no_miniwindows>=2
                            logic_separated_window  = 1;
                        end


                        lb_species_tmp                  = pos_tmp(i);
                        ub_species_tmp                  = pos_tmp(i);
                    end
                end
            end
        end
        if no_miniwindows<=1
            no_windows_new          = no_windows_new+1;
            all_lb_next(no_windows_new,:)        = lb_tmp;
            all_ub_next(no_windows_new,:)        = ub_tmp;
            window_new_origin(no_windows_new)    = window;
        else
            for i=1:no_miniwindows
                no_windows_new                      = no_windows_new+1;
                all_lb_next(no_windows_new,:)       = all_lb_next_species(i,:);
                all_ub_next(no_windows_new,:)       = all_ub_next_species(i,:);
                window_new_origin(no_windows_new)   = window;
            end
        end
    end
%--------------------------------------------Expand each window with SSA
    for i=1:no_windows_new
        lb_tmp              = all_lb_next(i,:);
        ub_tmp              = all_ub_next(i,:);
        [lb_new,ub_new]     = cme_update_ss_qttmw_SSA(lb_tmp,ub_tmp,t_step,propen_func);
        all_lb_next(i,:)    = lb_new;
        all_ub_next(i,:)    = ub_new;
    end
%----------------------------------Transform each window into QTT format
    for i=1:no_windows_new
        lb_tmp              = all_lb_next(i,:);
        ub_tmp              = all_ub_next(i,:);
        [lb_tmp,ub_tmp,l2size_tmp]  = transform_to_QTT(lb_tmp,ub_tmp);
        all_lb_next(i,:)    = lb_tmp;
        all_ub_next(i,:)    = ub_tmp;
        all_l2size_next(i,:)= l2size_tmp;
    end
%--------------------------------Check if any two windows are overlapped
    flag_overlapped = 1;
    while flag_overlapped==1
        flag_overlapped = 0;
        for window1=1:no_windows_new
            for window2=window1+1:no_windows_new
%               Check if these two windows are overlapped
                no_nonoverlapped_species    = 0;
                for species=1:no_species
                    if (all_ub_next(window1,species)<all_lb_next(window2,species))||(all_ub_next(window2,species)<all_lb_next(window1,species))
                        no_nonoverlapped_species    = no_nonoverlapped_species+1;
                    end
                end
                if no_nonoverlapped_species==0
%                   These two windows are overlapped
                    flag_overlapped = 1;
                    if window_new_origin(window1)==window_new_origin(window2)
%                       The overlapped windows are from the same source window...
%                       Combine the windows
                        lb_tmp  = min(all_lb_next(window1,:),all_lb_next(window2,:));
                        ub_tmp  = max(all_ub_next(window1,:),all_ub_next(window2,:));
%                       Transform the new window into QTT format
                        [lb_tmp,ub_tmp,l2size_tmp]  = transform_to_QTT(lb_tmp,ub_tmp);
%                       Save the new window as window 1
                        all_lb_next(window1,:)      = lb_tmp;
                        all_ub_next(window1,:)      = ub_tmp;
                        all_l2size_next(window1,:)  = l2size_tmp;
%                       Delete window 2
                        no_windows_new              = no_windows_new-1;
                        all_lb_next(window2,:)      = [];
                        all_ub_next(window2,:)      = [];
                        all_l2size_next(window2,:)  = [];
                        window_new_origin(window2)  = [];
                    else
%>>>>>>>>>>>>>>>>>>>>>>>
%>>>>>>>>>>>>>>>>>>>>>>>
%>>>>>>>>>>>>>>>>>>>>>>>
%>>>>>>>>>>>>>>>>>>>>>>>
%>>>>>>>>>>>>>>>>>>>>>>>
                    end
                    break;
                end
            end
            if flag_overlapped==1
                break;
            end
        end
    end
%--------------------Update the probability distribution for new windows
    all_v_new   = cell(no_windows_new,1);
    for window_new=1:no_windows_new
        window_old  = window_new_origin(window_new);

        lb_old      = all_lb_now(window_old,:);
        l2size_old  = all_l2size_now(window_old,:);

        lb_new      = all_lb_next(window_new,:);
        l2size_new  = all_l2size_next(window_new,:);

        v                   = all_v{window_old};
        v                   = cme_update_prob_qttmw(v,lb_old,l2size_old,lb_new,l2size_new);
        all_v_new{window_new}   = v;
    end

    all_lb_new      = all_lb_next;
    all_l2size_new  = all_l2size_next;

    fprintf('\nTotal number of windows = %d\n',no_windows_new);
    for i=1:no_windows_new
        fprintf('\n[%d %d] ---> [%d %d] from window %d\n\n',all_lb_next(i,1),all_lb_next(i,2),all_ub_next(i,1),all_ub_next(i,2),window_new_origin(i));
    end

end

function [lb_tmp,ub_tmp,l2size_tmp] = transform_to_QTT(lb_tmp,ub_tmp)
    global no_species
    l2size_tmp          = max(ceil(log2(ub_tmp-lb_tmp+1)),2);
    for species=1:no_species
%       Calculate how many states will be redundant because of QTT storage
        excess          = 2.^l2size_tmp(species)-ub_tmp(species)+lb_tmp(species)-1;
%       Distribute the excess states evenly on both sides of the bounds
        lb_tmp(species) = max(0,lb_tmp(species)-floor(excess/2));
        ub_tmp(species) = lb_tmp(species)+2^l2size_tmp(species)-1;
    end
end
