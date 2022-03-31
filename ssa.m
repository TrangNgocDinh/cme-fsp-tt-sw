function [vec_t,vec_x,state_max,state_min,final_state] = ssa(state_start,t_step,propen_func)
    global stoich_mat_total no_species no_reactions max_output_ssa
%   Initiate the system
    state           = state_start;
    t               = 0;
    vec_t           = zeros(max_output_ssa,1);
    vec_x           = zeros(max_output_ssa,no_species);
    count           = 1;
    vec_t(count)    = t;
    vec_x(count,:)  = state;
    state_max       = state;
    state_min       = state;
%   SSA main loop
    while t<t_step
%       Warning if the output vectors are exceeded
        if count+1>max_output_ssa
            warning('SSA:ExceededCapacity');
            return;
        end
%       Calculate reaction propensities and propensity sum
        vec_prop        = zeros(no_reactions,1);
        for reaction=1:no_reactions
            vec_prop(reaction)  = propen_func(reaction,state);

        end
        sum_prop        = sum(vec_prop);
%       Draw two uniformly random numbers
        r               = rand(1,2);
%       Find next time step
        tau             = -log(r(1))/sum_prop;
%       Find next reaction
        mu              = find((cumsum(vec_prop)>=r(2)*sum_prop),1,'first');
%       Update the system and save
        t               = t+tau;
        if t<t_step
            state       = state+stoich_mat_total(mu,:);
        end
        % count           = count+1;
        % vec_t(count)    = t;
        % vec_x(count,:)  = state;
%       Update the minimum and maximum for each species along the
%       trajectory
        for species=1:no_species
            if state_max(species)<state(species)
                state_max(species)  = state(species);
            end
            if state_min(species)>state(species)
                state_min(species)  = state(species);
            end
        end
    end
    final_state = state;
end
