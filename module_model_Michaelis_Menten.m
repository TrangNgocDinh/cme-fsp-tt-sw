    model_name                  = 'Michaelis-Menten';
%===Number of species and reactions in the model
    global no_species no_reactions species_name
    no_species                  = 4;
    no_reactions                = 3;
    species_name                = {'S' 'E' 'ES' 'P'};
%===Initial state
    global initial_state
    initial_state               = [30 70 0 0];
    % initial_state               = [3 7 0 0];
%===Final time
    global final_time
    final_time                  = 10;
    % final_time                  = 5;
%===Reaction rates
    global vec_rates
    vec_rates                   = [1 1 0.1];
%===Stoichiometric matrix
    global stoich_mat
    S                           = 1;
    E                           = 2;
    ES                          = 3;
    P                           = 4;
    stoich_mat                  = zeros(no_reactions,2*no_species);
%   S + E   -> ES
    stoich_mat(1,S)             = 1;
    stoich_mat(1,E)             = 1;
    stoich_mat(1,no_species+ES) = 1;
%   ES      -> S + E
    stoich_mat(2,ES)            = 1;
    stoich_mat(2,no_species+S)  = 1;
    stoich_mat(2,no_species+E)  = 1;
%   ES      -> P + E
    stoich_mat(3,ES)            = 1;
    stoich_mat(3,no_species+P)  = 1;
    stoich_mat(3,no_species+E)  = 1;
%===Left stoichiometric matrix (for reactants)
    stoich_mat_left             = stoich_mat(:,1:no_species);
%===Right stoichiometric matrix (for products)
    stoich_mat_right            = stoich_mat(:,no_species+1:2*no_species);
%===Total stoichiometric matrix (= right - left)
    global stoich_mat_total
    stoich_mat_total            = stoich_mat_right-stoich_mat_left;
%===Full propensity function
    propen_func                 = @function_propensity_full;
%===Partial propensity function
    propen_func_partial         = @function_propensity_partial;
%===============================================Full propensity function
function output = function_propensity_full(reaction,state)
    global no_species no_reactions vec_rates stoich_mat
    left_propen_vec = stoich_mat(reaction,1:no_species);
%   Check if the state is less than the left stoichiometric vector...
    if all(state>=left_propen_vec)
%       If no, then compute the propensity...
        output      = vec_rates(reaction);
        for species=1:no_species
            output  = output*nchoosek(state(species),stoich_mat(reaction,species));
        end
%       If yes, then propensity is 0...
    else
        output      = 0;
    end
end
%============================================Partial propensity function
function output = function_propensity_partial(species,QTTdeg,reaction,varargin)
    global tol_tt_tensor vec_rates
%   Initialize the output propensity vector
    output  = tt_tensor(ones(2^QTTdeg,1),tol_tt_tensor);
    nmax    = 2^QTTdeg-1;
%   Species indices
    S       = 1;
    E       = 2;
    ES      = 3;
    P       = 4;
%   If boundaries for the species are inputted, then create the
%   population vector accordingly, otherwise use the QTT degree
    if ~isempty(varargin)
        vec_initial = (varargin{1}:varargin{2})';
    else
        nmax        = 2^QTTdeg-1;
        vec_initial = (0:nmax)';
    end
%   Compute the partial propensity vectors
    switch reaction
        case 1
%           S + E   -> SE
            if (species==S)
                vec     = vec_rates(1)*vec_initial;
                output  = tt_tensor(vec,tol_tt_tensor);
            elseif (species==E)
                vec     = vec_initial;
                output  = tt_tensor(vec,tol_tt_tensor);
            end
        case 2
%           ES      -> S + E
            if (species==ES)
                vec     = vec_rates(2)*vec_initial;
                output  = tt_tensor(vec,tol_tt_tensor);
            end
        case 3
%           ES      -> P + E
            if (species==ES)
                vec     = vec_rates(3)*vec_initial;
                output  = tt_tensor(vec,tol_tt_tensor);
            end
    end
end
