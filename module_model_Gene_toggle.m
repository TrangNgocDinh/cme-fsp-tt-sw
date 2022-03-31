    model_name                  = 'Gene_toggle';
%===Number of species and reactions in the model
    global no_species no_reactions species_name
    no_species                  = 2;
    no_reactions                = 4;
    species_name                = {'U' 'V'};
%===Initial state
    global initial_state
    % initial_state               = [85 5];
    initial_state               = [100 100];
%===Final time
    global final_time
    % final_time                  = 1;
    final_time                    = 30;
%===Reaction rates
    global vec_rates
    scale                       = 100;
    d_1                         = 1.0;
    d_2                         = 1.0;
    s                           = 0.1;
    K_1                         = 1.0*scale;
    K_2                         = 1.0*scale;
    beta_1                      = 4.0*scale;
    beta_2                      = 4.0*scale;
    alpha_1                     = 0.2*scale;
    alpha_2                     = 0.2*scale;
    gamma                       = 1.0;
    eta                         = 1.0;
    vec_rates                   = [d_1 d_2 s K_1 K_2 beta_1 beta_2 alpha_1 alpha_2 gamma eta];
%===Stoichiometric matrix
    global stoich_mat
    U                           = 1;
    V                           = 2;
    stoich_mat                  = zeros(no_reactions,2*no_species);
%   Reaction 1 : U -> Empty
    stoich_mat(1,U)             = 1;
%   Reaction 2 : Empty > U
    stoich_mat(2,no_species+U)  = 1;
%   Reaction 3 : V -> Empty
    stoich_mat(3,V)             = 1;
%   Reaction 4 : Empty -> V
    stoich_mat(4,no_species+V)  = 1;
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
    global vec_rates
    left_propen_vec = stoich_mat(reaction,1:no_species);
%   Species indices
    U                           = 1;
    V                           = 2;
%   Constants
    d_1                         = vec_rates(1);
    d_2                         = vec_rates(2);
    s                           = vec_rates(3);
    K_1                         = vec_rates(4);
    K_2                         = vec_rates(5);
    beta_1                      = vec_rates(6);
    beta_2                      = vec_rates(7);
    alpha_1                     = vec_rates(8);
    alpha_2                     = vec_rates(9);
    gamma                       = vec_rates(10);
    eta                         = vec_rates(11);
%   Check if the state is less than the left stoichiometric vector...
    if all(state>=left_propen_vec)
%       If no, then compute the propensity...
        switch reaction
            case 1 % U -> Empty; function(U)
                output = (d_1+s*gamma/(1+s))*state(U);
            case 2 % Empty > U; function(V)
                output = eta*(alpha_1+beta_1*K_1^(3)/(K_1^(3)+state(V)^(3)));
            case 3 % V -> Empty; function(V)
                output = d_2*state(V);
            case 4 % Empty -> V; function(U)
                output = eta*(alpha_2+beta_2*K_2^(3)/(K_2^(3)+state(U)^(3)));
        end
%       If yes, then propensity is 0...
    else
        output      = 0;
    end
end
%============================================Partial propensity function
function output = function_propensity_partial(species,QTTdeg,reaction,varargin)
    global tol_tt_tensor vec_rates
%   Initialize the default output propensity vector
    output  = tt_tensor(ones(2^QTTdeg,1),tol_tt_tensor);
    nmax    = 2^QTTdeg-1;
%   Species indices
    U                           = 1;
    V                           = 2;
%   Constants
    d_1                         = vec_rates(1);
    d_2                         = vec_rates(2);
    s                           = vec_rates(3);
    K_1                         = vec_rates(4);
    K_2                         = vec_rates(5);
    beta_1                      = vec_rates(6);
    beta_2                      = vec_rates(7);
    alpha_1                     = vec_rates(8);
    alpha_2                     = vec_rates(9);
    gamma                       = vec_rates(10);
    eta                         = vec_rates(11);
%   If boundaries for the species are inputted, then create the
%   population vector accordingly, otherwise use the QTT degree
    if ~isempty(varargin)
        vec_initial = (varargin{1}:varargin{2})';
    else
        nmax        = 2^QTTdeg-1;
        vec_initial = (0:nmax)';
    end
%   Compute the output for each particular 'reaction' and 'species'
    switch reaction
        case 1 % U -> Empty
%   full_prop = function(U), choose 'U' for changing
            if (species==U)
                vec     = (d_1+s*gamma/(1+s))*vec_initial;
                output  = tt_tensor(vec,tol_tt_tensor);
            end
      case 2 % Empty -> U
%   full_prop = function(V), choose 'V' for changing
            if (species==V)
                vec     = eta*(alpha_1+beta_1*K_1^(3)./(K_1^(3)+vec_initial.^(3)));% Check '.operator'
                output  = tt_tensor(vec,tol_tt_tensor);
            end
        case 3 % V -> Empty
%   full_prop = function(V), choose 'V' for changing
            if (species==V)
                vec     = d_2*vec_initial;
                output  = tt_tensor(vec,tol_tt_tensor);
            end
        case 4 % Empty -> V
%   full_prop = function(U) choose 'U' for changing
            if (species==U)
                vec     = eta*(alpha_2+beta_2*K_2^(3)./(K_2^(3)+vec_initial.^(3)));% Check '.operator'
                output  = tt_tensor(vec,tol_tt_tensor);
            end
    end
end
