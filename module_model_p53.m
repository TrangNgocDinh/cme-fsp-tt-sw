    model_name                  = 'P53';
%===Number of species and reactions in the model
    global no_species no_reactions species_name
    % no_species                  = 7;
    no_species                  = 6;
    no_reactions                = 11;
    species_name                = {'RNA_nuc' 'RNA_cyt' 'MDM2_cyt' 'MDM2_nuc' 'ARF' 'p53' };
%===Species indices
    % p53                         = 1;
    % RNA_nuc                     = 2;
    % RNA_cyt                     = 3;
    % MDM2_cyt                    = 4;
    % MDM2_nuc                    = 5;
    % ARF                         = 6;
    % MDM2_nucARF                 = 7;

    RNA_nuc                     = 1;
    RNA_cyt                     = 2;
    MDM2_cyt                    = 3;
    MDM2_nuc                    = 4;
    ARF                         = 5;
    p53                         = 6;

%===Initial state
    global initial_state
    % initial_state               = zeros(1,no_species);
    % initial_state(MDM2_nuc)     = 100;
    % initial_state(ARF)          = 100;
    % initial_state(p53)          = 100;
     initial_state              = [5 30 548 152 323 326]
%===Final time
    global final_time
    % final_time                  = 120;
    final_time                  = 127.5;
    % final_time                  = 10;
%===Reaction rates
    global vec_rates
    k_p                         = 0.5;
    k_1                         = 9.963*10^(-6);
    d_p                         = 1.925*10^(-5);
    k_m                         = 1.5*10^(-3);
    k_2                         = 1.5*10^(-2);
    k_D                         = 740;
    k_0                         = 8.0*10^(-4);
    d_rc                        = 1.444*10^(-4);
    k_T                         = 1.66*10^(-2);
    k_i                         = 9.0*10^(-4);
    d_mn                        = 1.66*10^(-7);
    k_3                         = 9.963*10^(-6);
    k_a                         = 0.5;
    d_a                         = 3.209*10^(-5);
    vec_rates                   = [k_p k_1 d_p k_m k_2 k_D k_0 d_rc k_T ...
                                   k_i d_mn k_3 k_a d_a];
%===Stoichiometric matrix
    global stoich_mat
    stoich_mat                          = zeros(no_reactions,2*no_species);
%   Reaction 1 : Empty -> p53
    stoich_mat(1,no_species+p53)        = 1;
%   Reaction 2 : p53 -> Empty
    stoich_mat(2,p53)                   = 1;
%   Reaction 3 : Empty -> RNA_nuc
    stoich_mat(3,no_species+RNA_nuc)    = 1;
%   Reaction 4 : RNA_nuc -> RNA_cyt
    stoich_mat(4,RNA_nuc)               = 1;
    stoich_mat(4,no_species+RNA_cyt)    = 1;
%   Reaction 5 : RNA_cyt -> Empty
    stoich_mat(5,RNA_cyt)               = 1;
%   Reaction 6 : RNA_cyt -> RNA_cyt + MDM2_cyt
    stoich_mat(6,RNA_cyt)               = 1;
    stoich_mat(6,no_species+RNA_cyt)    = 1;
    stoich_mat(6,no_species+MDM2_cyt)   = 1;
%   Reaction 7 : MDM2_cyt -> MDM2_nuc
    stoich_mat(7,MDM2_cyt)              = 1;
    stoich_mat(7,no_species+MDM2_nuc)   = 1;
%   Reaction 8 : MDM2_nuc + MDM2_nuc -> MDM2_nuc
    stoich_mat(8,MDM2_nuc)              = 2;
    stoich_mat(8,no_species+MDM2_nuc)   = 1;
%   Reaction 9 : MDM2_nuc + ARF -> MDM2_nucARF
%   Reaction 9 : MDM2_nuc + ARF -> empty
    stoich_mat(9,MDM2_nuc)              = 1;
    stoich_mat(9,ARF)                   = 1;
    % stoich_mat(9,no_species+MDM2_nucARF)= 1;
%   Reaction 10 : Empty -> ARF
    stoich_mat(10,no_species+ARF)       = 1;
%   Reaction 11 : ARF -> Empty
    stoich_mat(11,ARF)                  = 1;
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
    % p53                         = 1;
    % RNA_nuc                     = 2;
    % RNA_cyt                     = 3;
    % MDM2_cyt                    = 4;
    % MDM2_nuc                    = 5;
    % ARF                         = 6;
    % MDM2_nucARF                 = 7;

    RNA_nuc                     = 1;
    RNA_cyt                     = 2;
    MDM2_cyt                    = 3;
    MDM2_nuc                    = 4;
    ARF                         = 5;
    p53                         = 6;
%   Constants
    k_p                         = vec_rates(1);
    k_1                         = vec_rates(2);
    d_p                         = vec_rates(3);
    k_m                         = vec_rates(4);
    k_2                         = vec_rates(5);
    k_D                         = vec_rates(6);
    k_0                         = vec_rates(7);
    d_rc                        = vec_rates(8);
    k_T                         = vec_rates(9);
    k_i                         = vec_rates(10);
    d_mn                        = vec_rates(11);
    k_3                         = vec_rates(12);
    k_a                         = vec_rates(13);
    d_a                         = vec_rates(14);
%   Check if the state is less than the left stoichiometric vector...
    if all(state>=left_propen_vec)
%       If no, then compute the propensity...
        switch reaction
            case 1 % Empty -> p53
                output = k_p;
            case 2 % p53 -> Empty
                output = state(p53)*(d_p+k_1*state(MDM2_nuc));
            case 3 % Empty -> RNA_nuc
                output = k_m+k_2*state(p53)^(1.8)/(k_D^(1.8)+state(p53)^(1.8));
            case 4 % RNA_nuc -> RNA_cyt
                output = k_0*state(RNA_nuc);
            case 5 % RNA_cyt -> Empty
                output = d_rc*state(RNA_cyt);
            case 6 % RNA_cyt -> RNA_cyt + MDM2_cyt
                output = k_T*state(RNA_cyt);
            case 7 % MDM2_cyt -> MDM2_nuc
                output = k_i*state(MDM2_cyt);
            case 8 % MDM2_nuc + MDM2_nuc -> MDM2_nuc
                output = d_mn*state(MDM2_nuc)*(state(MDM2_nuc)-1)/2;
            case 9 % MDM2_nuc + ARF -> MDM2_nucARF or empty in Huy's
                output = k_3*state(MDM2_nuc)*state(ARF);
            case 10 % Empty -> ARF
                output = k_a;
            case 11 % ARF -> Empty
                output = d_a*state(ARF);
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
    % p53                         = 1;
    % RNA_nuc                     = 2;
    % RNA_cyt                     = 3;
    % MDM2_cyt                    = 4;
    % MDM2_nuc                    = 5;
    % ARF                         = 6;
    % MDM2_nucARF                 = 7;

    RNA_nuc                     = 1;
    RNA_cyt                     = 2;
    MDM2_cyt                    = 3;
    MDM2_nuc                    = 4;
    ARF                         = 5;
    p53                         = 6;
%   Constants
    k_p                         = vec_rates(1);
    k_1                         = vec_rates(2);
    d_p                         = vec_rates(3);
    k_m                         = vec_rates(4);
    k_2                         = vec_rates(5);
    k_D                         = vec_rates(6);
    k_0                         = vec_rates(7);
    d_rc                        = vec_rates(8);
    k_T                         = vec_rates(9);
    k_i                         = vec_rates(10);
    d_mn                        = vec_rates(11);
    k_3                         = vec_rates(12);
    k_a                         = vec_rates(13);
    d_a                         = vec_rates(14);
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
        case 1 % Empty -> p53
%   full_prop = k_p, choose p53 for changing
            if (species==p53)
                vec     = k_p*ones(2^QTTdeg,1);
                output  = tt_tensor(vec,tol_tt_tensor);
            end
        case 2 % p53 -> Empty
%   full_prop = state(p53)*(d_p+k_1*state(MDM2_nuc)), change 'p53' and 'MDM2_nuc'
            if (species==p53)
                vec     = vec_initial;
                output  = tt_tensor(vec,tol_tt_tensor);
            elseif (species==MDM2_nuc)
                vec     = d_p+k_1*vec_initial;
                output  = tt_tensor(vec,tol_tt_tensor);
            end
        case 3 % Empty -> RNA_nuc
%   full_prop = k_m+k_2*state(p53)^(1.8)/(k_D^(1.8)+state(p53)^(1.8)), change 'p53'
            if (species==p53)
                vec     = k_m+k_2*(vec_initial.^1.8)./(k_D^1.8+vec_initial.^1.8);
                output  = tt_tensor(vec,tol_tt_tensor);
            end
        case 4 % RNA_nuc -> RNA_cyt
%   full_prop = k_0*state(RNA_nuc), change 'RNA_nuc'
            if (species==RNA_nuc)
                vec     = k_0*vec_initial;
                output  = tt_tensor(vec,tol_tt_tensor);
            end
        case 5 % RNA_cyt -> Empty
%   full_prop = d_rc*state(RNA_cyt), change 'RNA_cyt'
            if (species==RNA_cyt)
                vec     = d_rc*vec_initial;
                output  = tt_tensor(vec,tol_tt_tensor);
            end
        case 6 % RNA_cyt -> RNA_cyt + MDM2_cyt
%   full_prop = k_T*state(RNA_cyt), change 'RNA_cyt'
            if (species==RNA_cyt)
                vec     = k_T*vec_initial;
                output  = tt_tensor(vec,tol_tt_tensor);
            end
        case 7 % MDM2_cyt -> MDM2_nuc
%   full_prop = k_i*state(MDM2_cyt), change 'MDM2_cyt'
            if (species==MDM2_cyt)
                vec     = k_i*vec_initial;
                output  = tt_tensor(vec,tol_tt_tensor);
            end
        case 8 % MDM2_nuc + MDM2_nuc -> MDM2_nuc
%   full_prop = d_mn*state(MDM2_nuc)*(state(MDM2_nuc)-1), change 'MDM2_nuc'
            if (species==MDM2_nuc)
                vec     = 0.5*d_mn*vec_initial.*(vec_initial-1);
                output  = tt_tensor(vec,tol_tt_tensor);
            end
        case 9 % MDM2_nuc + ARF -> MDM2_nucARF
%   full_prop = k_3*state(MDM2_nuc)*state(ARF), change 'MDM2_nuc' and 'ARF'
            if (species==MDM2_nuc)
                vec     = k_3*vec_initial;
                output  = tt_tensor(vec,tol_tt_tensor);
            elseif (species==ARF)
                vec     = vec_initial;
                output  = tt_tensor(vec,tol_tt_tensor);
            end
        case 10 % Empty -> ARF
%   full_prop = k_a, choose 'ARF' for changing
            if (species==ARF)
                vec     = k_a*ones(2^QTTdeg,1);
                output  = tt_tensor(vec,tol_tt_tensor);
            end
        case 11 % ARF -> Empty
%   full_prop = d_a*state(ARF), change 'ARF'
            if (species==ARF)
                vec     = d_a*vec_initial;
                output  = tt_tensor(vec,tol_tt_tensor);
            end
    end
end
