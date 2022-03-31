model_name                  = 'Goutsias';
%===Number of species and reactions in the model
global no_species no_reactions species_name
no_species                  = 6;
no_reactions                = 10;
species_name                = {'M' 'D' 'RNA' 'DNA' 'DNA-D' 'DNA-2D'};
%===Species indices
M                           = 1;
D                           = 2;
RNA                         = 3;
DNA                         = 4;
DNA_D                       = 5;
DNA_2D                      = 6;
%===Initial state
global initial_state
% initial_state               = [20 60 0 20 0 0];
initial_state               = [2 6 0 2 0 0];
%===Final time
global final_time
% final_time                  = 300;
final_time                  = 100;
% final_time                  = 30;
% final_time                  = 1000;
%===Reaction rates
Avo                         = 6.0221415*10^(23);    % Avogadro's Number
% L                           = 1;                    % System constant????
Vol                           = 10^(-15);           % System volume

global vec_rates
vec_rates                   = zeros(1,no_reactions);

vec_rates(1)                = 0.043;
vec_rates(2)                = 0.0007;
vec_rates(3)                = 0.0715;
vec_rates(4)                = 0.0039;
vec_rates(5)                = 0.012*10^(9)/(Avo*Vol);
vec_rates(6)                = 0.4791;
vec_rates(7)                = 0.00012*10^(9)/(Avo*Vol);
vec_rates(8)                = 0.8765*10^(-11);
vec_rates(9)                = 0.05*10^(9)/(Avo*Vol);
vec_rates(10)               = 0.5;

%===Stoichiometric matrix
global stoich_mat
stoich_mat                  = zeros(no_reactions,2*no_species);
%   Reaction 1 : RNA -> RNA + M
  stoich_mat(1,RNA)                = 1;
  stoich_mat(1,no_species+RNA)     = 1;
  stoich_mat(1,no_species+M)       = 1;
%   Reaction 2 : M -> Empty
  stoich_mat(2,M)                  = 1;
%   Reaction 3 : DNA_D -> RNA + DNA_D
  stoich_mat(3,DNA_D)              = 1;
  stoich_mat(3,no_species+RNA)     = 1;
  stoich_mat(3,no_species+DNA_D)   = 1;
%   Reaction 4 : RNA -> Empty
  stoich_mat(4,RNA)                = 1;
%   Reaction 5 : DNA + D -> DNA_D
  stoich_mat(5,DNA)                = 1;
  stoich_mat(5,D)                  = 1;
  stoich_mat(5,no_species+DNA_D)   = 1;
%   Reaction 6 : DNA_D -> DNA + D
  stoich_mat(6,DNA_D)              = 1;
  stoich_mat(6,no_species+DNA)     = 1;
  stoich_mat(6,no_species+D)       = 1;
%   Reaction 7 : DNA_D + D -> DNA_2D
  stoich_mat(7,DNA_D)              = 1;
  stoich_mat(7,D)                  = 1;
  stoich_mat(7,no_species+DNA_2D)  = 1;
%   Reaction 8 : DNA_2D -> DNA_D + D
  stoich_mat(8,DNA_2D)             = 1;
  stoich_mat(8,no_species+DNA_D)   = 1;
  stoich_mat(8,no_species+D)       = 1;
%   Reaction 9 : M + M -> D
  stoich_mat(9,M)                  = 2;
  stoich_mat(9,no_species+D)       = 1;
%   Reaction 10 : D -> M + M
  stoich_mat(10,D)                 = 1;
  stoich_mat(10,no_species+M)      = 2;
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
  M            = 1;
  D            = 2;
  RNA          = 3;
  DNA          = 4;
  DNA_D        = 5;
  DNA_2D       = 6;
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
%   Reaction 1 : RNA -> RNA + M
        if (species==RNA)
            vec     = vec_rates(1)*vec_initial;
            output  = tt_tensor(vec,tol_tt_tensor);
        end
    case 2
%   Reaction 2 : M -> Empty
        if (species==M)
            vec     = vec_rates(2)*vec_initial;
            output  = tt_tensor(vec,tol_tt_tensor);
        end
    case 3
%   Reaction 3 : DNA_D -> RNA + DNA_D
        if (species==DNA_D)
            vec     = vec_rates(3)*vec_initial;
            output  = tt_tensor(vec,tol_tt_tensor);
        end
    case 4
%   Reaction 4 : RNA -> Empty
        if (species==RNA)
            vec     = vec_rates(4)*vec_initial;
            output  = tt_tensor(vec,tol_tt_tensor);
        end
    case 5
%   Reaction 5 : DNA + D -> DNA_D
        if (species==DNA)
            vec     = vec_rates(5)*vec_initial;
            output  = tt_tensor(vec,tol_tt_tensor);
      elseif (species==D)
            vec     = vec_initial;
            output  = tt_tensor(vec,tol_tt_tensor);
        end
  case 6
%   Reaction 6 : DNA_D -> DNA + D
        if (species==DNA_D)
            vec     = vec_rates(6)*vec_initial;
            output  = tt_tensor(vec,tol_tt_tensor);
        end
  case 7
%   Reaction 7 : DNA_D + D -> DNA_2D
        if (species==DNA_D)
            vec     = vec_rates(7)*vec_initial;
            output  = tt_tensor(vec,tol_tt_tensor);
      elseif (species==D)
            vec     = vec_initial;
            output  = tt_tensor(vec,tol_tt_tensor);
        end
  case 8
%   Reaction 8 : DNA_2D -> DNA_D + D
        if (species==DNA_2D)
            vec     = vec_rates(8)*vec_initial;
            output  = tt_tensor(vec,tol_tt_tensor);
        end
  case 9
%   Reaction 9 : M + M -> D
        if (species==M)
            vec     = max( (1/2)*vec_rates(9)*vec_initial.*(vec_initial-1),0 );
            output  = tt_tensor(vec,tol_tt_tensor);
        end
  case 10
%   Reaction 10 : D -> M + M
        if (species==D)
            vec     = vec_rates(10)*vec_initial;
            output  = tt_tensor(vec,tol_tt_tensor);
        end
end
end
