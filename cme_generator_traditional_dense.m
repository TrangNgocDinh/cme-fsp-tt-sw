function [A_full,state_space_list] = cme_generator_traditional_dense(propen_func,varargin)
    global vec_ss_minsize vec_ss_maxsize stoich_mat
%---Produces a dense CME matrix
%   [A_full,state_space_list] = cme_generator_traditional_dense(stoich_mat,propen_func,vec_minsize,vec_maxsize,varargin)
%---Input:
%*  stoich_mat          = stoichiometric matrix, by size (#reactions) X (2 #species)
%                         containing both left and right stoichiometric values
%*  propen_fun          = propensity function
%*  vec_minsize         = minimum bounds for each species in the state space
%*  vec_maxsize         = maximum bounds for each species in the state space
%*  varargin            = unnecessary parameters
%---Output:
%*  A_full              = dense CME matrix
%*  state_space_list    = state space list for the CME matrix
    [A,state_space_list]    = cme_generator_traditional_sparse(propen_func);
    A_full                  = full(A);
end
