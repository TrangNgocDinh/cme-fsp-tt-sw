%=====================================================PARAMETERS FOR SSA
    global max_output_ssa no_SSA_traj
%   Maximum number of entries in output
    max_output_ssa  = 100000000;
%   Number of SSA trajectories for producing the distribution
%     no_SSA_traj     = 1000;
    no_SSA_traj     = 10000;
    % no_SSA_traj     = 100000;
%===================PARAMETERS FOR ONE-STEP SPARSE/DENSE TRADITIONAL CME
    global vec_ss_minsize vec_ss_maxsize
%   Tolerance for Expokit
    tol_expokit     = 10^-10;
%   Krylov degree
    deg_expokit     = 50;
%   State space
    if strcmp(model_name,'Michaelis-Menten')
        vec_ss_minsize  = [0 0 0 0];
        vec_ss_maxsize  = [3 7 3 3];
    elseif strcmp(model_name,'Michaelis-Menten-big')
        vec_ss_minsize  = [0 0 0 0];
        vec_ss_maxsize  = [15 15 15 10];
    elseif strcmp(model_name,'Gene_toggle')
        vec_ss_minsize  = [0 0];
        vec_ss_maxsize  = [(2^9-1) (2^9-1)];
    elseif strcmp(model_name,'P53')
        vec_ss_minsize  = [0 0 0 0 0 0 0];
        vec_ss_maxsize  = [50 1 25 50 85 50 3];
    end
%   Total number of states
    no_states       = prod(vec_ss_maxsize-vec_ss_minsize+1);
%=================PARAMETERS FOR ADAPTIVE SPARSE TRADITIONAL SSA-FSP CME
    global vec_ss_minsize_SSA_initial vec_ss_maxsize_SSA_initial
    global no_SSA_ss_sparse tol_reduce_ss_sparse
    global t_step_initial_SSA
    global tol_local_fsp_SSA tol_expokit_SSA deg_expokit_SSA
%   Tolerance for reducing state space
    tol_reduce_ss_sparse        = 10^-15;
%   Number of SSA runs for each state on the state space boundary
%     if strcmp(model_name,'Michaelis-Menten')
% %         no_SSA_ss_sparse            = 5;
%           no_SSA_ss_sparse            = 10;
%     elseif strcmp(model_name,'Michaelis-Menten-big')
%         no_SSA_ss_sparse            = 10;
%     elseif strcmp(model_name,'Gene_toggle')
% %         no_SSA_ss_sparse            = 5;
%         no_SSA_ss_sparse            = 10;
%     elseif strcmp(model_name,'P53')
%         % no_SSA_ss_sparse            = 5;
%         no_SSA_ss_sparse            = 10;
%     end
    no_SSA_ss_sparse            = 10;
%   Tolerance for Expokit
    tol_expokit_SSA             = 10^-7;
%   Krylov degree for Expokit
    deg_expokit_SSA             = 50;
%   Tolerance for local error at each time step
    tol_local_fsp_SSA           = 10^-5;
%   Initial guess for the state space
    if strcmp(model_name,'Michaelis-Menten')
        vec_ss_minsize_SSA_initial  = initial_state;
        vec_ss_maxsize_SSA_initial  = initial_state;
    elseif strcmp(model_name,'Michaelis-Menten-big')
        vec_ss_minsize_SSA_initial  = initial_state;
        vec_ss_maxsize_SSA_initial  = initial_state;
    elseif strcmp(model_name,'Gene_toggle')
        vec_ss_minsize_SSA_initial  = initial_state;
        vec_ss_maxsize_SSA_initial  = initial_state;
    elseif strcmp(model_name,'P53')
        vec_ss_minsize_SSA_initial  = initial_state;
        vec_ss_maxsize_SSA_initial  = initial_state;
    elseif  strcmp(model_name,'Goutsias')
         vec_ss_minsize_SSA_initial  = initial_state;
         vec_ss_maxsize_SSA_initial  = initial_state;
    end
%   Initial guess for the time step
    if strcmp(model_name,'Michaelis-Menten')
        t_step_initial_SSA          = 0.5;
    elseif strcmp(model_name,'Michaelis-Menten-big')
        t_step_initial_SSA          = 0.5;
    elseif strcmp(model_name,'Gene_toggle')
        t_step_initial_SSA          = 0.5;
    elseif strcmp(model_name,'P53')
        t_step_initial_SSA          = 0.005;
    elseif strcmp(model_name,'Goutsias')
        t_step_initial_SSA          = 0.005;
    end
%=====================================================PARAMETERS FOR QTT
%==================================UNIFORMIZATION + AMEN, AND TT-EXPOKIT
    global tol_tt_tensor tol_fsp tol_unif tol_amen sc_factor tol_qtt_mat
    global FSP_qtt_size FSP_qtt_size_initial FSP_qtt_size
    global tol_expokit_tt deg_expokit_tt t_step_initial_qtt_expokit
    global tol_tt_mvk4
    global tol_fsp_probsum
%   FSP size for QTT (powers of two)
    if strcmp(model_name,'Michaelis-Menten')
        FSP_qtt_size_initial    = [2 3 2 2];
        % FSP_qtt_size_initial    = [5 8 5 5];
    elseif strcmp(model_name,'Michaelis-Menten-big')
        FSP_qtt_size_initial    = [4 4 4 3];
    elseif strcmp(model_name,'Gene_toggle')
        % FSP_qtt_size_initial    = [7 7];
        FSP_qtt_size_initial    = [9 9];
        % FSP_qtt_size_initial    = [7 6];
    elseif strcmp(model_name,'P53')
        % FSP_qtt_size_initial    = [6 2 5 6 7 6 2];
        % FSP_qtt_size_initial    = [7 2 2 2 7 7 2];
        % FSP_qtt_size_initial    = [2 2 2 7 7 7];
        FSP_qtt_size_initial    = [2 2 2 9 9 9];
    elseif strcmp(model_name,'Goutsias')
        % FSP_qtt_size_initial    = [3 3 2 3 2 2];
        % FSP_qtt_size_initial    = [4 4 3 4 3 3];
        % FSP_qtt_size_initial    = [3 3 3 3 3 3];
          % FSP_qtt_size_initial    = [5 6 2 5 2 2];
          FSP_qtt_size_initial    = [6 7 3 6 3 3];
        % [2 6 0 2 0 0]
        % [20 60 0 20 0 0];
    end
    FSP_qtt_size    = FSP_qtt_size_initial;
%   Tolerance for TT_tensor rounding
    % tol_tt_tensor   = 10^-10;
    tol_tt_tensor   = 10^-4;
%   Tolerance for FSP
    % tol_fsp         = 10^-4;
    tol_fsp         = 10^-4;
%   Tolerance for probability sum
    tol_fsp_probsum = 0.05;
    % tol_fsp_probsum = 0.01;
%   Tolerance for Uniformization
    % tol_unif        = 10^-4;
    tol_unif        = 10^-4;
%   Tolerance for AMEn
    % tol_amen        = 10^-5;
    tol_amen        = 10^-4;
%   Tolerance for TT-Expokit
    % tol_expokit_tt  = 10^-5;
    tol_expokit_tt  = 10^-4;
%   Tolerance for tt-mvk4 in TT-Expokit
    % tol_tt_mvk4     = 10^-10;
    tol_tt_mvk4     = 10^-4;
%   Krylov degree for TT-Expokit
    deg_expokit_tt  = 20;
%   Axiliary factor when choosing the step-size in Uniformiztion
%   nstep = [alpha*t/sc_factor], chosen as 100 in MacNamara et al.
    if strcmp(model_name,'P53')
        sc_factor   = 10;
    else
        sc_factor   = 100;
    end
%   Tolerance for building the CME matrix
    % tol_qtt_mat     = 10^-14;
    tol_qtt_mat     = 10^-4;
%   Initial guess for the time step for TT-Expokit
    if strcmp(model_name,'Michaelis-Menten')
        t_step_initial_qtt_expokit  = 2.0;
    elseif strcmp(model_name,'Michaelis-Menten-big')
        t_step_initial_qtt_expokit  = 2.0;
    elseif strcmp(model_name,'Gene_toggle')
        t_step_initial_qtt_expokit  = 2.0;
    elseif strcmp(model_name,'P53')
        t_step_initial_qtt_expokit  = 0.5;
    end
%======================================================PARAMETERS FOR TT
    global tol_reduce_ss_tt
%   Tolerance for reducing state space
    tol_reduce_ss_tt    = 10^-4;
