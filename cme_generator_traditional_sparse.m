function [A,state_space_list] = cme_generator_traditional_sparse(propen_func,varargin)
    global vec_ss_minsize vec_ss_maxsize stoich_mat
%---Produces a sparse CME matrix
%   [A,state_space_list] = cme_generator_traditional_sparse(stoich_mat,propen_func,vec_minsize,vec_maxsize,varargin)
%---Input:
%*  stoich_mat          = stoichiometric matrix, by size (#reactions) X (2 #species)
%                         containing both left and right stoichiometric values
%*  propen_fun          = propensity function
%*  vec_minsize         = minimum bounds for each species in the state space
%*  vec_maxsize         = maximum bounds for each species in the state space
%*  varargin            = unnecessary parameters
%---Output:
%*  A                   = sparse CME matrix
%*  state_space_list    = state space list for the CME matrix
    if nargin==1, maxsize=0; end;

    [i,j,aij,n,nz,state_space_list]=matrix_generator(stoich_mat,propen_func,vec_ss_minsize,vec_ss_maxsize,varargin{:});
    A           = sparse(i,j,aij);
end

function [i,j,aij,n,nz,state_space_list] = matrix_generator(s,c,vec_minsize,vec_maxsize,varargin)
%---Find total number of states, number of reactions (pd) and number of
%---species (sd)
    no_total_states = prod(vec_maxsize-vec_minsize+1);
    n               = no_total_states;
    [pd,sd]         = size(s);
    sd              = sd/2;
%---Build the state space list
    state_space_list    = zeros(no_total_states,sd);
    for i=1:no_total_states
        state_space_list(i,:)   = index_to_state(i,vec_minsize,vec_maxsize);
    end
%---Build the CME matrix
    nz  = 0;
    for current_index=1:no_total_states
        if (mod(current_index,1000)==0)
            fprintf('Bulding sparse CME matrix; at state %d out of %d\n',current_index,no_total_states);
        end
        current_state   = index_to_state(current_index,vec_minsize,vec_maxsize);
        [rs,props]      = reachable_states(current_state,s,c,varargin{:});
%       Compute the diagonal entry
        nz      = nz+1;
        i(nz)   = current_index;
        j(nz)   = current_index;
        aij(nz) = -sum(props);
%       Compute the off-diagonal entries
        for h=1:pd
            if (rs(h,1)>=0)
                next_state  = rs(h,:);
                next_index  = state_to_index(next_state,vec_minsize,vec_maxsize);
                if next_index>0
                    nz      = nz+1;
                    i(nz)   = next_index;
                    j(nz)   = current_index;
                    aij(nz) = props(h);
                end
            end
        end
    end
end

function index = state_to_index(state,vec_minsize,vec_maxsize)
%---Check if the state is within the bounds
    for i=1:length(state)
        if (state(i)<vec_minsize(i))|(state(i)>vec_maxsize(i))
            index   = -1;
            return
        end
    end
%---Conversion from state to index
    index   = 1;
    for i=1:length(state)
        C   = 1;
        for j=1:i-1
            C   = C*(vec_maxsize(j)-vec_minsize(j)+1);
        end
        index   = index + (state(i)-vec_minsize(i))*C;
    end
end

function state = index_to_state(index,vec_minsize,vec_maxsize)
%---Check if the index is within the bounds
    no_total_states = prod(vec_maxsize-vec_minsize+1);
    if (index<1)|(index>no_total_states)
        state   = -ones(1,length(vec_minsize));
        return
    end
%---Conversion from index to state
    state   = zeros(1,length(vec_minsize));
    K       = index;
    for i=length(vec_minsize):-1:1
        C   = 1;
        for j=1:i-1
            C   = C*(vec_maxsize(j)-vec_minsize(j)+1);
        end
        % state(i)    = idivide(int64(K),int64((C+vec_minsize(i))));
        state(i)    = idivide(int64(K),int64((C)));
        % if (mod(K,(C+vec_minsize(i)))==0)
        if (mod(K,(C))==0)
            state(i)    = state(i)-1;
        end
        % K           = K-(state(i)-vec_minsize(i))*C;
        K           = K-(state(i))*C;
    end
    state   = state+vec_minsize;
end

function [rs,props] = reachable_states(state,s,c,varargin)
%---Given a state, those states reachable in a single step are returned
    [pd,sd] = size(s);
    sd      = sd/2;
%---Check that size(state)==[sd,1]?
    rs      = -ones(pd,sd);
    props   = zeros(pd,1);
    sv      = s(:,sd+1:2*sd)-s(:,1:sd);
    for i=1:pd
        if (length(find(state>=s(i,1:sd)))==sd)
            if isa(c,'function_handle')
                props(i)        = c(i,state,varargin{:});
            else
                props(i)        = c(i);
                for j=1:sd
                    props(i)    = props(i)*nchoosek(state(j),s(i,j));
                end
            end
            rs(i,:) = state+sv(i,:);
        end
    end
end
