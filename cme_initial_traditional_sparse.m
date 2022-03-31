function p_0 = cme_initial_traditional_sparse(state_space_list,initial_state)
%   Find the position of the initial state in the state space list
    pos = -1;
    for row=1:size(state_space_list,1)
        if isequal(initial_state,state_space_list(row,:))
            pos = row;
            break
        end
    end
%   Initial probability vector will be 0 everywhere and 1 at the initial state
    if pos>0
        p_0 = zeros(size(state_space_list,1),1);
        p_0(pos)    = 1;
    end
end
