function A = cme_generator_qttmw(lb,l2size,propen_func_partial)
    global no_species no_reactions stoich_mat_total tol_qtt_mat FSP_qtt_size
%   Do we need lower_bond = lb and upper_bound = ub as global variables ?
%   Do we need update FSP_qtt_size = l2size ?

    ub  = lb + 2.^l2size -1; % Set the corresponding upper bound

    Sk  = cell(no_species,1);

    for s=1:no_reactions
%       Define the shift operator
        for k=1:no_species
            Sk{k}   = tt_matrix(tt_qshift(l2size(k),stoich_mat_total(s,k)));
        end
        Shift   = Sk{1};
        for k=2:no_species
%           Unlike kron, tkron will preserve the multi-index...
            Shift   = tkron(Shift,Sk{k});
        end
%       Find the propensity tensor
        for k=1:no_species
            Sk{k}   = propen_func_partial(k,l2size(k),s,lb(k),ub(k));
            Sk{k}   = tt_reshape(Sk{k},2*ones(l2size(k),1));
        end
        omega=Sk{1};
        for k=2:no_species
%           Unlike kron, tkron will preserve the multi-index...
            omega   = tkron(omega,Sk{k});
        end
        Momega      = diag(omega);
        if s==1
            A       = Shift*Momega-Momega;
        else
            A       = A+Shift*Momega-Momega;
        end;
    end
    A   = round(A,tol_qtt_mat);
%   Do we need to update the global variable lower_bound and upper_bound ?
end
