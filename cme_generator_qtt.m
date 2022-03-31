function A = cme_generator_qtt(propen_func_partial)
    global no_species no_reactions stoich_mat_total FSP_qtt_size tol_qtt_mat

    Sk  = cell(no_species,1);

    for s=1:no_reactions
%       Define the shift operator
        for k=1:no_species
            Sk{k}   = tt_matrix(tt_qshift(FSP_qtt_size(k),stoich_mat_total(s,k)));
        end
        Shift   = Sk{1};
        for k=2:no_species
%           Unlike kron, tkron will preserve the multi-index...
            Shift   = tkron(Shift,Sk{k});
        end
%       Find the propensity tensor
        for k=1:no_species
            Sk{k}   = propen_func_partial(k,FSP_qtt_size(k),s);
            Sk{k}   = tt_reshape(Sk{k},2*ones(FSP_qtt_size(k),1));
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
end
