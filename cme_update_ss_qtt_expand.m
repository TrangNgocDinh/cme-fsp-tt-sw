function v1 = cme_update_ss_qtt_expand(v, lfspsize, expand_ind)
% Padd the QTT-formatted tensor v with zeros when the FSP is doubled in ALL
% directions
% expand_ind(1:N) is a logical vector that specify the species to be
% expanded

    v1 = core(v);
    N = length(lfspsize);
    jc = cumsum(lfspsize);  % indices of cores corresponding to the edges of the hyper-rectangle
    for i = 1:N
        if (expand_ind(i))
            v1{jc(i)} = [v1{jc(i)};zeros(size( v1{jc(i)}) )];
        end
    end
    v1 = tt_tensor(v1);
    v1 = tt_reshape(v1, 2*ones(sum(lfspsize) + sum(expand_ind),1));
end
