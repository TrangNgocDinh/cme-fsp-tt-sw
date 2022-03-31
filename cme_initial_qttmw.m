function v = cme_initial_qttmw(lb,l2size)
    global initial_state no_species
    vtmp                        = zeros(2^l2size(1),1);
    vtmp(initial_state(1)-lb(1)+1)    = 1;
    v                           = tt_tensor(vtmp);
    for i=2:no_species
        vtmp                        = zeros(2^l2size(i),1);
        vtmp(initial_state(i)-lb(i)+1)    = 1;
        v                           = tkron(v,tt_tensor(vtmp));
    end
    v   = reshape(v,2*ones(1,sum(l2size)));
end
