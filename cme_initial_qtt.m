function v = cme_initial_qtt()
    global FSP_qtt_size initial_state no_species
    vtmp                        = zeros(2^FSP_qtt_size(1),1);
    vtmp(initial_state(1)+1)    = 1;
    v                           = tt_tensor(vtmp);
    for i=2:no_species
        vtmp                        = zeros(2^FSP_qtt_size(i),1);
        vtmp(initial_state(i)+1)    = 1;
        v                           = tkron(v,tt_tensor(vtmp));
    end
    v   = reshape(v,2*ones(1,sum(FSP_qtt_size)));
end
