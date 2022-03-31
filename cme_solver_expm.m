function p_out = cme_solver_expm(t,A,p_in)
    p_out   = expm(t*A)*p_in;
end
