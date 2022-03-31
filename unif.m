function [w,W,d,m,t_amen_elapsed] = unif(A,v,t_step,nstep)
    global tol_unif tol_amen

    t_amen_elapsed = 0;
%   Alpha = bound on the sum of propensities
    alpha   = tt_max(-diag(A));
    I       = tt_eye(size(A));
    P       = I + (1/alpha)*A;
%   Taylor degree
    beta    = exp(-alpha*t_step);
    p2      = beta;
    p1      = 1 - p2;
    m       = 1;
    while (p1>(tol_unif/nstep))
        p2  = p2*alpha*t_step/m;
        p1  = p1 - p2;
        m   = m + 1;
    end
    m       = 2^(ceil(log2(m)))-1;
%   Left hand side for the linear system in QTT format
    msize   = [2*ones(ceil(log2(m+1)),1),2*ones(ceil(log2(m+1)),1)];
    vsize   = 2*ones(ceil(log2(m+1)),1);
    St      = tt_matrix(diag(ones(m,1),-1));
    St      = tt_reshape(St,msize);
    bigI    = tkron(I,tt_eye(vsize));
    LHS     = bigI - tkron(P,St);
%   Coefficients for Taylor polynomial
    tcoeff      = zeros(m+1,1);
    tcoeff(1)   = exp(-alpha*t_step);
    for i=2:m+1
        tcoeff(i)   = tcoeff(i-1)*alpha*t_step/(i-1);
    end
    tcoeff  = tt_reshape(tt_tensor(tcoeff),vsize);
%   Right hand side for the linear system in QTT format
    e1      = zeros(m+1,1);
    e1(1)   = 1;
    e1      = tt_tensor(e1);
    e1      = tt_reshape(e1,vsize);
    tones   = tt_ones(vsize);
    d       = length(size(v));
    d1      = d + length(size(e1));

    RHS     = tkron(v,e1);
%   Initial guess for AMEn
    W0      = tkron(v,tones);
%   Solve the Uniformization equation using AMEn
    t_amen  = tic;

    W       = amen_solve2(LHS,RHS,tol_amen,'x0',W0,'verb',0);

    t_amen_elapsed      = t_amen_elapsed + toc(t_amen);
%   Collapse W to find w
    w       = dot(tcoeff,W,d+1,d1);
end
