%  [w, err, hump] = expv( t, A, v, tol, m )
%  EXPV computes an approximation of w = exp(t*A)*v for a
%  general matrix A using Krylov subspace  projection techniques.
%  It does not compute the matrix exponential in isolation but instead,
%  it computes directly the action of the exponential operator on the
%  operand vector. This way of doing so allows for addressing large
%  sparse problems. The matrix under consideration interacts only
%  via matrix-vector products (matrix-free method).
%
%  w = expv( t, A, v )
%  computes w = exp(t*A)*v using a default tol = 1.0e-7 and m = 30.
%
%  [w, err] = expv( t, A, v )
%  renders an estimate of the error on the approximation.
%
%  [w, err] = expv( t, A, v, tol )
%  overrides default tolerance.
%
%  [w, err, hump] = expv( t, A, v, tol, m )
%  overrides default tolerance and dimension of the Krylov subspace,
%  and renders an approximation of the `hump'.
%
%  The hump is defined as:
%          hump = max||exp(sA)||, s in [0,t]  (or s in [t,0] if t < 0).
%  It is used as a measure of the conditioning of the matrix exponential
%  problem. The matrix exponential is well-conditioned if hump = 1,
%  whereas it is poorly-conditioned if hump >> 1. However the solution
%  can still be relatively fairly accurate even when the hump is large
%  (the hump is an upper bound), especially when the hump and
%  ||w(t)||/||v|| are of the same order of magnitude (further details in
%  reference below).
%
%  Example 1:
%  ----------
%    n = 100;
%    A = rand(n);
%    v = eye(n,1);
%    w = expv(1,A,v);
%
%  Example 2:
%  ----------
%    % generate a random sparse matrix
%    n = 100;
%    A = rand(n);
%    for j = 1:n
%        for i = 1:n
%            if rand < 0.5, A(i,j) = 0; end;
%        end;
%    end;
%    v = eye(n,1);
%    A = sparse(A); % invaluable for a large and sparse matrix.
%
%    tic
%    [w,err] = expv(1,A,v);
%    toc
%
%    disp('w(1:10) ='); disp(w(1:10));
%    disp('err =');     disp(err);
%
%    tic
%    w_matlab = expm(full(A))*v;
%    toc
%
%    disp('w_matlab(1:10) ='); disp(w_matlab(1:10));
%    gap = norm(w-w_matlab)/norm(w_matlab);
%    disp('||w-w_matlab|| / ||w_matlab|| ='); disp(gap);
%
%  In the above example, n could have been set to a larger value,
%  but the computation of w_matlab will be too long (feel free to
%  discard this computation).
%
%  See also MEXPV, EXPOKIT.

%  Roger B. Sidje (rbs@maths.uq.edu.au)
%  EXPOKIT: Software Package for Computing Matrix Exponentials.
%  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
function [w,err,hump] = tt_expv(t,A,v,tol,m)
    global tol_tt_mvk4

    T_total = tic;
    t_ttmvk4_elapsed    = 0;
    t_norm_elapsed      = 0;
    t_dot_elapsed       = 0;
    t_round_elapsed     = 0;

%---Find size of the matrix
    % [n,n] = size(A);
    %---T: A.n: row size of A in tt-format
    n = prod(A.n);
%---Set up default tolerance and Krylov degree, unless given by user
    if nargin == 3,
        tol = 1.0e-7;
        m   = min(n,30);
    end;
    if nargin == 4,
        m   = min(n,30);
    end;
%---Initialize internal parameters
    t_norm  = tic;

    anorm   = norm(A,2);

    t_norm_elapsed  = t_norm_elapsed+toc(t_norm);

    mxrej   = 10;
    btol    = 1.0e-7;
    gamma   = 0.9;
    delta   = 1.2;
    mb      = m;
    t_out   = abs(t);
    nstep   = 0;
    t_new   = 0;
    t_now   = 0;
    s_error = 0;
    rndoff  = anorm*eps;

    k1      = 2;
    xm      = 1/m;

    t_norm  = tic;

    normv   = norm(v);

    t_norm_elapsed  = t_norm_elapsed+toc(t_norm);

    beta    = normv;
    fact    = (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1));
    t_new   = (1/anorm)*((fact*tol)/(4*beta*anorm))^xm;
    s       = 10^(floor(log10(t_new))-1);
    t_new   = ceil(t_new/s)*s;
    sgn     = sign(t);
    nstep   = 0;
%---Initialize w to be input vector
    w       = v;
%---Initialize hump to be ||v||
    hump    = normv;
%---
    max_r   = 2^(floor(log2(n)/2));
%=============================Time integration to compute w = exp(t*A)*v
    while t_now<t_out
        nstep   = nstep + 1;
        t_step  = min(t_out-t_now,t_new);

        V       = cell(m+1,1);

        H       = zeros(m+2,m+2);
%-------Arnoldi process to find V & H
%       First column in V is the normalized v
        V{1}    = (1/beta)*w;
%       For every next column...
        for j=1:m
%           Initialize next column = A * previous column
            t_ttmvk4    = tic;

            A_d         = size(A.m,1);
            A_N         = 2*ones(A_d,1);
            A_r         = 2*ones(A_d+1,1);
            A_r(1)      = 1;
            A_r(end)    = 1;
            A_y0        = tt_rand(A_N,A_d,A_r);
            mtimes_parameters = {'nswp',10,'rmax',max_r,...
                                 'kicktype','mr','pcatype', 'svd',...
                                 'y0',A_y0,...
                                 'verb',0,'kickrank',10,'step_dpow',0.1,...
                                 'min_dpow',1,'bot_conv',0.1,...
                                 'top_conv',0.99,'block_order',[],...
                                 'als_tol_low',5,'als_iters',4};
            p = tt_mvk4(A, V{j},tol_tt_mvk4,mtimes_parameters{:});

            t_ttmvk4_elapsed    = t_ttmvk4_elapsed+toc(t_ttmvk4);
%           For every past column...
            for i=1:j
%               h_ij = inner product of this past and the next column
                t_dot   = tic;

                H(i,j)  = dot(V{i},p);

                t_dot_elapsed   = t_dot_elapsed+toc(t_dot);
%               Orthogonize p against this past column
                p       = p - H(i,j)*V{i};

                t_round = tic;

                p       = round(p, 1e-14);

                t_round_elapsed = t_round_elapsed+toc(t_round);
            end;
%           HAPPY BREAKDOWN if norm of the next column is less than the
%           tolerance, in which case the next time point is the final
%           time point
            t_norm  = tic;

            s   = norm(p);

            t_norm_elapsed  = t_norm_elapsed+toc(t_norm);

            if s<btol,
                k1      = 0;
                mb      = j;
                t_step  = t_out-t_now;
                break;
            end;
%           Finalize next column in V & H
            H(j+1,j)    = s;
            V{j+1} = (1/s)*p;
        end;
%       Extend the H matrix if happy breakdown did not occur
        if k1 ~= 0,
            H(m+2,m+1)  = 1;

            t_ttmvk4    = tic;

            av          = tt_mvk4(A,V{m+1},tol_tt_mvk4,mtimes_parameters{:});

            t_ttmvk4_elapsed    = t_ttmvk4_elapsed+toc(t_ttmvk4);

            t_norm      = tic;

            avnorm      = norm(av);

            t_norm_elapsed  = t_norm_elapsed+toc(t_norm);
        end;
%-------Estimate error and find time step
        ireject = 0;
        while ireject<=mxrej,
%           Compute F = exp(tau*H) for the actual Krylov degree (+2)
            mx  = mb+k1;
            F   = expm(sgn*t_step*H(1:mx,1:mx));
%           Find the local error estimation
            if k1==0,
    	        err_loc    = btol;
                break;
            else
                phi1    = abs( beta*F(m+1,1) );
                phi2    = abs( beta*F(m+2,1) * avnorm );
                if phi1>10*phi2,
                    err_loc = phi2;
                    xm      = 1/m;
                elseif phi1>phi2,
                    err_loc = (phi1*phi2)/(phi1-phi2);
                    xm      = 1/m;
                else
                    err_loc = phi1;
                    xm      = 1/(m-1);
                end;
            end;
%           Check if the error is below tolerance...
            if err_loc<=delta*t_step*tol,
%               If yes, then move on
                break;
            else
%               If not, then decrease the time step and try again
                t_step      = gamma*t_step*(t_step*tol/err_loc)^xm;
                s           = 10^(floor(log10(t_step))-1);
                t_step      = ceil(t_step/s)*s;
%               If too many rejections, then give up
                if ireject==mxrej,
                    error('The requested tolerance is too high.');
                end;
                ireject     = ireject + 1;
            end;
        end;
%       Find the Krylov approximation for this time step
        mx      = mb+max(0,k1-1);
        w = F(1,1)*V{1};
        for i = 2:mx
            w   = w + F(i,1)*V{i};

            t_round         = tic;

            w   = round(w, 1e-14);

            t_round_elapsed = t_round_elapsed+toc(t_round);
        end
        w   = beta*w;
%       Update the hump
        t_norm  = tic;

        beta    = norm(w);

        t_norm_elapsed  = t_norm_elapsed+toc(t_norm);

        hump    = max(hump,beta);
%       Update the time and guess the next time step
        t_now   = t_now + t_step;
        t_new   = gamma * t_step * (t_step*tol/err_loc)^xm;
        s       = 10^(floor(log10(t_new))-1);
        t_new   = ceil(t_new/s) * s;
%       Estimate global error as sum of local errors
        err_loc = max(err_loc,rndoff);
        s_error = s_error + err_loc;
    end;
    err     = s_error;
    hump    = hump / normv;

    T_total_elapsed = toc(T_total);
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
    fprintf('   Total time for TT-EXPOKIT: %f\n',T_total_elapsed);
    fprintf('   Percentage of runtime for tt-mvk4 is %d%%\n',round(100*t_ttmvk4_elapsed/T_total_elapsed));
    fprintf('   Percentage of runtime for norm is %d%%\n',round(100*t_norm_elapsed/T_total_elapsed));
    fprintf('   Percentage of runtime for dot is %d%%\n',round(100*t_dot_elapsed/T_total_elapsed));
    fprintf('   Percentage of runtime for round is %d%%\n',round(100*t_round_elapsed/T_total_elapsed));

    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
