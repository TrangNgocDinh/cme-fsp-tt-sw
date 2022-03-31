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

function [w, err, hump] = cme_solver_expokit(t,A,v,tol,m)
%---Find size of the matrix
    [n,n] = size(A);
%---Set up default tolerance and Krylov degree, unless given by user
    if nargin == 3,
        tol = 1.0e-7;
        m   = min(n,30);
    end;
    if nargin == 4,
        m   = min(n,30);
    end;
%---Initialize internal parameters
    anorm   = norm(A,'inf');
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
    normv   = norm(v);
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
%=============================Time integration to compute w = exp(t*A)*v
    while t_now<t_out
        nstep   = nstep + 1;
        t_step  = min( t_out-t_now,t_new );
        V       = zeros(n,m+1);
        H       = zeros(m+2,m+2);
%-------Arnoldi process to find V & H
%       First column in V is the normalized v
        V(:,1)  = (1/beta)*w;
%       For every next column...
        for j=1:m
%           Initialize next column = A * previous column
            p   = A*V(:,j);
%           For every past column...
            for i=1:j
%               h_ij = inner product of this past and the next column
                H(i,j)  = V(:,i)'*p;
%               Orthogonize p against this past column
                p       = p-H(i,j)*V(:,i);
            end;
%           HAPPY BREAKDOWN if norm of the next column is less than the
%           tolerance, in which case the next time point is the final
%           time point
            s   = norm(p);
            if s<btol,
                k1      = 0;
                mb      = j;
                t_step  = t_out-t_now;
                break;
            end;
%           Finalize next column in V & H
            H(j+1,j)    = s;
            V(:,j+1)    = (1/s)*p;
        end;
%       Extend the H matrix if happy breakdown did not occur
        if k1 ~= 0,
            H(m+2,m+1)  = 1;
            avnorm      = norm(A*V(:,m+1));
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
        w       = V(:,1:mx)*(beta*F(1:mx,1));
%       Update the hump
        beta    = norm(w);
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
