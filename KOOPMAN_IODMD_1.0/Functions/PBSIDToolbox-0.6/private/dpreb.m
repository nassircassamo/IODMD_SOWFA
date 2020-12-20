function [X,K] = dpreb(A,B,Q,R,S,tol,maxit)
%DPRE Discrete-time Periodic Riccati Equation
%  [X,K]=DPRE(A,B,Q,R,S) computes the unique stabilizing solution X{k},
%  k = 1:P, of the discrete-time periodic Riccati equation
%
%   X{k} = A{k}'X{k+1}A{k} - (A{k}'X{k+1}B{k} + S{k})*...
%                 (B{k}'X{k+1}B{k} + R{k})\(A{k}'X{k+1}B{k} + S{k})' + Q{k}
%
%  When omitted, R, S and E are set to the default values R{k}=I, S{k}=0,
%  and E{k}=I. Beside the solution X{k}, PDRE also returns the gain matrix
%
%   K{k} = (B{k}'X{k+1}B{k} + R{k})\(B{k}'X{k+1}A{k} + S{k}'),
%
%  All input matrices have to be multidimensional arrays, like matrix 
%  A(N,N,P) and B(N,R,P). Output matrices are also multidimensional arrays
%  in the size of X(N,N,P) and K(R,N,P).
%  
%  [X,K]=DPRE(A,B,Q,R,S,TOL) specifies the tolerance of the backward
%  iteration method. If tol is [] then PPRE uses the default, 1e-6.
%
%  [X,K]=DPRE(A,B,Q,R,S,TOL,MAXIT) specifies the maximum number of 
%  iterations. If MAXIT is [] then DPRE uses the default, 1000. Warning is
%  given when the maximum iterations is reached.
%
%  See also DARE.

%  This version uses a backward iteration 

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% assign default values to unspecified parameters
[m,n,p] = size(A);
[mb,r,pb] = size(B);
if (nargin < 7) || isempty(maxit)
    maxit = 1000;
end
if (nargin < 5) || isempty(S)
    S = zeros(mb,r,pb);
end
if (nargin < 4) || isempty(R)
    R = zeros(r,r,pb);
    for i = 1:pb
        R(:,:,i) = eye(r);
    end
end
if (nargin < 6) || isempty(tol)
    tol = 1e-6;
end

% check input arguments
if nargin < 2 
    error('DPRE requires at least three input arguments')
end
[mq,nq,pq] = size(Q);
[mr,nr,pr] = size(R);
[ms,ns,ps] = size(S);
if ~isequal(p,pb,pq,pr,ps)
    error('The number of periods must be the same for A, B, Q, R and S.')
end
if ~isequal(m,mq,mb)
    error('The number of rows of matrix A, B, and Q must be the same size.')
end
if ~isequal(n,nq)
    error('The number of columns of matrix A, and Q must be the same size.')
end
if ~isequal(mb,ms)
    error('The number of rows of matrix B and S must be the same size.')
end
if ~isequal(r,ns)
    error('The number of columns of matrix B and S must be the same size.')
end
if ~isequal(r,nr)
    error('The number of columns of matrix R must be the same size as the number of columns of matrix B.')
end
if ~isequal(r,mr)
    error('The number of rows of matrix R must be the same size as the number of columns of matrix B.')
end

% allocate matrices
X = zeros(n,n,p); 
K = zeros(r,n,p);

% initial guess X1
X(:,:,1) = dare(A(:,:,1),B(:,:,1),Q(:,:,1),R(:,:,1),S(:,:,1));

% backwards iteration
k = 1;
ok = true;
res = ones(1,p);
while ok == true && k <= maxit
    for i = p:-1:1
        if i == p
            K(:,:,i) = (B(:,:,i)'*X(:,:,1)*B(:,:,i) + R(:,:,i))\(A(:,:,i)'*X(:,:,1)*B(:,:,i) + S(:,:,i))';
            X1 = Q(:,:,i) + A(:,:,i)'*X(:,:,1)*A(:,:,i) - (A(:,:,i)'*X(:,:,1)*B(:,:,i) + S(:,:,i))*K(:,:,i);
        else
            K(:,:,i) = (B(:,:,i)'*X(:,:,i+1)*B(:,:,i) + R(:,:,i))\(A(:,:,i)'*X(:,:,i+1)*B(:,:,i) + S(:,:,i))';
            X1 = Q(:,:,i) + A(:,:,i)'*X(:,:,i+1)*A(:,:,i) - (A(:,:,i)'*X(:,:,i+1)*B(:,:,i) + S(:,:,i))*K(:,:,i);
        end
        
        res(i) = norm(X1 - X(:,:,i));
        X(:,:,i) = X1;
    end
       
    if all(res <= tol)
        ok = false;
    end    

    k = k + 1;
end

% return warning if k exceeds maxit
if ok == true
    warning('DPRE:maxIt','Maximum number of iterations exceeded.')
end


