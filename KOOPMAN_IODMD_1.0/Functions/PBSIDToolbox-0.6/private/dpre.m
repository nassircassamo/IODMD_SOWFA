function [X,K] = dpre(A,B,Q,R,S,E,tol,maxit)
%DPRE Discrete-time Periodic Riccati Equation
%  [X,K]=DPRE(A,B,Q,R,S,E) computes the unique stabilizing solution X{k},
%  k = 1:P, of the discrete-time periodic Riccati equation
%
%   E{k}'X{k}E{k} = A{k}'X{k+1}A{k} - (A{k}'X{k+1}B{k} + S{k})*...
%                 (B{k}'X{k+1}B{k} + R{k})\(A{k}'X{k+1}B{k} + S{k})' + Q{k}
%
%  When omitted, R, S and E are set to the default values R{k}=I, S{k}=0,
%  and E{k}=I. Beside the solution X{k}, DPRE also returns the gain matrix
%
%   K{k} = (B{k}'X{k+1}B{k} + R{k})\(B{k}'X{k+1}A{k} + S{k}'),
%
%  All input matrices have to be multidimensional arrays, like matrix 
%  A(N,N,P) and B(N,R,P). Output matrices are also multidimensional arrays
%  in the size of X(N,N,P) and K(R,N,P).
%  
%  [X,K]=DPRE(A,B,Q,R,S,E,TOL) specifies the tolerance of the cyclic qz 
%  method. If tol is [] then DPRE uses the default, 1e-6.
%
%  [X,K]=DPRE(A,B,Q,R,S,E,TOL,MAXIT) specifies the maximum number of 
%  iterations. If MAXIT is [] then DPRE uses the default, 1000. Warning is
%  given when the maximum iterations is reached.
%
%  See also DARE.

%  This version uses a cyclic qz method, see references.

%  References:
%    [1] J.J. Hench and A.J. Laub, Numerical solution of the discrete-time 
%        periodic Riccati equation, IEEE Trans. on automatic control, 1994 
%    [2] Varga, A., On solving discrete-time periodic Riccati equations,
%        Proc. of 16th IFAC World Congress 2005, Prague, July 2005.
 
%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% assign default values to unspecified parameters
[m,n,p] = size(A);
[mb,r,pb] = size(B);
if (nargin < 8) || isempty(maxit)
    maxit = 1000;
end
if (nargin < 6) || isempty(E)
    E = zeros(m,n,p);
    for i = 1:p
        E(:,:,i) = eye(m,n);
    end
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
if (nargin < 7) || isempty(tol)
    tol = 1e-6;
end

% check input arguments
if nargin < 2 
    error('DPRE requires at least three input arguments')
end
[mq,nq,pq] = size(Q);
[mr,nr,pr] = size(R);
[ms,ns,ps] = size(S);
[me,ne,pe] = size(E);
if ~isequal(p,pb,pq,pr,ps,pe)
    error('The number of periods must be the same for A, B, Q, R, S and E.')
end
if ~isequal(m,me,mq,mb)
    error('The number of rows of matrix A, B, E and Q must be the same size.')
end
if ~isequal(n,ne,nq)
    error('The number of columns of matrix A, E and Q must be the same size.')
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
M = zeros(2*n+r,2*n+r,p);
L = zeros(2*n+r,2*n+r,p);
V = zeros(2*n+r,2*n+r,p);
Z = zeros(2*n+r,2*n+r,p);
Y = zeros(2*n+r,2*n+r,p);
T = zeros(n+r,n,p);

% build the periodic matrix pairs
for i = 1:p
    L(:,:,i) = [A(:,:,i) zeros(n) B(:,:,i);
        -Q(:,:,i) E(:,:,i) -S(:,:,i);
        S(:,:,i)' zeros(r,n) R(:,:,i)];
    M(:,:,i) = [E(:,:,i) zeros(n,n+r);
        zeros(n) A(:,:,i)' zeros(n,r);
        zeros(r,n) -B(:,:,i)' zeros(r)];
    V(:,:,i) = eye(2*n+r);
    Z(:,:,i) = eye(2*n+r);
end

% cyclic qz decomposition
k = 1;
ok = true;
res = ones(1,p);
while ok == true && k <= maxit
    for j = 1:p    
        % QR decomposition
        [Q,R] = qr(M(:,:,j));
        V(:,:,j) = V(:,:,j)*Q';
        M(:,:,j) = Q'*M(:,:,j);
        
        % RQ decomposition
        [Q,R] = qr(fliplr((Q'*L(:,:,j))'));
        Q = fliplr(Q);
        R = fliplr(flipud(R))';

        Z(:,:,j) = Z(:,:,j)*Q;
        Y(:,:,j) = Q;
        L(:,:,j) = R;
    end
    
    for j = 1:p
        if j == p
            M(:,:,p) = M(:,:,p)*Y(:,:,1);
        else
            M(:,:,j) = M(:,:,j)*Y(:,:,j+1);
        end

        T1 = Z(n+1:2*n+r,1:n,j)/Z(1:n,1:n,j);
            
        % calculate residue
        res(j) = norm(T1 - T(:,:,j)); 
        T(:,:,j) = T1;
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

% retrieve K and X matrices
X = zeros(n,n,p);
K = zeros(r,n,p);
for i = 1:p
    X(:,:,i) = T(1:n,1:n,i);
    K(:,:,i) = -T(n+1:n+r,1:n,i);
end





