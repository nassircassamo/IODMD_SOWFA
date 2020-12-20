function [A,B,C,D,K] = px2abcdk(x,u,y,mu,f,p,c,pind)
%PX2ABCDK Estimates the matrices A, B, C and D of the state space model
%  [A,B,C,D,K]=dx2abcdk(x,u,y,f,p,c,pind) estimates the matrices A, B, C,
%  and D of the state space model:
%
%       x(k+1) = A kron(mu(k),x(k)) + B kron(mu(k),u(k)) + K(j) e(k)
%       y(k)   = C kron(mu(k),x(k)) + D kron(mu(k),u(k)) + e(k)
%
%  using the knowledge of the state vector x, the input vector u, the
%  output vector u, the scheduling vector mu. The past window size p is
%  recomended to be higher then the expected system order n. Future window
%  size f must equal or smaller then past window size p. The input pind
%  contains the period size and repeating index. The outputs are the linear
%  parameter-varying matrices A, B, C and D, where matrices have the form
%  of A=[A(1) A(2) ... A(m)].
%
%  [A,B,C,D,K]=dx2abcdk(x,u,y,f,p,c,pind,c) specifies which of the system
%  matrices are constant and not parameter-varing. For each of the matrices
%  A, B, C, and D an 1 or 0 can be given in the vector c. If the A matrix
%  is assumed constant the C matrix cannot be assumed constant and also
%  visa versa. (default C=[0 0 0 0 0])
%
%  References:
%    [1] van Wingerden, J.W., Houtzager, I., Verhaegen, M., Closed-loop
%    identification of the time-varying dynamics of variable-speed wind
%    turbines, Int. J. Robust Nonlinear Control 2008
%
%  See also: pvordarx, pmodx.m, px2abck.m.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology
%  The Netherlands, 2010

% check number if input arguments
if nargin < 8
    error('PX2ABCD requires at least eight input arguments.')
end

% check dimensions of inputs
if size(y,2) < size(y,1)
    y = y';
end
if size(mu,2) < size(mu,1)
    mu = mu';
end
N = size(y,2);
l = size(y,1);
n = size(x{1},1);
s = size(mu,1);
if l == 0
    error('PX2ABCDK requires an output vector y.')
end
if n == 0
    error('PX2ABCDK requires an state vector x.')
end
if ~isequal(N,length(u))
    error('The number of rows of vectors/matrices u and y must be the same.')
end
if isempty(u);
    r = 0;
    u = zeros(0,N);
else
    if size(u,2) < size(u,1)
        u = u';
    end
    r = size(u,1);
    if ~isequal(N,length(u))
        error('The number of rows of vectors/matrices u and y must be the same.')
    end
end

% check the size of the windows
if f > p
    error('Future window size f must equal or smaller then past window p. (f <= p)')
end

% allocate matrices
if c(3) == 1 && c(4) == 1
    CD = zeros(l,n+r);
elseif c(3) == 1 && c(4) == 0
    CD = zeros(l,n+s*r);
elseif c(3) == 0 && c(4) == 1
    CD = zeros(l,s*n+r);
else
    CD = zeros(l,s*n+s*r);
end
if c(1) >= 1 && c(2) == 1 && c(5) == 1
    ABK = zeros(n,n+r+l);
elseif c(1) >= 1 && c(2) == 0 && c(5) == 0
    ABK = zeros(n,n+s*r+s*l);
elseif c(1) >= 1 && c(2) == 1 && c(5) == 0
    ABK = zeros(n,n+r+s*l);
elseif c(1) >= 1 && c(2) == 0 && c(5) == 1
    ABK = zeros(n,n+s*r+l);
elseif c(1) == 0 && c(2) == 0 && c(5) == 0
    ABK = zeros(n,s*n+s*r+s*l);
elseif c(1) == 0 && c(2) == 1 && c(5) == 0
    ABK = zeros(n,s*n+r+s*l);
elseif c(1) == 0 && c(2) == 0 && c(5) == 1
    ABK = zeros(n,s*n+s*r+l);
end

% check the number of repetitions (must be larger then past window)
q = 0;
ind = [];
xx = zeros(n,N);
for k = 1:size(pind,1)
    if length(pind{k,1}) >= p
        q = q + 1;
        ind = [ind pind{k,1}];
        xx(:,p+pind{k,1}) = x{k,1};
    end
end
ind = sort(ind);

% recluster periods
k = 1;
add = 0;
vind = {};
for i = 1:length(ind)
    if i == 1 || add == 1
        vind{k,1} = ind(i);
        add = 0;
    elseif ind(i) == ind(i-1)+1
        vind{k,1} = [vind{k,1} ind(i)];
    else
        k = k + 1;
        add = 1;
    end
end

% start for periodic identification
for k = 1:size(vind,1)
    
    % store data in new matrices
    X = xx(:,p+vind{k,1});
    U = u(:,p+vind{k,1});
    Y = y(:,p+vind{k,1});
    MU = mu(:,p+vind{k,1});
    MUX = khatrirao(MU,X);
    MUU = khatrirao(MU,U);
    
    % obtaining the LPV matrices C and D
    if c(3) == 1 && c(4) == 1
        CD = CD + (Y - CD*vertcat(X,U))*pinv(vertcat(X,U));
        E = Y - CD*vertcat(X,U);
    elseif c(3) == 1 && c(4) == 0
        CD = CD + (Y - CD*vertcat(X,MUU))*pinv(vertcat(X,MUU));
        E = Y - CD*vertcat(X,MUU);
    elseif c(3) == 0 && c(4) == 1
        CD = CD + (Y - CD*vertcat(MUX,U))*pinv(vertcat(MUX,U));
        E = Y - CD*vertcat(MUX,U);
    else
        CD = CD + (Y - CD*vertcat(MUX,MUU))*pinv(vertcat(MUX,MUU));
        E = Y - CD*vertcat(MUX,MUU);
    end
    if c(5) == 0
        MUE = khatrirao(MU,E);
    end
    
    % obtaining the LPV matrices A, B and K
    if c(1) >= 1 && c(2) == 1 && c(5) == 1
        ABK = ABK + (X(:,2:end) - ABK*vertcat(X(:,1:end-1),U(:,1:end-1),E(:,1:end-1)))*pinv(vertcat(X(:,1:end-1),U(:,1:end-1),E(:,1:end-1)));
    elseif c(1) >= 1 && c(2) == 0 && c(5) == 0
        ABK = ABK + (X(:,2:end) - ABK*vertcat(X(:,1:end-1),MUU(:,1:end-1),MUE(:,1:end-1)))*pinv(vertcat(X(:,1:end-1),MUU(:,1:end-1),MUE(:,1:end-1)));
    elseif c(1) >= 1 && c(2) == 1 && c(5) == 0
        ABK = ABK + (X(:,2:end) - ABK*vertcat(X(:,1:end-1),U(:,1:end-1),MUE(:,1:end-1)))*pinv(vertcat(X(:,1:end-1),U(:,1:end-1),MUE(:,1:end-1)));
    elseif c(1) >= 1 && c(2) == 0 && c(5) == 1
        ABK = ABK + (X(:,2:end) - ABK*vertcat(X(:,1:end-1),MUU(:,1:end-1),E(:,1:end-1)))*pinv(vertcat(X(:,1:end-1),MUU(:,1:end-1),E(:,1:end-1)));
    elseif c(1) == 0 && c(2) == 0 && c(5) == 0
        ABK = ABK + (X(:,2:end) - ABK*vertcat(MUX(:,1:end-1),MUU(:,1:end-1),MUE(:,1:end-1)))*pinv(vertcat(MUX(:,1:end-1),MUU(:,1:end-1),MUE(:,1:end-1)));
    elseif c(1) == 0 && c(2) == 1 && c(5) == 0
        ABK = ABK + (X(:,2:end) - ABK*vertcat(MUX(:,1:end-1),U(:,1:end-1),MUE(:,1:end-1)))*pinv(vertcat(MUX(:,1:end-1),U(:,1:end-1),MUE(:,1:end-1)));
    elseif c(1) == 0 && c(2) == 0 && c(5) == 1
        ABK = ABK + (X(:,2:end) - ABK*vertcat(MUX(:,1:end-1),MUU(:,1:end-1),E(:,1:end-1)))*pinv(vertcat(MUX(:,1:end-1),MUU(:,1:end-1),E(:,1:end-1)));
    end
end

% obtaining the LPV matrices C and D
if c(3) == 1 && c(4) == 1
    C = [CD(:,1:n) zeros(l,(s-1)*n)];
    D = [CD(:,n+1:n+r) zeros(l,(s-1)*r)];
elseif c(3) == 1 && c(4) == 0
    C = [CD(:,1:n) zeros(l,(s-1)*n)];
    D = CD(:,n+1:n+s*r);
elseif c(3) == 0 && c(4) == 1
    C = CD(:,1:s*n);
    D = [CD(:,s*n+1:s*n+r) zeros(l,(s-1)*r)];
else
    C = CD(:,1:s*n);
    D = CD(:,s*n+1:s*n+s*r);
end

% obtaining the LPV matrices A, B and K
if c(1) >= 1 && c(2) == 1 && c(5) == 1
    A = [ABK(:,1:n) zeros(n,(s-1)*n)];
    B = [ABK(:,n+1:n+r) zeros(n,(s-1)*r)];
    K = [ABK(:,n+r+1:n+r+l) zeros(n,(s-1)*l)];
elseif c(1) >= 1 && c(2) == 0 && c(5) == 0
    A = [ABK(:,1:n) zeros(n,(s-1)*n)];
    B = ABK(:,n+1:n+s*r);
    K = ABK(:,n+s*r+1:n+s*r+s*l);
elseif c(1) >= 1 && c(2) == 1 && c(5) == 0
    A = [ABK(:,1:n) zeros(n,(s-1)*n)];
    B = [ABK(:,n+1:n+r) zeros(n,(s-1)*r)];
    K = ABK(:,n+r+1:n+r+s*l);
elseif c(1) >= 1 && c(2) == 0 && c(5) == 1
    A = [ABK(:,1:n) zeros(n,(s-1)*n)];
    B = ABK(:,n+1:n+s*r);
    K = [ABK(:,n+s*r+1:n+s*r+l) zeros(n,(s-1)*l)];
elseif c(1) == 0 && c(2) == 0 && c(5) == 0
    A = ABK(:,1:s*n);
    B = ABK(:,s*n+1:s*n+s*r);
    K = ABK(:,s*n+s*r+1:s*n+s*r+s*l);
elseif c(1) == 0 && c(2) == 1 && c(5) == 0
    A = ABK(:,1:s*n);
    B = [ABK(:,s*n+1:s*n+r) zeros(n,(s-1)*r)];
    K = ABK(:,s*n+r+1:s*n+r+s*l);
elseif c(1) == 0 && c(2) == 0 && c(5) == 1
    A = ABK(:,1:s*n);
    B = ABK(:,s*n+1:s*n+s*r);
    K = [ABK(:,s*n+s*r+1:s*n+s*r+l) zeros(n,(s-1)*l)];
end