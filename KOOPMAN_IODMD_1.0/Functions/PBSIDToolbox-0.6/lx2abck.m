function [A,B,C,K,gamma] = lx2abck(x,u,y,mu,f,p,c,pred,stable)
%LX2ABCK Estimates the matrices A, B, and C of the LPV state space model
%  [A,B,C,K]=lx2abck(x,u,y,mu,f,p) estimates the matrices A, B, and C
%  of the state space model:
%
%     x(k+1) = A kron(mu(k),x(k)) + B kron(mu(k),u(k)) + K kron(mu(k),e(k))
%     y(k)   = C x(k) + e(k)
%
%  using the knowledge of the state vector x, the input vector u, the
%  output vector u, the scheduling vector mu. The past window size p is
%  recomended to be higher then the expected system order n. Future window
%  size f must equal or smaller then past window size p. The outputs are
%  the linear parameter-varying matrices A, B, K, where matrices have the
%  form of A=[A(1) A(2) ... A(m)].
%
%  [A,B,C,K]=lx2abck(x,u,y,mu,f,p,c) specifies which of the system
%  matrices are constant and not parameter-varing. For each of the matrices
%  A, B, and K an 1 or 0 can be given in the vector c.(default C=[0 0 0])
%
%  [A,B,C,D,K,gamma]=lx2abcdk(x,u,y,mu,f,p,c,pred,stable) with
%
%  pred = 1 (default: pred=0) returns the system matrices At, B, K, C 
%  of the predictor form:
%     x_est(k+1) = At kron(mu(k),x(k)) + B kron(mu(k),u(k)) +
%                 + K kron(mu(k),y(k))
%     y_est(k)   = C x(k)
%   At = [A(1)-K(1)*C, ..., A(m)-K(m)*C]
%
%  stable = 1 (default: stable=0) designs K such that the predictor form is 
%  stable in the range of the scheduling mu, and that the predictor is Hinf
%  optimal w.r.t the noise in the model (see also the function
%  lobshinfsyn). In that case, gamma returns the Hinf norm of the transfer
%  function between noise and state error i.e. a smaller indicates better 
%  estimation of the states given the noise properties. This is only
%  implemented for c(3) = 0.
%
%  References:
%    [1] J.W. van Wingerden, and M. Verhaegen, ``Subspace identification
%    of Bilinear and LPV systems for open- and closed-loop data'',
%    Automatica 45, pp 372--381, 2009.
%
%  See also: lvordarx, lmodx.m, lx2abcdk.m, lobshinfsyn.m.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology
%  The Netherlands, 2010

%  Pieter Gebraad
%  Delft Center of Systems and Control
%  Delft University of Technology
%  The Netherlands, 2011

gamma = [];

% check number of input arguments
if nargin < 6
    error('LX2ABCK requires at least six input arguments.')
end

if nargin < 9 || isempty(stable)
    stable = 0;
end
if nargin < 8 || isempty(pred)
    pred = 0;
end
if nargin < 7 || isempty(c)
    c = [0 0 0];
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
    error('LX2ABCK requires an output vector y.')
end
if n == 0
    error('LX2ABCK requires an state vector x.')
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

% remove the window sizes from input and output vector
U = u(:,p+1:p+size(x,2));
Y = y(:,p+1:p+size(x,2));
MU = mu(:,p+1:p+size(x,2));
MUX = khatrirao(MU,x);
MUU = khatrirao(MU,U);
    
C = Y*pinv(x);
E = Y - C*x;
if c(3) == 0
    MUE = khatrirao(MU,E);
end
C = [C zeros(l,(s-1)*n)];

% obtaining the LPV matrices A, B and K
if c(1) >= 1 && c(2) == 1 && c(3) == 1
    ABK = x(:,2:end)*pinv(vertcat(x(:,1:end-1),U(:,1:end-1),E(:,1:end-1)));
elseif c(1) >= 1 && c(2) == 0 && c(3) == 0
    ABK = x(:,2:end)*pinv(vertcat(x(:,1:end-1),MUU(:,1:end-1),MUE(:,1:end-1)));
elseif c(1) >= 1 && c(2) == 1 && c(3) == 0
    ABK = x(:,2:end)*pinv(vertcat(x(:,1:end-1),U(:,1:end-1),MUE(:,1:end-1)));
elseif c(1) >= 1 && c(2) == 0 && c(3) == 1
    ABK = x(:,2:end)*pinv(vertcat(x(:,1:end-1),MUU(:,1:end-1),E(:,1:end-1)));
elseif c(1) == 0 && c(2) == 0 && c(3) == 0
    ABK = x(:,2:end)*pinv(vertcat(MUX(:,1:end-1),MUU(:,1:end-1),MUE(:,1:end-1)));
elseif c(1) == 0 && c(2) == 1 && c(3) == 0
    ABK = x(:,2:end)*pinv(vertcat(MUX(:,1:end-1),U(:,1:end-1),MUE(:,1:end-1)));
elseif c(1) == 0 && c(2) == 0 && c(3) == 1
    ABK = x(:,2:end)*pinv(vertcat(MUX(:,1:end-1),MUU(:,1:end-1),E(:,1:end-1)));
end

% obtaining the LPV matrices A, B and K
if c(1) >= 1 && c(2) == 1 && c(3) == 1
    A = [ABK(:,1:n) zeros(n,(s-1)*n)];
    B = [ABK(:,n+1:n+r) zeros(n,(s-1)*r)];
    K = [ABK(:,n+r+1:n+r+l) zeros(n,(s-1)*l)];
elseif c(1) >= 1 && c(2) == 0 && c(3) == 0
    A = [ABK(:,1:n) zeros(n,(s-1)*n)];
    B = ABK(:,n+1:n+s*r);
    K = ABK(:,n+s*r+1:n+s*r+s*l);
elseif c(1) >= 1 && c(2) == 1 && c(3) == 0
    A = [ABK(:,1:n) zeros(n,(s-1)*n)];
    B = [ABK(:,n+1:n+r) zeros(n,(s-1)*r)];
    K = ABK(:,n+r+1:n+r+s*l);
elseif c(1) >= 1 && c(2) == 0 && c(3) == 1
    A = [ABK(:,1:n) zeros(n,(s-1)*n)];
    B = ABK(:,n+1:n+s*r);
    K = [ABK(:,n+s*r+1:n+s*r+l) zeros(n,(s-1)*l)];
elseif c(1) == 0 && c(2) == 0 && c(3) == 0
    A = ABK(:,1:s*n);
    B = ABK(:,s*n+1:s*n+s*r);
    K = ABK(:,s*n+s*r+1:s*n+s*r+s*l);
elseif c(1) == 0 && c(2) == 1 && c(3) == 0
    A = ABK(:,1:s*n);
    B = [ABK(:,s*n+1:s*n+r) zeros(n,(s-1)*r)];
    K = ABK(:,s*n+r+1:s*n+r+s*l);
elseif c(1) == 0 && c(2) == 0 && c(3) == 1
    A = ABK(:,1:s*n);
    B = ABK(:,s*n+1:s*n+s*r);
    K = [ABK(:,s*n+s*r+1:s*n+s*r+l) zeros(n,(s-1)*l)];
end

if stable
       %estimate noise covariance matrices
       if c(1)==0 && c(2)==0
       	WV = [x(:,2:end); Y(:,1:end-1)] - [A B; C zeros(l,s*r)]*vertcat(MUX(:,1:end-1),MUU(:,1:end-1));
       elseif c(1)==1 && c(2)==0
        WV = [x(:,2:end); Y(:,1:end-1)] - [A(:,1:n) B; C(:,1:n) zeros(l,s*r)]*vertcat(x(:,1:end-1),MUU(:,1:end-1));
       elseif c(1)==0 && c(2)==1
        WV = [x(:,2:end); Y(:,1:end-1)] - [A B(:,1:r); C zeros(l,r)]*vertcat(MUX(:,1:end-1),U(:,1:end-1));
       elseif c(1)==1 && c(2)==1
        WV = [x(:,2:end); Y(:,1:end-1)] - [A(:,1:n) B(:,1:r); C(:,1:n) zeros(l,r)]*vertcat(x(:,1:end-1),U(:,1:end-1));
       end
       QSR = (WV*WV')./size(WV,2);
       Q = QSR(1:n,1:n);
       S = QSR(1:n,n+1:n+l);
       R = QSR(n+1:n+l,n+1:n+l);
       % find the range of the scheduling
       murange = [min(mu(2:end,:),[],2),max(mu(2:end,:),[],2)]; 
       % find the hinf optimal observer, guaranteed stable in the range of the scheduling
       try
         [K,gamma] = lobshinfsyn(A,C,Q,R,S,murange,c(3));
       catch exception
           warning('lx2abck:lobshinfsynerror',['In LX2ABCK, the function LOBSHINFSYN returned an error: \n ---------------- \n', exception.message, '\n ----------------\n LX2ABCK returns the K calculated by least squares, not guaranteed to stabilize the predictor form.'])
           gamma = [];
       end       
end

if pred
    for i = 1:s
        A(:,(i-1)*n+1:i*n) = A(:,(i-1)*n+1:i*n) - K(:,(i-1)*l+1:i*l)*C(:,1:n);
    end
end