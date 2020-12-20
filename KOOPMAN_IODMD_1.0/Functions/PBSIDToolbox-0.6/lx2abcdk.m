function [A,B,C,D,K,gamma] = lx2abcdk(x,u,y,mu,f,p,c,pred,stable)
%LX2ABCDK Estimates the matrices A, B, C and D of the LPV state space model
%  [A,B,C,D,K] = lx2abcdk(x,u,y,mu,f,p) estimates the matrices A, B, C, and
%  D of the state space model:
%
%     x(k+1) = A kron(mu(k),x(k)) + B kron(mu(k),u(k)) + K kron(mu(k),e(k))
%     y(k)   = C x(k) + D u(k) + e(k)
%
%  using the knowledge of the state vector x, the input vector u, the
%  output vector y and the scheduling vector mu. The past window size p is
%  recomended to be higher then the expected system order n. Future window
%  size f must equal or smaller then past window size p. The outputs are
%  the linear parameter-varying matrices A, B, K, where matrices have the
%  form of A=[A(1) A(2) ... A(m)].
%
%  [A,B,C,D,K]=lx2abcdk(x,u,y,mu,f,p,c) specifies which of the system
%  matrices are constant and not parameter-varing. For each of the matrices
%  A, B, and K an 1 or 0 can be given in the vector c. (default C=[0 0 0])
%
%  [A,B,C,D,K,gamma]=lx2abcdk(x,u,y,mu,f,p,c,pred,stable) with
%
%  pred = 1 (default: pred=0) returns the system matrices At, Bt, K, C, D 
%  of the predictor form:
%     x_est(k+1) = At kron(mu(k),x(k)) + Bt kron(mu(k),u(k)) +
%                 + K kron(mu(k),y(k))
%     y_est(k)   = C x(k) + D u(k)
%   At = [A(1)-K(1)*C, ..., A(m)-K(m)*C]
%   Bt = [B(1)-K(1)*D, ..., B(m)-K(m)*D]
%
%  stable = 1 (default: stable=0) designs K such that the predictor form is 
%  stable in the range of the scheduling mu, and that the predictor is Hinf
%  optimal w.r.t the noise in the model (see also the function
%  lobshinfsyn). In that case, gamma returns the Hinf norm of the transfer
%  function between noise and state error i.e. a smaller indicates better 
%  estimation of the states given the noise properties.
%
%  If stable = [mumin,mumax], a two-column matrix then K is designed such 
%  that the predictor form is stable in the range  mumin <= mu <= mumax.
%
%  References:
%    [1] J.W. van Wingerden, and M. Verhaegen, ``Subspace identification
%    of Bilinear and LPV systems for open- and closed-loop data'',
%    Automatica 45, pp 372--381, 2009.
%
%  See also: lvordarx, lmodx.m, lx2abck.m, lobshinfsyn.m.

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
    error('LX2ABCDK requires at least six input arguments.')
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
n = size(x,1);
s = size(mu,1);

% find the range of the scheduling for which the predictor should be stable
if (size(stable,2) == 2 && size(stable,1) == s)
    murange = stable(2:end,:);
    stable = 1;
elseif stable==1
    murange = [min(mu(2:end,:),[],2),max(mu(2:end,:),[],2)];
elseif stable~=0
    error('Invalid size of matrix STABLE (the range of the scheduling for which the predictor should be stable. First row of this matrix should be ones.')
end
    
if l == 0
    error('LX2ABCDK requires an output vector y.')
end
if n == 0
    error('LX2ABCDK requires an state vector x.')
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
    error('Future window size f must be equal or smaller than past window p. (f <= p)')
end

% remove the window sizes from input and output vector
U = u(:,p+1:p+size(x,2));
Y = y(:,p+1:p+size(x,2));
MU = mu(:,p+1:p+size(x,2));
MUX = khatrirao(MU,x);
MUU = khatrirao(MU,U);
    
CD = Y*pinv(vertcat(x,U));
if pred
    E = Y;
else
    E = Y - CD*vertcat(x,U);
end
if c(3) == 0
    MUE = khatrirao(MU,E);
end
    
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

% obtaining the LPV matrices C and D
C = [CD(:,1:n) zeros(l,(s-1)*n)];
D = [CD(:,n+1:n+r) zeros(l,(s-1)*r)];

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
elseif c(1) == 0 && c(2) == 1 && c(3) == 1
    A = ABK(:,1:s*n);
    B = [ABK(:,s*n+1:s*n+r) zeros(n,(s-1)*r)];
    K = [ABK(:,s*n+s*r+1:s*n+s*r+l) zeros(n,(s-1)*l)];
elseif c(1) == 0 && c(2) == 0 && c(3) == 1
    A = ABK(:,1:s*n);
    B = ABK(:,s*n+1:s*n+s*r);
    K = [ABK(:,s*n+s*r+1:s*n+s*r+l) zeros(n,(s-1)*l)];
end

if stable
       %estimate noise covariance matrices
       if c(1)==0 && c(2)==0
       	WV = [x(:,2:end); Y(:,1:end-1)] - [A B; C D]*vertcat(MUX(:,1:end-1),MUU(:,1:end-1));
       elseif c(1)==1 && c(2)==0
        WV = [x(:,2:end); Y(:,1:end-1)] - [A(:,1:n) B; C(:,1:n) D]*vertcat(x(:,1:end-1),MUU(:,1:end-1));
       elseif c(1)==0 && c(2)==1
        WV = [x(:,2:end); Y(:,1:end-1)] - [A B(:,1:r); C D(:,1:r)]*vertcat(MUX(:,1:end-1),U(:,1:end-1));
       elseif c(1)==1 && c(2)==1
        WV = [x(:,2:end); Y(:,1:end-1)] - [A(:,1:n) B(:,1:r); C(:,1:n) D(:,1:r)]*vertcat(x(:,1:end-1),U(:,1:end-1));
       end
       QSR = (WV*WV')./size(WV,2);
       Q = QSR(1:n,1:n);
       S = QSR(1:n,n+1:n+l);
       R = QSR(n+1:n+l,n+1:n+l);
       % find the hinf optimal observer, guaranteed stable in the range of the scheduling
%       try
         [K,gamma] = lobshinfsyn(A,C,Q,R,S,murange,c(3));
%        catch exception
%            warning('lx2abcdk:lobshinfsynerror',['In LX2ABCDK, the function LOBSHINFSYN returned an error: \n ---------------- \n', exception.message, '\n ----------------\n LX2ABCDK returns the K calculated by least squares, not guaranteed to stabilize the predictor form.'])
%            gamma = [];
%        end       
end

if pred
    for i = 1:s
        A(:,(i-1)*n+1:i*n) = A(:,(i-1)*n+1:i*n) - K(:,(i-1)*l+1:i*l)*C(:,1:n);
        B(:,(i-1)*r+1:i*r) = B(:,(i-1)*r+1:i*r) - K(:,(i-1)*l+1:i*l)*D(:,1:r);
    end
end
end