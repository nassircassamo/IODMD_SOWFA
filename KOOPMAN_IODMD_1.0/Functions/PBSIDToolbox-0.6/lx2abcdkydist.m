function [A,B,C,D,K,gamma] = lx2abcdkydist(x,u,d,y,mu,f,p,c,pred,stable)
%LX2ABCDKYDIST is LX2ABCDK for a special case with output disturbances
%   x(k+1) = A kron(mu(k),x(k)) + B kron(mu(k),u(k)) + K kron(mu(k),e(k))
%   y(k)   = C x(k) + Du u(k) + Dd kron(mu(k),d(k)) + e(k)
%   
%  if c(4)=1 then Dd is not varying with the scheduling mu
%
%  See also: lvordarx, lmodx.m, lx2abck.m, lobshinfsyn.m.

%  Pieter Gebraad
%  Delft Center of Systems and Control
%  Delft University of Technology
%  The Netherlands, 2011

gamma = [];

% check number of input arguments
if nargin < 7
    error('LX2ABCDK requires at least seven input arguments.')
end
if nargin < 10 || isempty(stable)
    stable = 0;
end
if nargin < 9 || isempty(pred)
    pred = 0;
end
if nargin < 8 || isempty(c)
    c = [0 0 0 0];
end

% check dimensions of inputs
if size(y,2) < size(y,1)
    y = y';
end
if size(mu,2) < size(mu,1)
    mu = mu';
end
if size(d,2) < size(d,1)
    d = d';
end
N = size(y,2);
l = size(y,1);
n = size(x,1);
s = size(mu,1);

rd = size(d,1);
if ~isequal(N,length(d))
    error('The number of rows of vectors/matrices d and y must be the same.')
end

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
d = d(:,p+1:p+size(x,2));
Y = y(:,p+1:p+size(x,2));
MU = mu(:,p+1:p+size(x,2));
if c(4)==0
    d = khatrirao(MU,d);
    rd = rd*s;
end
MUX = khatrirao(MU,x);
MUU = khatrirao(MU,U);
    
CD = Y*pinv(vertcat(x,U,d));
if pred
    E = Y;
else
    E = Y - CD*vertcat(x,U,d);
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
D = [CD(:,n+1:n+(r+rd)) zeros(l,(s-1)*(r+rd))];

% obtaining the LPV matrices A, B and K
if c(1) >= 1 && c(2) == 1 && c(3) == 1
    A = [ABK(:,1:n) zeros(n,(s-1)*n)];
    B = [ABK(:,n+1:n+r) zeros(n,(s-1)*(r))];
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

Bu = B;
B = zeros(n,(r+rd)*s);
for i = 1:s 
    B(:,(i-1)*(r+rd)+1:(i-1)*(r+rd)+r) = Bu(:,(i-1)*r+1:i*r);
end

MUU = khatrirao(MU,vertcat(U,d));

if stable
       %estimate noise covariance matrices
       if c(1)==0 && c(2)==0
       	WV = [x(:,2:end); Y(:,1:end-1)] - [A B; C D]*vertcat(MUX(:,1:end-1),MUU(:,1:end-1));
       elseif c(1)==1 && c(2)==0
        WV = [x(:,2:end); Y(:,1:end-1)] - [A(:,1:n) B; C(:,1:n) D]*vertcat(x(:,1:end-1),MUU(:,1:end-1));
       elseif c(1)==0 && c(2)==1
        WV = [x(:,2:end); Y(:,1:end-1)] - [A B(:,1:r); C D(:,1:r+rd)]*vertcat(MUX(:,1:end-1),U(:,1:end-1),d(:,1:end-1));
       elseif c(1)==1 && c(2)==1
        WV = [x(:,2:end); Y(:,1:end-1)] - [A(:,1:n) B(:,1:r); C(:,1:n) D(:,1:r+rd)]*vertcat(x(:,1:end-1),U(:,1:end-1),d(:,1:end-1));
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
        B(:,(i-1)*r+1:i*r) = B(:,(i-1)*r+1:i*r) - K(:,(i-1)*l+1:i*l)*D(:,1:r+rd);
    end
end
end