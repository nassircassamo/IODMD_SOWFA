function [P,sigma,dA,dB,dC,dK] = dvar4abck(x,u,y,f,p,A,B,C,K,U,Zps)
%DVAR4ABCK Asymptotic variance of the PBSIDopt (VARX only) estimation
%  P=dvar4abck(x,u,f,p,A,B,C,K,U,Zps) returns the covariance of the
%  estimated state space matrices and acts as a pre-processor for dvar2frd.
%  The latter is used to calculate the probalistic error bounds around the
%  identified bode diagrams. The data matrices U and Zps can be obtained
%  from dordvarx. 
%
%  [P,sigma]=dvar4abck(x,u,f,p,A,B,C,K,U,Zps) also returns the covariance
%  matrix of the innovation noise. 
%
%  [P,sigma,dA,dB,dC,dK]=dvar4abck(x,u,f,p,A,B,C,K,U,Zps)

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check number if input arguments
if nargin < 11
    error('DVAR4ABCK requires eleven input arguments.')
end

% check the size of the windows
if f > p
    error('Future window size f must equal or smaller then past window p. (f <= p)')
end

% check dimensions of inputs
if size(y,2) < size(y,1)
    y = y';
end
if size(x,2) < size(x,1)
    x = x';
end
N = size(y,2);
l = size(y,1);
n = size(x,1);
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
if l == 0
    error('DVAR4ABCK requires an output vector y.')
end

% store the past and future vectors
m = r+l;
z = [u; y];
Z = zeros(p*m,N-p);
for i = 1:p
    Z((i-1)*m+1:i*m,:) = z(:,i:N+i-p-1);
end

% select the only the system order
U = U(1:n,:);
if size(Zps,2)/m > p
    Zps = Zps(:,1:p*m);
end

% remove the window sizes from input and output vector
u = u(:,p+1:p+size(x,2));
y = y(:,p+1:p+size(x,2));

% calculate the innovation sequence
e = y - C*x;
sigma = (e*e')/length(e);

% asymptotic variance
LL  = pinv([x(:,1:end-1); u(:,1:end-1); e(:,1:end-1)]);
LL2 = pinv(x);
Term1 = ObsContSum(Zps,speye(N-p),l,r,f,p,LL,Z(:,2:end),U);
Term2 = ObsContSum(Zps,speye(N-p),l,r,f,p,LL,Z(:,1:end-1),A*U);
Term3 = ObsContSum(Zps,speye(N-p),l,r,f,p,LL2,Z,-C*U);
alpha1P = Term1-Term2;
alpha2P = Term3;
beta1 = -kron(LL'*Z(:,1:end-1)'*Zps',K);
beta2 = -kron(LL2'*Z'*Zps',eye(l));
if l==1
    P = sigma*[alpha1P+beta1;alpha2P+beta2]*[alpha1P+beta1;alpha2P+beta2]';
else
    P = [alpha1P+beta1;alpha2P+beta2]*sparse(kron(speye(N-p),sigma))*[alpha1P+beta1;alpha2P+beta2]';
end

if nargout > 2
    dTh = [alpha1P+beta1;alpha2P+beta2]*e(:);
    dA = reshape(dTh(1:n*n,1),n,n);
    dB = reshape(dTh(n*n+1:n*n+n*r,1),n,r);
    dK = reshape(dTh(n*n+n*r+1:n*n+n*r+n*l,1),n,l);
    dC = reshape(dTh(n*n+n*r+n*l+1:n*n+n*r+n*l+n*l,1),l,n);
end

end

function SumKron = ObsContSum(Zps,Y,l,r,f,p,LL,Z,S)
q = size(Y,1);

for i = 1:p
    CK(:,1+(l+r)*(i-1):(l+r)*i) = Y*Zps(:,1+(l+r)*(i-1):(l+r)*i);
end

SumKron = zeros(size(LL,2)*size(S,1),size(Y,2)*l);
for i = 1:f
    GammaK = zeros(q,(l+r)*p);
    GammaK(:,1+(l+r)*(i-1):(l+r)*p) = CK(:,1:(l+r)*(p+1-i));
    SumKron = SumKron + kron(LL'*Z'*GammaK',S(:,1+l*(i-1):l*i));
end
end



