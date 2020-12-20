function [sys,u,y] = idafflpvA2idss(sys,u,y,mu)
%IDAFFLPV/IDAFFLPVA2IDSS Convert IDAFFLPV (constant A) model to single IDSS model.
%  [M] = IDAFFLPVA2IDSS(M) converts IDAFFLPV model M (with constant A) to
%  single IDSS models. Usefull for discretization of LPV models with
%  constant A.

% determine lpv system sizes
[l,r,n,m] = size(sys);

% get system matrices
[Alpv,Blpv,Clpv,Dlpv,Klpv] = getABCDK(sys);
Ts = sys.Ts;
Name = sys.Name;
X0 = sys.x0;
lambda = sys.NoiseVariance;
InputName = sys.InputName;
OutputName = sys.OutputName;
StateName = {sys.StateName{1:n}};

% allocate cells
A = Alpv(1:n,1:n);
B = Blpv;
C = zeros(l*(m+1),n);
D = zeros(l*(m+1),r*(m+1));
K = Klpv;

% construct ss models
for j = 1:(m+1)
    C((j-1)*l+1:j*l,:) = Clpv(:,(j-1)*n+1:j*n);
    D((j-1)*l+1:j*l,(j-1)*r+1:j*r) = Dlpv(:,(j-1)*r+1:j*r);
end
sys = idss(A,B,C,D,K,X0,Ts,'Name',Name,'InputName',InputName,'OutputName',OutputName,'StateName',StateName);
sys.x0 = X0;
lambda = kron(eye(m+1),lambda);
sys = pvset(sys,'NoiseVariance',lambda);


if nargin > 1
    if size(u,2) < size(u,1)
        u = u';
    end
    if size(y,2) < size(y,1)
        y = y';
    end
    if size(mu,2) < size(mu,1)
        mu = mu';
    end
    mu = [ones(1,size(mu,2)); mu];
    u = khatrirao(mu,u)';
    y = khatrirao(mu,y)';
end