function [sys] = idss2idafflpvA(sys,m)
%IDAFFLPV/IDSS2IDAFFLPVA Convert single IDSS model to IDAFFLPV (constant A) model
%  [M] = IDSS2IDAFFLPVA(M,p) converts single IDSS model to IDAFFLPV model M
%  (with constant A). Usefull for discretization of LPV models with
%  constant A.

% determine lpv system sizes
n = size(sys.a,1);
[l,r] = size(sys);
l = l/(m+1);
r = r/(m+1);

% get system matrices
[A,B,C,D,K,X0] = ssdata(sys);
Ts = sys.Ts;
Name = sys.Name;
X0 = sys.x0;
lambda = sys.NoiseVariance;
InputName = sys.InputName{1:r};
OutputName = sys.OutputName{1:l};
StateName = sys.StateName;
for i = 1:m
    StateName = [StateName; sys.StateName];
end

% allocate cells
Alpv = [A zeros(n,m*n)];
Blpv = B;
Klpv = K;
Clpv = zeros(l,n*(m+1));
Dlpv = zeros(l,r*(m+1));

% construct ss models
for j = 1:(m+1)
    Clpv(:,(j-1)*n+1:j*n) = C((j-1)*l+1:j*l,:);
    Dlpv(:,(j-1)*r+1:j*r) = D((j-1)*l+1:j*l,(j-1)*r+1:j*r);
end
sys = idafflpv(Alpv,Blpv,Clpv,Dlpv,Klpv,X0,Ts,'Name',Name,'InputName',InputName,'OutputName',OutputName,'StateName',StateName);
sys.x0 = X0;
sys.NoiseVariance = lambda(1:l,1:l);