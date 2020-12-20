function [sys] = idafflpv2idss(sys,p)
%IDAFFLPV/IDAFFLPV2IDSS Convert IDAFFLPV model to multiple IDSS models.
%  [M] = IDAFFLPV2IDSS(M,MU) converts IDAFFLPV model M to multiple IDSS
%  models for every time step in the scheduling vector MU.

% determine lpv system sizes
[l,r,n,m] = size(sys);
N = size(p,1);

% get system matrices
[Alpv,Blpv,Clpv,Dlpv,Klpv] = getABCDK(sys);
Ts = sys.Ts;
Name = sys.Name;
X0 = sys.x0;
lambda = sys.NoiseVariance;
InputName = sys.InputName;
OutputName = sys.OutputName;
StateName = sys.StateName;

% allocate cells
A = zeros(n,n,N);
B = zeros(n,r,N);
C = zeros(l,n,N);
D = zeros(l,r,N);
K = zeros(n,l,N);

% construct ss models
for i = 1:N
    A(:,:,i) = Alpv(1:n,1:n);
    B(:,:,i) = Blpv(1:n,1:r);
    C(:,:,i) = Clpv(1:l,1:n);
    D(:,:,i) = Dlpv(1:l,1:r);
    K(:,:,i) = Klpv(1:n,1:l);
    for j = 1:m
        A(:,:,i) = A(:,:,i) + Alpv(1:n,(j-1)*n+1+n:j*n+n).*p(i,j);
        B(:,:,i) = B(:,:,i) + Blpv(1:n,(j-1)*r+1+r:j*r+r).*p(i,j);
        C(:,:,i) = C(:,:,i) + Clpv(1:l,(j-1)*n+1+n:j*n+n).*p(i,j);
        D(:,:,i) = D(:,:,i) + Dlpv(1:l,(j-1)*r+1+r:j*r+r).*p(i,j);
        K(:,:,i) = K(:,:,i) + Klpv(1:n,(j-1)*l+1+l:j*l+l).*p(i,j);
    end
end
sys = idss(A,B,C,D,K,X0,Ts,'Name',Name,'InputName',InputName,'OutputName',OutputName,'StateName',StateName);
sys.x0 = X0;
sys = pvset(sys,'NoiseVariance',lambda);