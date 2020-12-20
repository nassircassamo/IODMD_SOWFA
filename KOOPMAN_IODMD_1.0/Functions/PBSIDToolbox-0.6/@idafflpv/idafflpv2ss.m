function [sys] = idafflpv2ss(sys,p)
%IDAFFLPV/IDAFFLPV2SS Convert affine LPV model to multiple LTI SS models.
%  [M] = IDAFFLPV2SS(M,MU) converts IDAFFLPV model M to multiple SS
%  models for every time step in the scheduling vector MU.

% determine lpv system sizes
[l,r,n,m] = size(sys);
N = size(p,1);

% get system matrices
[Alpv,Blpv,Clpv,Dlpv,Klpv] = getABCDK(sys);
Ts = sys.Ts;
Name = sys.Name;
InputName = sys.InputName;
OutputName = sys.OutputName;
for i = 1:l
    InputName = [InputName; 'Innovation'];
end
StateName = sys.StateName;
L = chol(sys.NoiseVariance,'lower');

% allocate cells
A = zeros(n,n,N);
B = zeros(n,r+l,N);
C = zeros(l,n,N);
D = zeros(l,r+l,N);

% construct ss models
for i = 1:N
    A(:,:,i) = Alpv(1:n,1:n);
    B(:,:,i) = [Blpv(1:n,1:r) Klpv(1:n,1:l)*L];
    C(:,:,i) = Clpv(1:l,1:n);
    D(:,:,i) = [Dlpv(1:l,1:r) eye(l,l)*L];
    for j = 1:m
        A(:,:,i) = A(:,:,i) + Alpv(1:n,(j-1)*n+1+n:j*n+n).*p(i,j);
        B(:,:,i) = B(:,:,i) + [Blpv(1:n,(j-1)*r+1+r:j*r+r) Klpv(1:n,(j-1)*l+1+l:j*l+l)*L].*p(i,j);
        C(:,:,i) = C(:,:,i) + Clpv(1:l,(j-1)*n+1+n:j*n+n).*p(i,j);
        D(:,:,i) = D(:,:,i) + [Dlpv(1:l,(j-1)*r+1+r:j*r+r) zeros(l)].*p(i,j);
    end
end
warning off;
sys = ss(A,B,C,D,Ts,'Name',Name,'InputName',InputName,'OutputName',OutputName,'StateName',StateName);
warning on;