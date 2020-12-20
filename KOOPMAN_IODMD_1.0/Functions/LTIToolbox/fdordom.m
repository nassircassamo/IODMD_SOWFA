function [S,Rnew] = fdordom(H,w,s,Rold)
%FDORDOM    Delivers information about the order of the discrete 
%           time LTI state space model from frequency data. 
%           It acts as a pre-processor function for 
%           fdmodom. The latter actually estimates 
%           the  system matrices A and C. 
%
% Model structure: 
%                          -1 
%           H = C (w I - A)   B  +  D 
% 
% Syntax: 
%           [Sn,Rnew]=fdordom(u,y,s); 
%           [Sn,Rnew]=fdordom(u,y,s,Rold); 
% 
% Input: 
%  H,w     The measured FRF and the complex frequencies at which it is 
%          measured. The FRF should be a matrix of size l x m x N, in 
%          which H(:,:,i) is the complex FRF at the i-th complex frequency. 
%          w should be a vector of size N x 1 containing the complex 
%          frequencies. 
%  s       The dimension parameter that determines the number 
%          of block rows in the processed Hankel matrices. 
%          This parameter should be chosen larger than the expected 
%          system order. The optimal value has to be found by trial 
%          and error. Generally twice as large is a good starting value. 
%  Rold    Should not be there on the first call of 
%          the routine. When a second (or third, ...) 
%          I/O data batch is processed a present R 
%          matrix is used to store the information 
%          from these different data batches. 
% 
% Output: 
% Sn       Singular values bearing information on the order 
%          of the system. 
% Rnew     Data structure used by fdmodom for the estimation of A 
%          and C or by a next call to fdordom. 
% 
% See also: fdmodom, fcordom 
 
% The structure of this function is based loosely on 
% the SMI-2.0 function 'dordom' by B. Haverkamp and 
% ff2ss by T. McKelvey. 
% Revised by Niek Bergboer, 2001 
% Revised by Ivo Houtzager, 2007
% Copyright (c) 2001-2007, Delft Center of Systems and Control 

if (nargin ~= 3) && (nargin ~= 4)
    error('FDORDOM requires 3 or 4 input arguments.');
end
if nargin == 3
    Rold = [];
end

if ndims(H) ~= 3
    error('H must be a 3D array.');
end
l = size(H,1);
m = size(H,2);
N = size(H,3);

if size(w,1)~=N
    error('w must contain the same number of frequencies as H.');
end
if size(w,2)~=1
    error('w must have one column.');
end
if s < 2
    error('s must be at least 2.');
end
if l < 1
    error('FDMODOM requires an output.');
end
if m < 1
    error('FDMODOM requires an input.');
end
if 2*N*m < l*s
    error('The number of samples is too small.');
end

if max(abs(abs(w)-1))>10*eps
    error('At least one complex frequency is not unit-magnitude.');
end

if ~isempty(Rold)
    if (size(Rold,1)~=s * (m+l)) || (size(Rold,2)~=s * (m+2 * l))
        error('The size of the old R-factor is incorrect.');
    end
    if Rold(1,2) ~= m
        error('The number of inputs does not correspond to the old R-factor.');
    end
    if Rold(1,3) ~= l
        error('The number of outputs does not correspond to the old R-factor.');
    end
    if Rold(1,4) ~= s
        error('The blocksize does not correspond to the old R-factor.');
    end
    if Rold(2,3) ~= 0
        error('The previous R-factor seems to belong to a continuous-time model.');
    end
end

G = zeros(s*l,m*N);
Wm = sparse([],[],[],s*m,m*N,m*s*N);
G(1:l,:) = reshape(H,l,m*N);
Wm(1:m,:) = kron(ones(1,N),speye(m));

for j = 1:ceil(log2(s))
    if 2^j<=s
        % Full case
        rowsrc = 1:m*2^(j-1);
        rowdst = m*2^(j-1)+1:m*2^j;
        rowsrcG = 1:l*2^(j-1);
        rowdstG = l*2^(j-1)+1:l*2^j;
    else
        % Last few block-rows
        rowsrc = 1:m*(s-2^(j-1));
        rowdst = m*(2^(j-1))+1:m*s;
        rowsrcG = 1:l*(s-2^(j-1));
        rowdstG = l*(2^(j-1))+1:l*s;
    end
    % Complex frequency multiplication
    Xim = kron(spdiags(w.^(2^(j-1)),0,N,N),eye(m));
    Wm(rowdst,:) = Wm(rowsrc,:)*Xim;
    G(rowdstG,:) = G(rowsrcG,:)*Xim;
end
clear Xim;

if isempty(Rold)
    M = zeros(2*N*m,s*(l+m));
else
    M = zeros(2*N*m+s*(l+m),s*(l+m));
    M(2*N*m+1:2*N*m+s*(m+l),:) = tril(Rold(:,1:s * (m+l)))';
end
M(1:N*m,1:s*m) = real(Wm)';
Wm = imag(Wm);
M(N*m+1:2*N*m,1:s*m) = Wm';
clear Wm;
M(1:N*m,s*m+1:s*(l+m)) = real(G)';
G = imag(G);
M(N*m+1:2*N*m,s*m+1:s*(m+l)) = G';
clear G;

R = qr(M);
clear M;
R = triu(R(1:s*(m+l),:))';
R22 = R(s*m+1:s*(m+l),s*m+1:s*(m+l));

Wl = sparse([],[],[],s*l,l*N,l*s*N);
Wl(1:l,:) = kron(ones(1,N),speye(l));
for j = 1:ceil(log2(s))
    if 2^j<=s
        % Full case
        rowsrc = 1:l * 2^(j-1);
        rowdst = l * 2^(j-1)+1:l * 2^j;
    else
        % Last few block-rows
        rowsrc = 1:l * (s-2^(j-1));
        rowdst = l * (2^(j-1))+1:l * s;
    end
    % Complex frequency multiplication
    Xil = kron(spdiags(w.^(2^(j-1)),0,N,N),eye(l));
    Wl(rowdst,:) = Wl(rowsrc,:) * Xil;
end
clear Xil;

KTK = real(Wl*Wl');
clear Wl
if ~isempty(Rold)
    KTK = KTK+kron(tril(Rold(1:m:m*s,1:m:m*s))*tril(Rold(1:m:m*s,1:m:m*s))',speye(l));
end
[K,msg] = chol(KTK);
clear KTK;
if msg~=0
    error('FORDOM: Frequency Vandermonde matrix is non positive definite.');
end

[Un,Sn] = svd(full(K\R22));
clear Vn;
Un = full(K*Un);
clear K;
S = diag(Sn);
S = S(1:s,1);

Rnew = zeros(s*(m+l),s*(m+2*l));
Rnew(:,1:s*(m+l)) = R;
Rnew(s*m+1:s*(m+l),s*(m+l)+1:s*(m+2*l)) = Un;
Rnew(1,2) = m;
Rnew(1,3) = l;
Rnew(1,4) = s;
