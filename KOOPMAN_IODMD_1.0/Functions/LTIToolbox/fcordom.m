function [S,Rnew] = fcordom(H,w,s) 
%FCORDOM    Delivers information about the order of the continuous 
%           time LTI state space model from frequency data. 
%           It acts as a pre-processor function for 
%           fcmodom. The latter actually estimates 
%           the  system matrices A and C. 
%
% Model structure: 
%                          -1 
%           H = C (w I - A)   B  +  D 
% 
% Syntax: 
%           [Sn,Rnew]=fcordom(H,w,s); 
% 
% Input: 
%  H,w     The measured FRF and the complex frequencies at which it is 
%          measured. 
%  s       The dimension parameter that determines the number 
%          of block rows in the processed Hankel matrices. 
%          This parameter should be chosen larger than the expected 
%          system order. The optimal value has to be found by trial 
%          and error. Generally twice as large is a good starting value. 
% 
% Output: 
% Sn       Singular values bearing information on the order 
%          of the system. 
% Rnew     Data structure used by fcmodom for the estimation of A 
%          and C. 
% 
% See also: fcmodom, fdordom 
 
% The structure of this function is based loosely on 
% the SMI-2.0 function 'dordom' by B. Haverkamp and 
% ff2ss by T. McKelvey. 
% Revised by Niek Bergboer, 2001 
% Revised by Ivo Houtzager, 2007
% Copyright (c) 2001-2007, Delft Center of Systems and Control 

if nargin ~= 3
    error('FCORDOM requires 3 input arguments.');
end
Rold = [];
if ndims(H) ~= 3
    error('H must be a 3D array.');
end
l = size(H,1);
m = size(H,2);
N = size(H,3);

if size(w,1) ~= N
    error('w must contain the same number of frequencies as H.');
end
if size(w,2) ~= 1
    error('w must have one column.');
end
if s < 2
    error('s must be at least 2.');
end
if l < 1
    error('FCMODOM requires an output.');
end
if m < 1
    error('FCMODOM requires an input.');
end
if 2*N*m < l*s
    error('The number of samples is too small.');
end
if max(abs(real(w))) > 10*eps
    error('At least one complex frequency is not imaginary.');
end

wscale = (max(imag(w)) + min(imag(w)))/2;
w = w/wscale;
[G,Zo] = forsythe(reshape(H,l,m*N),w,s);
Wm = forsythe(kron(ones(1,N),speye(m)),w,s);

if isempty(Rold),
    M = zeros(2*N*m,s*(l+m));
else
    M = zeros(2*N*m+s*(l+m),s*(l+m));
    M(2*N*m+1:2*N*m+s*(m+l),:) = tril(Rold(:,1:s*(m+l)))';
end;
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

[Un,Sn] = svd(R22);
clear Vn;
S = diag(Sn);
S = S(1:s,1);

Rnew = zeros(s*(m+l),s*(m+2*l));
Rnew(:,1:s*(m+l)) = R;
Rnew(s*m+1:s*(m+l),s*(m+l)+1:s*(m+2*l)) = Un;
Rnew(1:s,s*(m+l)+1:s*(m+l)+l) = Zo;
Rnew(1,2) = m;
Rnew(1,3) = l;
Rnew(1,4) = s;
Rnew(2,3) = wscale;
end

function [M,Z] = forsythe(BRow,MFac,r)
MFac = MFac(:);
l = size(BRow,1);
N = size(MFac,1);
m = size(BRow,2)/N;

% Allocate space
M = zeros(l * r,N * m);
Z = zeros(r,l);

% Calculate multiplication matrix for each block-row
Dw = spdiags(kron(MFac(:),ones(m,1)),0,N*m,N*m);

M(1:l,:) = BRow;
for j = 1:l
    Z(1,j) = norm(M(j,:),2);
    M(j,:) = M(j,:)/Z(1,j);
end

M(l+1:2*l,:) = M(1:l,:)*Dw;
for j = 1:l
    Z(2,j) = norm(M(l+j,:),2);
    M(l+j,:) = M(l+j,:)/Z(2,j);
end

for k = 2:r-1
    M(k * l+1:(k+1)*l,:) = M((k-1)*l+1:k*l,:)*Dw;
    for j = 1:l
        M(k * l+j,:) = M(k * l+j,:) + M((k-2)*l+j,:)*Z(k,j);
        Z(k+1,j) = norm(M(k * l+j,:),2);
        M(k * l+j,:) = M(k * l+j,:)/Z(k+1,j);
    end
end
end






