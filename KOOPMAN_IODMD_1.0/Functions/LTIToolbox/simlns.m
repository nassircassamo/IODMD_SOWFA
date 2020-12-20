function U2 = simlns(A,B,C,K,fD,fx) 
%SIMLNS     This function returns the left null-space of a Lee-Poolla 
%           similarity map for a LTI state-space model. 
% 
% Syntax: 
%           U2 = simlns(A,B,C,[],[],[]); 
%           U2 = simlns(A,B,C,K,fD,fx); 
% 
% Input: 
%  A,B,C    System matrices describing the state space system. 
%  K        Kalman gain for the state-space model. This matrix 
%           is optional and should be passed as an empty matrix 
%           if not needed. 
%  fD       Indicates whether rows corresponding to the D-matrix 
%           should be present in the left null-space matrix. 
%  fx       Indicates whether rows corresponding to the initial state 
%           should be present in the left null-space matrix. 
% 
% Output: 
%  U2       The left null-space of the similarity map. 
% 
% Remarks: 
%           Specifying fx only causes an n x n identify-matrix to 
%           be appended to the lower right of the left null space; 
%           in a non-linear optimization, applying the left null-space 
%           ensures that the state-basis does not change. It thus does 
%           not have to be projected. 
% 
%           simlns is implemented both as an M-file (a reference 
%           implementation) and as a MEX-file. The MEX-file uses 
%           the QR-factorization in LAPACK. 
% 
% Used by: 
%           dfunlti, ffunlti 
% 
 
% Niek Bergboer, 2001 
% Revised by Ivo Houtzager, 2007
% Copyright (c) 2001-2007, Delft Center of Systems and Control 
 
if nargin < 6
    error('SIMLNS requires six input arguments');
end

% Get dimensions
n = size(A,1);
m = size(B,2);
l = size(C,1);

% Check dimensions of A, B, C and K
if size(A,2)~=n
    error('A matrix must be square.');
end
if size(B,1)~=n
    error('B matrix has wrong number of rows.');
end
if size(C,2)~=n
    error('C matrix has wrong number of columns.');
end
if ~isempty(K)
    if size(K,1)~=n
        error('K matrix has wrong number of rows.');
    end
    if size(K,2)~=l
        error('K matrix has wrong number of columns.');
    end
    fK = 1;
else
    fK = 0;
end
if ~isempty(fD)
    if fD~=1 && fD~=0
        error('fD should be empty, 0 or 1.');
    end
else
    fD = 0;
end

if ~isempty(fx)
    if fx~=1 && fx ~=0,
        error('fx should be empty, 0 or 1.');
    end
else
    fx = 0;
end

% Set number of rows and columns in M and allocate memory.
SizeR = n*n + n*m + l*n + fD*l*m + fK*n*l;
SizeC = n*n;
M = zeros(SizeR,SizeC);
U2 = zeros(SizeR + fx*n,SizeR - SizeC + fx*n);

M(1:n*(n+l),:) = kron(eye(n),[A;C]);
M(1:n*(n+l),:) = M(1:n*(n+l),:)-kron(A',[eye(n);zeros(l,n)]);
if fD
    M(n*(n+l)+1:(n+m)*(n+l),:) = -kron(B',[eye(n);zeros(l,n)]);
else
    M(n*(n+l)+1:n*(n+l)+n*m,:) = -kron(B',eye(n));
end
if fK
    StartRowK = n*n + n*m + l*n + fD*l*m;
    M(StartRowK+1:StartRowK+n*l,:) = -kron(K',eye(n));
end

[Q,R] = qr(M);
clear M;
clear R;
U2(1:SizeR,1:SizeR-SizeC) = Q(:,SizeC+1:SizeR);
clear Q;

if fx
    U2(SizeR+1:SizeR+n,SizeR-SizeC+1:SizeR-SizeC+n) = eye(n);
end
 
 
 
 
 

