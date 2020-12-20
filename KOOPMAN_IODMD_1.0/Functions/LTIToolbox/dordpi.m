function [Sn,Rnew] = dordpi(u,y,s,Rold)
%DORDPI     Delivers information about the order of the LTI state
%           space model and acts as a pre-processor for
%           dmodpi. The latter actually estimates the
%           system matrices A and C. Concatentation
%           of data sets is possible by use of the extra 
%           parameter R (Rold,Rnew).
%           Model structure: 
%                x(k+1) = Ax(k) + Bu(k)
%                y(k)   = Cx(k) + Du(k) + v(k)
%           where v(k) is zero-mean noise of arbitary color, 
%           independent of the noise-free input u(k).
%
% Syntax:
%           [Sn,Rnew]=dordpi(u,y,s);
%           [Sn,Rnew]=dordpi(u,y,s,Rold);
%
% Input:
%   u,y     The input respectively output data of the system to be 
%           identified.
%   s       the dimension parameter that determines the number
%           of block rows in the processed Hankel matrices.
%           This parameter should be chosen larger than the expected
%           system order. The optimal value has to be found by trial
%           and error. Generally twice as large is a good starting value.
%   Rold    should not be there on the first call of
%           the routine. When a second (or third, ...)
%           i/o data batch is processed a present R
%           matrix is used to save the information
%           from these different data batches.
% 
% Output:
%   Sn      Singular values bearing information on the order
%           of the system.
%   Rnew    The compressed lower triangular factor with
%           additional information (such as i/o dimension 
%           stored in the zeros.
%
% See also: dmodpi, dordpo, dmodpo

% The dordpi (dmodpi) routine corresponds to the PI scheme 
% derived and analyzed in VERHAEGEN: "Subspace Model
% Identification. Part 3" Int. J. Control, Vol. 57.
%
% Michel Verhaegen 11-01-1990
% Revised by Michel Verhaegen, 1994
% Revised by Ivo Houtzager, 2007
% Copyright (c) 1990-2007, Delft Center of Systems and Control 

% Check number of arguments
if nargin < 3
    error('DORDPI requires at least three input arguments.');
end
if size(y,2) > size(y,1)
    y = y';
end
if size(u,2) > size(u,1)
    u = u';
end

N = size(y,1);
l = size(y,2);
m = size(u,2);

if l == 0
    error('DORDPI requires an output')
end
if m == 0
    error('DORDPI requires an input')
end
if ~(size(u,1) == N)
    error('Input and output should have same lenght')
end
if s < 2
    error('s should be at least 2')
end
if (2*m+l)*s >= N-2*s+1
    error('s is chosen too large or number of datapoints is too small')
end
if nargin < 4,
    Rold = [];
elseif ~((size(Rold,1)==s*(2*m+l)) && (size(Rold,2)==2*s*(m+l)))
    error('R-matrix has unexpected size.');
else
    Rold = tril(Rold(:,1:(2*m+l)*s));
end

% construction of Hankel matrices
NN = N-2*s+1;
Uf = zeros(NN,m*s);
Up = zeros(NN,m*s);
Yf = zeros(NN,l*s);
for i = (1:s),
    Up(:,(i-1)*m+1:i*m) = u(i:NN+i-1,:);
    Uf(:,(i-1)*m+1:i*m) = u(s+i:NN+s+i-1,:);
    Yf(:,(i-1)*l+1:i*l) = y(s+i:NN+s+i-1,:);
end 
Rnew = triu(qr([Rold'; [Uf Up Yf]]));
Rnew = Rnew(1:(2*m+l)*s,1:(2*m+l)*s)';
R32 = Rnew(2*m*s+1:(2*m+l)*s,m*s+1:2*m*s);
[Un,Sn] = svd(R32);
Sn = diag(Sn);
Sn = Sn(1:s);

Rnew(1,2) = m;
Rnew(1,3) = l;
Rnew(2,3) = s;
Rnew(2*m*s+1:(2*m+l)*s,(2*m+l)*s+1:2*(m+l)*s) = Un(:,1:l*s);


