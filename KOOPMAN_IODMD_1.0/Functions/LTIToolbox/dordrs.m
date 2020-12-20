function [Sn,Rnew]=dordrs(u,y,x,s,Rold)
%DORDRS     An iterative method to extract information about
%           the LTI state space model, based on the
%           reconstructed state. With the updated state
%           space model, a new reconstructed state can be
%           computed, using dmodrs. 
%           Then dordrs can again try to find an R matrix.
%           To initialize the iteration, the reconstructed
%           state x can be obtained from a model,
%           estimated for instance with dordpi.
%           Concatentation of data sets is possible by use
%           of the extra parameter R (Rold,Rnew). 
%           Model structure: 
%                 x(k+1) = Ax(k) + Bu(k) 
%                 y(k)   = Cx(k) + Du(k) + v(k)
%           where v(k) is zero-mean noise of arbitary color, 
%           independent of the noise-free input u(k).
%
% Syntax:
%           [Sn,Rnew]=dordrs(u,y,x,s);
%           [Sn,Rnew]=dordrs(u,y,x,s,Rold);
%
% Input:
%   u,y     The input respectively output data of the system to be 
%           identified.
%   x       Reconstructed state.
%   s       The dimension parameter that determines the number
%           of block rows in the processed Hankel matrices
%   Rold    Should not be there on the first call of
%           the routine. When a second (or third, ...)
%           i/o data batch is processed a present R
%           matrix is used to store the information
%           from these different data batches.
% 
% Output:
%   Sn      Singular values bearing information on the order
%           of the system
%   Rnew    The compressed lower triangular factor with
%           additional information (such as i/o dimension etc.)
%           stored in the zeros.
%
% See also: dmodpi, dordpi, dmodpo, dordpo, dmodrs

% Michel Verhaegen 11-01-1990
% Revised by Ivo Houtzager, 2007
% Copyright (c) 1990-2007, Delft Center of Systems and Control 

% Check number of arguments
if nargin < 4
    error('DORDRS requires at least four input arguments.');
end
if size(y,2) > size(y,1)
    y = y';
end
if size(u,2) > size(u,1)
    u = u';
end
if size(x,2) > size(x,1)
    x = x';
end

N = size(y,1);
l = size(y,2);
m = size(u,2);
n = size(x,2);

if m == 0
    error('DORDRS requires an input.')
end
if n == 0
    error('DORDRS requires an estimated state vector.')
end
if l == 0
    error('DORDRS requires an output.')
end
if ~(size(u,1) == N)
    error('Input and output should have same lenght.')
end
if ~(size(x,1) == N)
    error('State-vector should have same length as input.')
end
if s < 2
    error('s should be at least 2.')
end
if s*(m+l)+n > N-2*s-1
    error('s is chosen too large or number of datapoints is too small.')
end
if nargin < 5,
    Rold = [];
elseif ~((size(Rold,1)==s*(2*m+n+l)) && (size(Rold,2)==2*(m+l)*s+n*s))
    error('R-matrix has unexpected size.')
else
    Rold = tril(Rold(:,1:(2*m+l)*s));
end

% construction of Hankel matrices
NN = N-2*s+1;

Yf = zeros(NN,l*s);
Uf = zeros(NN,m*s);
Up = zeros(NN,m*s);
Xe = zeros(NN,n*s);
for i = (1:s),
    Up(:,(i-1)*m+1:i*m) = u(i:NN+i-1,:);
    Uf(:,(i-1)*m+1:i*m) = u(s+i:NN+s+i-1,:);
    Yf(:,(i-1)*l+1:i*l) = y(s+i:NN+s+i-1,:);
    Xe(:,(i-1)*n+1:i*n) = x(s+i:NN+s+i-1,:);
end

Rnew = triu(qr([Rold'; [Uf Up Xe Yf]]));
Rnew = Rnew(1:(2*m+n+l)*s,1:(2*m+n+l)*s)';
R32 = Rnew(2*m*s+s*n+1:(2*m+n+l)*s,m*s+1:(2*m+n)*s);
[Un,Sn] = svd(R32);
Sn = diag(Sn);
Sn = Sn(1:s);
Rnew(1,2) = m;
Rnew(1,3) = l;
Rnew(1,4) = n;
Rnew(2,3) = s;
Rnew(2*m*s+n*s+1:(2*m+n+l)*s,(2*m+n+l)*s+1:2*(m+l)*s+n*s) = Un(:,1:l*s);
