function [Sn,Rnew]=dordpo(u,y,s,Rold)
%DORDPO     Delivers information about the order of the LTI
%           state space model and acts as a pre-processor for
%           dmodpo. The latter actually estimates
%           the  system matrices A and C  and Kalman gain.
%           Model structure: 
%                  x(k+1) = Ax(k) + Bu(k) + w(k)
%                  y(k)   = Cx(k) + Du(k) + v(k)
%           where w(k), v(k) are zero-mean white noise sequences,
%           independent of the noise-free input u(k).
%
% Syntax:
%           [Sn,Rnew]=dordpo(u,y,s);
%           [Sn,Rnew]=dordpo(u,y,s,Rold);
%
% Input:
%   u, y    The input respectively output data of the system to be 
%           identified. DORDPO does handle empty inputs.
%   s       The dimension parameter that determines the number
%           of block rows in the processed Hankel matrices.
%           This parameter should be chosen larger than the expected
%           system order. The optimal value has to be found by trial
%           and error. Generally twice as large is a good starting value.
%   Rold    Should not be there on the first call of
%           the routine. When a second (or third, ...)
%           i/o data batch is processed a present R
%           matrix is used to store the information
%           from these different data batches. 
%
% Output:
%   Sn      Singular values bearing information on the order
%           of the system.
%   Rnew    The compressed lower triangular factor with
%           additional information (such as i/o dimension etc.)
%           stored in the zeros.
%
% See also: dmodpi, dordpi, dmodpo

% Michel Verhaegen December 1994
% Revised by Michel Verhaegen, 1996
% Revised by Michel Verhaegen, 1997
% Revised by Ivo Houtzager, 2007
% Copyright (c) 1994-2007, Delft Center of Systems and Control 

% Check number of arguments
if nargin < 3
    error('DORDPO requires at least three input arguments.');
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

% po-moesp does handle empty inputs
if l == 0
    error('DORDPO requires an output')
end
if s < 2
    error('s must be at least 2')
end
if (~(size(u,1) == N) && ~isempty(u))
    error('Input and output should have same lenght')
end
if 2*(m+l)*s >= N-2*s+1
    error('s is chosen too large or number of datapoints is too small')
end
if nargin < 4
    Rold = [];
end
if ~isempty(Rold)
    if ~((size(Rold,1)==s*(2*m+2*l)) && (size(Rold,2)==max(4,s*(2*m+3*l))))
        error('R-matrix has unexpected size.')
    else
        Rold = tril(Rold(:,1:2*(m+l)*s));
    end
end

% construction of Hankel matrices
NN = N-2*s+1;
Up = zeros(NN,m*s);
Uf = zeros(NN,m*s);
Yf = zeros(NN,l*s);
Yp = zeros(NN,l*s);
for i = (1:s)
    if m > 0
        Up(:,(i-1)*m+1:i*m) = u(i:NN+i-1,:);
        Uf(:,(i-1)*m+1:i*m) = u(s+i:NN+s+i-1,:);
    end
    Yp(:,(i-1)*l+1:i*l) = y(i:NN+i-1,:);
    Yf(:,(i-1)*l+1:i*l) = y(s+i:NN+s+i-1,:);
end
Rnew = triu(qr([Rold'; [Uf Up Yp Yf]]))';
Rnew = Rnew(1:2*(m+l)*s,1:2*(m+l)*s);

R32 = Rnew((2*m+l)*s+1:2*(m+l)*s,m*s+1:(2*m+l)*s);
[Un,Sn] = svd(R32);
Sn = diag(Sn);
Sn = Sn(1:s);

Rnew(1,2) = m;
Rnew(1,3) = l;
Rnew(1,4) = s;
Rnew((2*m+l)*s+1:(2*m+2*l)*s,(2*m+2*l)*s+1:(2*m+3*l)*s) = Un;







