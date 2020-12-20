function [S,W,Rnew] = dordpe(u,y,s,Rold)
%DORDPE     Delivers information about the persistency of excitation of the
%           input signal u(k) and acts as a pre-processor for DSIGPE. The
%           latter actually estimates a excitation signal r(k), which
%           improves the persistency of excitation.
%
% Syntax:
%           [S,W,Rnew]=dordpe(u,y,s); [S,W,Rnew]=dordpe(u,y,s,Rold);
%
% Input:
%   u, y    The input respectively output data of the system to be 
%           identified. DORDPE does not handle empty inputs.
%   s       The dimension parameter that determines the number
%           of block rows in the processed Hankel matrices. This parameter
%           should be chosen larger than the expected system order. The
%           optimal value has to be found by trial and error. Generally
%           twice as large is a good starting value.
%   Rold    Should not be there on the first call of
%           the routine. When a second (or third, ...) i/o data batch is
%           processed a present R matrix is used to store the information
%           from these different data batches. 
%
% Output:
%   S       Singular values bearing information on the persistency of
%           excitation of the input signal u(k). 
%   W       Orthogonal matrix containing the excited directions of the 
%           input signal u(k).
%   Rnew    The compressed lower triangular factor with additional
%           information (such as i/o dimension etc.) stored in the zeros.
%
% See also: DSIGPE

% References:
%   R. Hallouzi, "Multiple-Model Based Diagnosis for Adaptive 
%   Fault-Tolerant Control" PhD thesis, April 17, 2008.

% Revised by Ivo Houtzager, 2008
% Copyright (c) 2008, Delft Center of Systems and Control 

% check number of arguments
if nargin < 3
    error('DORDPE requires at least three input arguments.');
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

% does not handle empty inputs
if m == 0
    error('DORDPE requires an input')
end
if l == 0
    error('DORDPE requires an output')
end
if s < 2
    error('s must be at least 2')
end
if (~(size(u,1) == N) && ~isempty(u))
    error('Input and output should have same length')
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
    Up(:,(i-1)*m+1:i*m) = u(i:NN+i-1,:);
    Uf(:,(i-1)*m+1:i*m) = u(s+i:NN+s+i-1,:);
    Yp(:,(i-1)*l+1:i*l) = y(i:NN+i-1,:);
    Yf(:,(i-1)*l+1:i*l) = y(s+i:NN+s+i-1,:);
end
Rnew = triu(qr([Rold'; [Yp Up Uf Yf]]))';
Rnew = Rnew(1:2*(m+l)*s,1:2*(m+l)*s);

% return W and S
R2122 = Rnew((m+l)*s+1:(2*m+l)*s,1:(2*m+l)*s);
[W,S] = svd(R2122);
S = diag(S);
S = S(1:s);