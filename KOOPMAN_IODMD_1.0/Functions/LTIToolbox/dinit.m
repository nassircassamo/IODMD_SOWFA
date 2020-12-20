function x0 = dinit(A,B,C,D,u,y)
%DINIT	    Estimate the initial state, given the estimated system
%           matrices, and a  set of input/output data.
%
% Syntax:
%           x0=dinit(A,B,C,D,u,y);
%
% Input:
%   A,B,C,D System matrices.
%   u, y    The input respectively output data of the system to be 
%           identified.
%
% Output:
%   x0      Estimated initial state.
% 
% See also: dac2bd, dac2b

% Bert Haverkamp, April 1996
% Revised by Ivo Houtzager, 2007
% Copyright (c) 1996-2007, Delft Center of Systems and Control 

% Check number of arguments
if nargin < 6
    error('DFUNLTI requires at least six input arguments.');
end
 if size(y,2)>size(y,1)
     y = y';
 end
if size(u,2)>size(u,1)
    u = u';
end

error(abcdchk(A,B,C,D));
n = size(A,1);
l = size(C,1);
N = size(y,1);

if (~(size(u,1)==N) && ~isempty(u))
    error('Input and output should have the same length.')
end

%Gamma
temp = zeros(n*l,N);
Gamma = C;
An = A;
for i = 1:ceil(log(N)/log(2)),
    Gamma(2^(i-1)*l+1:2^i*l,:) = Gamma*An;
    An = An*An;
end
Gamma = Gamma(1:N*l,:);
temp(:) = Gamma';
for j = 1:l
    Gamma(N*(j-1)+1:N*j,:) = temp(n*(j-1)+1:n*j,:)';
end

if ~isempty(u)
    x = ltiitr(A,B,u,[],zeros(n,1));
    ye = x*C.'+u*D' ;
    Y = zeros(N*l,1);
    Y(:) = y - ye;
else
    Y(:) = y;
end
% equation to solve:  Y=Gamma*x0
x0 = Gamma\Y;

