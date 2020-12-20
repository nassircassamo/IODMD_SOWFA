function [A,C] = dmodpi(R,n,stable)
%DMODPI     Estimates the system matrices  A and C of a 
%           LTI state space model using the output of the 
%           dordpi routine.
%
% Model structure:
%               x(k+1) = Ax(k) + Bu(k)
%               y(k)   = Cx(k) + Du(k) + v(k)
%           where v(k) is zero-mean noise of arbitrary color,
%           independent of the noise-free input u(k).
%
% Syntax:
%           [A,C]=dmodpi(R,n);
%           [A,C]=dmodpi(R,n,'stable');
%
% Input:
%   R       Triangular factor from dordpi
%   n       Order of system to be estimated
%   stable  Estimates a stable A matrix.
%
% Output:
%   [A,C]   System matrices.
%
% See also: dordpi, dordpo, dmodpo

% The modpi routine corresponds to the PI scheme
% derived and analyzed in VERHAEGEN: "Subspace Model
% Identification. Part 3" Int. J. Control, Vol. 57.
%
% Michel Verhaegen 11-01-1990
% Revised by Ivo Houtzager, 2010
% Copyright (c) 1990-2010, Delft Center of Systems and Control 

% Check number of arguments
if nargin < 2
    error('DMODPI requires at least two input arguments.');
end
if (nargin < 3) || (isempty(stable))
    stable = 0;
elseif strcmpi(stable,'stable')
    stable = 1;
else
    stable = 0;
end

if (size(R,1)<2) || (size(R,2)<3)
    error('Matrix R has wrong size')
end
m = R(1,2);
l = R(1,3);
s = R(2,3);

if m < 1
    error('Illegal value for number of inputs')
end
if l < 1
    error('Illegal value for number of outputs')
end
if s < 2
    error('Illegal value  for ''s'' in R matrix')
end
if ~((size(R,1)==s*(2*m+l)) && (size(R,2)==2*s*(m+l)))
    error('R-matrix has unexpected size.')
end

if n < 1
    error('System order of zero or lower does not make sense!')
end
if  n >= s
    error('n chosen too large, it should be smaller  than ''s''.')
end

Un = R(2*m*s+1:(2*m+l)*s,(2*m+l)*s+1:2*(m+l)*s);
if stable
    un1 = Un(1:s*l,1:n);
    un2 = [Un(l+1:s*l,1:n); zeros(l,n)];
else
    un1 = Un(1:(s-1)*l,1:n);
    un2 = Un(l+1:s*l,1:n);
end
A = un1\un2;
C = un1(1:l,:);

