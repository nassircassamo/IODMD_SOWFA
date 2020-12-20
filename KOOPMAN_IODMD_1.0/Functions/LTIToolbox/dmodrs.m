function [A,C] = dmodrs(R)
%DMODRS     Estimates the system matrices A and C of a 
%           LTI state space model using the output of the 
%           dordrs routine.
%
% Model structure:
%               x(k+1) = Ax(k) + Bu(k)
%               y(k)   = Cx(k) + Du(k) + v(k)
%           where v(k) is zero-mean noise of arbitrary color,
%           independent of the noise-free input u(k).
%
% Syntax:
%           [A,C]=dmodrs(R);
%
% Input:
%  R        Compressed data matrix.
%
% Output:
%  [A,C]    System matrices.
%
% See also: dmodpi, dordpi, dmodpo, dordpo, dordrs

% Michel Verhaegen 11-01-1990
% Revised by Ivo Houtzager, 2007
% Copyright (c) 1990-2007, Delft Center of Systems and Control 

% Check number of arguments
if nargin ~= 2
    error('DMODRS requires one input argument.');
end

if (size(R,1)<2 || size(R,2)<4)
    error('Matrix R has wrong size')
end

m = R(1,2);
l = R(1,3);
n = R(1,4);
s = R(2,3);

if m < 1
    error('Illegal value for number of inputs in R matrix')
end
if l < 1
    error('Illegal value for number of outputs in R matrix ')
end
if s < 1
    error('Illegal value  for ''s'' in R matrix')
end

if ~((size(R,1)==s*(2*m+n+l)) && (size(R,2)==s*(2*m+2*l+n)))
    error('R-matrix has unexpected size.')
end

if n >= s
    error('n chosen too large, it should be smaller than ''s''.')
end

Un = R(2*m*s+n*s+1:(2*m+n+l)*s,(2*m+n+l)*s+1:2*(m+l)*s+n*s);

un1 = Un(1:(s-1)*l,1:n);
un2 = Un(l+1:s*l,1:n);
A = un1\un2;
C = un1(1:l,:);