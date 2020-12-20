function [A,C] = fdmodom(R,n) 
%FDMODOM   Estimates the A and C matrices of an discrete-time LTI 
%          state space model form using the result of the preprocessor 
%          routine fdordom. 
%
% General model structure: 
%                          -1 
%           H = C (w I - A)   B  +  D 
% 
% Syntax: 
%          [A,C] = fdmodom(R,n); 
% 
% Input: 
%   R      Data structure obtained from fdordom, containing the 
%          triangular factor and additional information 
%          (such as I/O dimension etc.) 
%   n      Order of system to be estimated. 
% 
% Output: 
%   A,C    Estimated system matrices. 
% 
% See also: fdordom, fcmodom 
% 
 
% Written by Niek Bergboer, 2001 
% Revised by Ivo Houtzager, 2007
% Copyright (c) 1995-2007, Delft Center of Systems and Control
 
if nargin ~= 2
    error('FDMODOM requires 2 input arguments.');
end
if size(R,2)<2 || size(R,2)<4,
    error('The data matrix is too small to contain any information.');
end

m = R(1,2);
l = R(1,3);
s = R(1,4);
wscale = R(2,3);
Un = R(s*m+1:s*(m+l),s*(m+l)+1:s*(m+2*l));

if m < 1
    error('FDMODOM requires an input.');
end
if l < 1
    error('FDMODOM requires an output.');
end
if s < 2
    error('Illegal value for the block-size in R.');
end
if n > s*l
    error('The order is too large, it must be smaller than s l.');
end

if wscale ~= 0
    error('The data matrix belongs to a continuous-time problem.');
end
if size(R,1) ~= s*(m+l) || size(R,2) ~= s*(m+2*l),
    error('The data matrix has is incorrectly sized.');
end

A = Un(1:(s-1)*l,1:n)\ Un(l+1:s*l,1:n);
C = Un(1:l,1:n);




