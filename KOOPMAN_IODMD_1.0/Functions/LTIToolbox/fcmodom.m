function [A,C] = fcmodom(R,n) 
%FCMODOM   Estimates the A and C matrices of an continuous-time LTI 
%          state space model form using the result of the preprocessor 
%          routine fcordom. 
%
% General model structure: 
%                          -1 
%           H = C (w I - A)   B  +  D 
% 
% Syntax: 
%          [A,C] = fcmodom(R,n); 
% 
% Input: 
%   R      Data structure obtained from fcordom, containing the 
%          triangular factor and additional information 
%          (such as I/O dimension etc.) 
%   n      Order of system to be estimated. 
% 
% Output: 
%   A,C    Estimated system matrices. 
% 
% See also: fcordom, fdmodom 
% 
 
% Written by Niek Bergboer, 2001 
% Revised by Ivo Houtzager, 2007
% Copyright (c) 2001-2007, Delft Center of Systems and Control 

if nargin ~= 2
    error('FCMODOM requires 2 input arguments.');
end
if size(R,2)<2 || size(R,2)<4,
    error('The data matrix is too small to contain any information.');
end

m = R(1,2);
l = R(1,3);
s = R(1,4);
wscale = R(2,3);
Un = R(s * m+1:s * (m+l),s * (m+l)+1:s * (m+2 * l));

if m < 1
    error('FCMODOM requires an input.');
end
if l < 1
    error('FCMODOM requires an output.');
end
if s < 2
    error('Illegal value for the block-size in R.');
end 
if n > s*l
    error('The order is too large, it must be smaller than s l.');
end

if wscale == 0,
    error('The data matrix belongs to a discrete-time problem.');
end

if size(R,1)~=s*(m+l) || size(R,2)~=s*(m+2*l)
    error('The data matrix has is incorrectly sized.');
end

Zo = R(1:s,s*(m+l)+1:s*(m+l)+l);
Os = Un(:,1:n);
dZ = reshape(Zo',s*l,1);
D1 = spdiags(1./dZ(l+1:s*l,1),0,(s-1)*l,(s-1)*l);
D2 = diag(dZ(l+1:(s-1)*l)./dZ(2*l+1:s*l,1));
b = [zeros(l,n); D2*Os(1:(s-2)*l,:)];

A = (D1*Os(1:(s-1)*l,:))\(Os(l+1:s*l,:)-b);
C = full(diag(dZ(1:l,:)))*Os(1:l,:);
A = wscale*A;
C = wscale*C;
