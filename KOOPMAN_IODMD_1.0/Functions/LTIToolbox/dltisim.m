function [y,x] = dltisim(A,B,C,D,u,x0) 
%DLTISIM    Simulates a discrete-time LTI state-space system. 
%           The model structure is the following: 
% 
%              x(k+1) = Ax(k) + Bu(k) + Ke(k) 
%              y(k)   = Cx(k) + Du(k) +  e(k) 
% 
% Syntax: 
%               y = dltisim(A,B,C,D,u) 
%           [y,x] = dltisim(A,B,C,D,u,x0) 
% 
% Input: 
%   A,B,C,D The state-space system matrices. 
%   u       The input sequence to the system. 
%   x0      (optional) The initial state. 
% 
% Output: 
%   y       The simulated output sequence. 
%   x       (optional) The simulated state sequence. 
% 
% See also: 
%           lpvitr, dlsim 
% 

% Revised by Ivo Houtzager, 2007
% Copyright (c) 1996-2007, Delft Center of Systems and Control 

% Check number of arguments
if ~(nargin == 5 || nargin == 6)
    error('DLTISIM requires at 5 or 6 input arguments.');
end
if nargin == 5
    x0 = zeros(size(A,1),1);
end

error(abcdchk(A,B,C,D));
[mx,nx] = size(x0);
if mx ~= size(A,1)
    error('X0 vector must have as many rows as A.')
end
if nx ~= 1
    error('X0 vector must consists of 1 column.')
end

x = ltiitr(A,B,u,[],x0);
y = x*C' + u*D';




