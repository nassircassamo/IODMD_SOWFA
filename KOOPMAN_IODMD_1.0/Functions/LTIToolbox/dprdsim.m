function [y,x] = dprdsim(A,B,C,D,K,u,y,x0) 
%DPRDSIM  Simulates a predictor discrete-time LTI state-space system. 
% 
% The model structure is the following: 
%           x(k+1) = (A-KC)x(k) + (B-KD)u(k) + Ky(k) 
%           y(k)   = Cx(k) + Du(k)
% 
% Syntax: 
%               y = dprdsim(A,B,C,D,K,u,y) 
%           [y,x] = dprdsim(A,B,C,D,K,u,y,x0) 
% 
% Input: 
%   A,B,C,D The state-space system matrices. 
%   K       Kalman gain
%   u       The measured input sequence to the system.
%   y       The measured output sequence of the system 
%   x0      (optional) The initial state. 
% 
% Output: 
%   y       The estimated output sequence. 
%   x       (optional) The estimated state sequence. 
% 
% See also: 
%           dltisim, dlsim 
% 

% Revised by Ivo Houtzager, 2008
% Copyright (c) 1996-2008, Delft Center of Systems and Control 

% Check number of arguments
if ~(nargin == 7 || nargin == 8)
    error('DPRDSIM requires at 5 or 6 input arguments.');
end
if nargin == 7
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
N = size(u,2);

% Simulation
x = zeros(size(A,1),N);
for k = 1:N
    x(:,k) = x0;
    x0 = (A-K*C)*x0 + (B-K*D)*u(:,k) + K*y(:,k);
end
y = C*x + D*u;