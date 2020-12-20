function y = prbn(N,rate)
%PRBN       Produce a binary sequence, with values 0 and 1. 
%           The chance of switching from level is given by the
%           parameter 'rate'. Rate=0 gives a constant value
%           0. Rate=1 gives an signal that changes
%           constantly between 0 and 1. Any value in
%           between results in a random binary sequence.
%           This kind of testsignal has been described in 
%           H.J.A.F Tulleken, 'Generalized binary noise test-signal 
%           concept for improved identification-experiment design',
%           Automatica 1990 vol 26.
% Syntax:
%           y=prbn(N,rate);
% Input:
%  N        Number of points.
%  rate     Chance of the signal changing level at every
%           sample. Default is 0.5.
%  
% Output:
%  y        Random binary noise.

% Bert Haverkamp, April 1996
% Revised by Ivo Houtzager, 2007
% Copyright (c) 1996-2007, Delft Center of Systems and Control 

if nargin < 1
  error('PRBN requires at least one input parameters')
end

if N < 1
  error('number of points should be larger than zero')
end  
if nargin == 1
  rate = 0.5;
end

if rate<  0 || rate>1
  error('value of rate  should be between 0 and 1')
end

dy = floor(rand(N,1)+rate);
y = rem(cumsum(dy),2);







