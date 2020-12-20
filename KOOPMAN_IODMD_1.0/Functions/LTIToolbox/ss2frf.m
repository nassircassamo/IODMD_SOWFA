function H=ss2frf(A,B,C,D,w)
%SS2FRF     This function calculates the frequency-response function
%           of a LTI state-space system that has the following model
%           structure:
%              . 
%              x(t) = Ax(t) + Bu(t)
%              y(t) = Cx(t) + Du(t)
%              x(0) = x0 
%           or:
%              x(k+1) = Ax(k) + Bu(k)
%                y(k) = Cx(k) + Du(k)
% 
% Syntax: 
%           H = ss2frf(A,B,C,D,w)
%
% Input: 
%  A,B,C,D  System matrices describing the state-space system. 
%           If there is no D-matrix, pass an empty matrix.
%  w        Vector of complex frequencues at which the frequency
%           response function is to be calculated. Typically, this
%           is (j.omega) for continuous-time systems and
%           exp(j.omega) for discrete-time systems.
% 
% Output: 
%   H       The measured frequency response function (FRF). This is
%           a matrix that follows the MATLAB 5 FRF Convention:
%           It is a 3D-array of size l x m x N, in which H(:,:,j)
%           is the FRF at the j-th complex frequency.
 
% Niek Bergboer, 2001
% Revised by Ivo Houtzager, 2007
% Copyright (c) 2001-2007, Delft Center of Systems and Control 

if nargin < 5
    error('SS2FRF requires five input arguments');
end
if size(w,2)>size(w,1),
    w = w.';
end

% Get data from ltifrf
H = ltifrf(A,B,C,D,[],w,[]);

