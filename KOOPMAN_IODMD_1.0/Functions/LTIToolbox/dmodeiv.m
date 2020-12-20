function [A,C] = dmodeiv(R,n)
%DMODEIV    Estimates the A and C matrices of a LTI state 
%           space model in innovations form using the output 
%           of the dordeiv routine. 
% 
% Model structure for the errors-in-variables problem: 
%              x(k+1) = Ax(k) + Bu~(k) + f(k)
%              y~(k)  = Cx(k) + Du~(k) 
%           with measurements 
%              u(k) = u~(k) + w(k)
%              y(k) = y~(k) + v(k)
%           where f(k), w(k) and v(k) are zero-mean white noise sequences
%           independent of the input u~(j) for k >= j. 
%
%           Note that the plant can be operated under either open-loop
%           or closed-loop. 
%
% Syntax:
% 	    [A,C] = dmodeiv(R,n);
% 					
% Input:
%   R       Triangular factor from order detection step dordeiv
%   n       Order of system to be estimated
% 			
% Output:
%   A, C    Estimated system matrices
% 					
% See also: dordeiv, dac2bd_eiv, dac2b_eiv, dac2bd_cl, dac2b_cl 
%
% Reference: C.T. Chou and M. Verhaegen 
%            Subspace algorithms fpr the identification of
%            multivariable dynamic errors-in-variables models
%            Automatica vol 33, no 10, pp. 1857-1869, October, 1997
% 

% C.T. Chou, Oct 1997
% Revised by Ivo Houtzager, 2009
% Copyright (c) 1997-2009, Delft Center of Systems and Control 

%m = R(1,2);
L = R(1,3);
s = R(2,3);

sL = s*L;

% find A and C 
[Un,Sn] = svd(tril(R));
A = Un(1:sL-L,1:n)\Un(L+1:sL,1:n);
C = Un(1:L,1:n);

