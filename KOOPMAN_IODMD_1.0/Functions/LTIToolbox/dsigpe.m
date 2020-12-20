function r = dsigpe(S,W,s,npe,rho,option)
%DSIGPE     Returns signal r(k) that has a persistency of excitation close
%           to the npe non-persistently excited directions of input u(k) as
%           possible.
% 					
% Syntax:
% 	        r = dsigpe(S,W,s,npe);
% 	        r = dsigpe(S,W,s,npe,rho);
% 	        r = dsigpe(S,W,s,npe,rho,option);
% 					
% Input:
%   S       Singular values bearing information on the persistency of
%           excitation of the input signal u(k). 
%   W       Orthogonal matrix containing the excited directions of the 
%           input signal u(k).
%   n       Order of system to be estimated.
%   npe     Number of least non-persistently excited directions that should
%           be additionally excited.
%   rho     Excitation level.
%   option  Option for weighting the excited directions: 'normal' for no 
%           weighting (default), 'weight' for using singular values S as 
%           weighting factor, 'sqrt' for using sqrt(S) as weighting factor.
% 			
% Output:
%   r       Additional excitation signal.
% 					
% See also: dordpe, dordpe_cl

% References:
%   R. Hallouzi, "Multiple-Model Based Diagnosis for Adaptive 
%   Fault-Tolerant Control" PhD thesis, April 17, 2008.

% Revised by Ivo Houtzager, 2008
% Copyright (c) 2008, Delft Center of Systems and Control 

% Check number of arguments
if nargin < 3
    error('DSIGPE requires at least three input arguments.');
end
if npe < 1
    error('NPE of zero or lower does not make sense!')
end
if rho <= 0
    error('RHO of zero or lower does not make sense!')
end

% assign default values to unspecified parameters
if (nargin < 6) || isempty(option)
    option = 'normal';
end
if (nargin < 5) || isempty(rho)
    rho = 1;
end

if strcmpi(option,'normal')
    P = W(:,end-npe+1:end);
elseif strcmpi(option,'sqrt')
    P = diag(sqrt(S))*W(:,end-npe+1:end);    
elseif strcmpi(option,'weight')
    P = diag(S)*W(:,end-npe+1:end);
else
    error('Option not recognized!')
end
P = P(:);
I = kron(ones(npe,1),eye(s));
r = (I\(rho.*P))';