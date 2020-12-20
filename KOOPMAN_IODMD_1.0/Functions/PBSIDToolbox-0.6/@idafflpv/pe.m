function [e,x0] = pe(sys,u,y,t,p,type)
%PE Computes the prediction errors of an IDAFFLPV.
%  [E,X0] = PE(M,U,Y,T,MU) returns the residue response and initial state
%  of the IDAFFLPV model M to the input and scheduling signal described by
%  U, Y, MU and T. The time vector T consists of regularly spaced time
%  samples, U, Y, and MU is are matrices with as many columns as inputs and
%  scheduling variables and whose i-th row specifies the input value at
%  time T(i). For discrete-time models, U, Y, and MU should be sampled at
%  the same rate as M.
%
%  [E,X0] = PE(M,U,Y,T,MU,'Type') specifies the type of LPV predictor. 'K'
%  specifies that K is not dependent on scheduling, or 'CD' specifies
%  that C and D is not dependent on scheduling. Default is Type = 'CD'.
%
%  NOTE: This is only possible if K or C and D are not dependent on the
%  scheduling sequence.

% Assign values to unspecified parameters
if nargin < 6 || isempty(type)
    type = 'CD';
end
x0 = findstates(sys,u,y,t,p,type);
e = resid(sys,u,y,t,p,x0,type);
end
