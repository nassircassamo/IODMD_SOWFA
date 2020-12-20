
%% DTH2SS
% Converts a parameter vector into a discrete-time LTI state-space model.

%% Syntax
% |[A,C] = dth2ss(theta,params)|
%%
% |[A,B,C] = dth2ss(theta,params)|
%%
% |[A,B,C,D] = dth2ss(theta,params)|
%%
% |[A,B,C,D,x0] = dth2ss(theta,params)|
%%
% |[A,B,C,D,x0,K] = dth2ss(theta,params)| 

%% Description
% his function converts a parameter vector that describes a continuous-time
% state space model into the state space matrices of that model.
%
% $$x(k+1) = Ax(k) + Bu(k) + Ke(k)$$
%
% $$y(k) = Cx(k) + Du(k) + e(k)$$
%

%% Inputs
% |theta| is the parameter vector describing the system.
%%
% |params| is a structure that contains the dimension parameters of the
% system, such as the order, the number of inputs and whether |D|, |x0| or
% |K| is present.
%%
% |T| is the transformation matrix between the input state space system and
% the state space system in the form described by |theta|.

%% Outputs
% |A,B,C,D| are system matrices describing the state space system. If
% |theta| does not contain  parameters for |D|, this matrix will be
% returned as an empty matrix.
%%
% |x0| is the initial state. If |theta| does not contain parameters for
% |x0|, this vector will be returned as an empty matrix.
%%
% |K| is the Kalman gain. If |theta| does not contain parameters for |K|,
% this vector will be returned as an empty matrix.

%% Remarks
% This function is based on the SMI Toolbox 2.0 function |dth2ss|,
% copyright 1996 Johan Bruls. Support for the omission of |D|, |x0| and/or
% |K| has been added, as well as support for the full parametrization.

%% Algorithm
% The model parametrization for the output normal form and the tridiagonal
% parametrization is carried out according to [1]. The full model
% parametrization is a simple vectorization of the system matrices. In its
% most general form, the parameter vector is given by
%%
% <<param.jpg>>

%% Used By
% <doptlti.html |doptlti|>, <dfunlti.html |dfunlti|>, <foptlti.html
% |foptlti|>, <ffunlti.html |ffunlti|>

%% See Also
% <css2th.html |css2th|>, <dss2th.html |dss2th|>

%% References
% [1] B. Haverkamp, _Subspace Method Identification, Theory and Practice._
% PhD thesis, Delft University of Technology, Delft, The Netherlands, 2000.