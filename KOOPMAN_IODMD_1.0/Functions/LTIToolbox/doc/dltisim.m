
%% DLTISIM
% Simulates a discrete-time LTI state-space system.

%% Syntax
% |y = dltisim(A,B,C,D,u)|
%%
% |[y,x] = dltisim(A,B,C,D,u,x0)|

%% Description
% This function simulates a discrete-time LTI state-space system. The model
% structure is the following:
%
% $$ x(k+1) = Ax(k) + Bu(k) $$
%
% $$ y(k) = Cx(k) + Du(k) $$
%
% An optional initial state can be specified.

%% Inputs
% |A,B,C,D| is the discrete-time LTI state-space system matrices.
%%
% |u| is an _N_ x _m_ matrix conatining _N_ samples of the _m_ inputs.
%%
% |x0| is the (optional) initial stae. an _n_ x _1_ vector.

%% Outputs
% |y| is the computed output sequence, an _N_ x _l_ matrix.
%%
% |x| is the (optional) computed state, an _N_ x _n_ matrix.

%% Algorithm
% A direct iteration of the system's state-transition equation is used to 
% obtain the state-trajectory for all time-instants. The function
% <ltiitr.html |ltiitr|> is used to this end.

%% Uses Functions
% <ltiitr.html |ltiitr|>

%% See Also
% |dlsim|, <ltiitr.html |ltiitr|>
