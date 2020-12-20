
%% DINIT
% Estimates the initial state, given estimated discrete-time state-space
% system matrices and a batch of measured input-output data.

%% Syntax
% |x0 = dinit(A,B,C,D,u,y)|

%% Description
% This function estimates the initial state for a measured input-output
% batch of a discrete-time LTI state-space model. The estimate is based on
% the measured input-output data sequences, and on the _A_, _B_, _C_ and
% _D_ matrices, which are possibly estimated using any of the subspace
% identification functions.

%% Inputs
% |A,B,C,D| is the discrete-time LTI state-space model.
%%
% |u,y| is the measured input-output data from the system to be identified.

%% Outputs
% |x0| is the estimated initial state.

%% Algorithm
% Estimating the initial state |x0| from input-output data and the system
% matrices is a linear regression [1]:
%
% $$ x_0 = \Phi^\dagger \theta $$
%
% The regression matrix |Phi| and data matrix |theta| are given by:
%%
% <<dinit_pic1.jpg>>
%%
% in which _yhatk)_ is simulated using the estimated system matrices and 
% the measured input _u(k)_. The function <ltiitr.html |ltiitr|> is used to
% efficiently calculate _yhat(k)_.

%% Used By
% This a top-level function that is used directly by the user.

%% Uses Functions
% <ltiitr.html |ltiitr|>

%% See Also
% <dac2b.html |dac2b|>, <dac2bd.html |dac2bd|>, <ltiitr.html |ltiitr|>

%% References
% [1] B. Haverkamp, _Subspace Method Identification, Theory and Practice._
% PhD thesis, Delft University of Technology, Delft, The Netherlands, 2000.