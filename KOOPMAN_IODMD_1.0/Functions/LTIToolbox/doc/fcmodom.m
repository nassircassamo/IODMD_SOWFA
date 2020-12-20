
%% FCMODOM
% Estimates the _A_ and _C_ matrix in a continuous-time state-space model
% from frequency response function (FRF) data that was preprocessed by
% <fcordom.html |fcordom|>.

%% Syntax
% |[A,C] = fcmodom(R,n)|

%% Description
% This function estimates the _A_ and _C_ matrices corresponding to an _n_
% th order discrete-time LTI state-space model. The compressed data matrix
% _R_ from the preprocessor function <fcordom.html |fcordom|> is used to
% this end.

%% Inputs
% |R| is a compressed data matrix containing information about the measured
% data, as well as information regarding the system dimensions.
%% 
% |n| is the desired model order _n_.

%% Outputs
% |A| is the state-space model's _A_ matrix.
%%
% |C| is the state-space model's _C_ matrix.

%% Algorithm
% The data matrix obtained with <fcordom.html |fcordom|> contains the
% weighted left singular vectors of a matrix similar to the |R22| matrix
% (see <fdordom.html |fdordom|>). Unlike in the discrete-time case, the
% first _n_ of these vectors do not form a direct estimate _Os_ of the
% extended observability matrix. Rather, a generalized matrix _Og_ is
% estimated because of the Forsythe-recursions in the data-compression
% step. The _Ahat_ and _Chat_ estimates are extracted such that this
% generalized shift-structure is taken into account [1].

%% Used By
% This a top-level function that is used directly by the user.

%% See Also
% <fcordom.html |fcordom|>, <fdmodom.html |fdmodom|>

%% References
% [1] R. Pintelon, "Frequency domain subspace system identfication using
% non-parametric noise models", in _Proceedings of the 40th IEEE Conference
% on Decision and Control_, Orlando, Florida, pp. 3916-3921, Dec 2001.