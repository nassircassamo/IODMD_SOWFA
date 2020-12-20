
%% FCMODOM
% Estimates the _A_ and _C_ matrix in a discrete-time state-space model
% from frequency response function (FRF) data that was preprocessed by
% <fdordom.html |fdordom|>.

%% Syntax
% |[A,C] = fdmodom(R,n)|

%% Description
% This function estimates the _A_ and _C_ matrices corresponding to an _n_
% th order discrete-time LTI state-space model. The compressed data matrix
% _R_ from the preprocessor function <fdordom.html |fdordom|> is used to
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
% (see <fdordom.html |fdordom|>). The first _n_ of these vectors form an
% estimate _Os_ of the system's extended observability matrix:
%%
% <<extobs.jpg>>
%%
% The estimates Ahat and Chat are obtained by linear regression:
%
% $$ \hat{C} = \hat{\mathcal{O}}_s(1:\ell,:) $$
%
% $$ \hat{A} = \hat{\mathcal{O}}_s(1:(s-1)\ell,:)^\dagger
% \hat{\mathcal{O}}_s(\ell+1:s\ell,:) $$
%

%% Used By
% This a top-level function that is used directly by the user.

%% See Also
% <fdordom.html |fdordom|>, <fcmodm.html |fcmodom|>
