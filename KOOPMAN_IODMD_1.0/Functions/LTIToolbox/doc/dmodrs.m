
%% DMODRS
% Estimates the _A_ and _C_ matrix in a discrete-time state-space model 
% from time-domain data that was preprocessed by <dordrs.html |dordrs|>.

%% Syntax
% |[A,C] = dmodrs(R)|

%% Description
% This function estimates the _A_ and _C_ matrices corresponding to an _n_
% th order discrete-time LTI state-space model. The compressed data matrix
% |R| from the preprocessor function <dordrs.html |dordrs|> is used to this
% end. As _n_ is determined from the _x_ matrix that was passed to
% <dordrs.html |dordrs|>, it does not have to be specified here.

%% Inputs
% |R| is a compressed data matrix containing information about the measured
% data, as well as information regarding the system dimensions.

%% Outputs
% |A| is the state-space model's _A_ matrix.
%%
% |C| is the state-space model's _C_ matrix.

%% Algorithm
% The data matrix obtained with <dordrs.html |dordrs|> contains the
% weighted left singular vectors of the _R_ matrix. The first _n_ of these
% vectors form an estimate _Os_ of the system's extended observability
% matrix:
%%
% <<extobs.jpg>>
%%
% The estimates |Ahat| and |Chat| are obtained by linear regression:
%
% $$ \hat{C} = \hat{\mathcal{O}}_s(1:\ell,:) $$
%
% $$ \hat{A} = \hat{\mathcal{O}}_s(1:(s-1)\ell,:)^\dagger
% \hat{\mathcal{O}}_s(\ell+1:s\ell,:) $$
%

%% Used By
% This a top-level function that is used directly by the user.

%% See Also
% <dordrs.html |dordrs|>, <dordpo.html |dordpo|>, <dmodpo.html |dmodpo|>,
% <dordpi.html |dordpi|>, <dmodpi.html |dmodpi|>