
%% DMODPI
% Estimates the _A_ and _C_ matrix in a discrete-time state-space model 
% from time-domain data that was preprocessed by <dordpi.html |dordpi|>.

%% Syntax
% |[A,C] = dmodpi(R,n)|
% |[A,C] = dmodpi(R,n,'stable')|

%% Description
% This function estimates the _A_ and _C_ matrices corresponding to an _n_ 
% th order discrete-time LTI state-space model. The compressed data matrix
% |R| from the preprocessor function <dordpi.html |dordpi|> is used to this
% end.

%% Inputs
% |R| is a compressed data matrix containing information about the measured
% data, as well as information regarding the system dimensions.
%%
% |n| is the desired model order _n_.
%%
% |stable| estimates a stable A matrix, see[1].

%% Outputs
% |A| is the state-space model's _A_ matrix.
%%
% |C| is the state-space model's _C_ matrix.

%% Algorithm
% The data matrix obtained with <dordpi.html |dordpi|> contains the
% weighted left singular vectors of the _R32_ matrix. The first _n_ of
% these vectors form an estimate _Os_ of the system's extended
% observability matrix:
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
% <dordpo.html |dordpo|>, <dmodpo.html |dmodpo|>, <dordpi.html |dordpi|>,
% <dordrs.html |dordrs|>, <dmodrs.html |dmodrs|>

%% References
%  [1] J.M. Maciejowski, "Guaranteed Stability with Subspace Methods",
%  Submitted to Systems and Control Letters, 1994.