
%% DORDPI
% Preprocesses time-domain data for PI-MOESP subspace identification of
% discrete-time LTI state-space models. Delivers an order-estimate.

%% Syntax
% |[S,R] = dordpi(u,y,s)|
%%
% |[S,R] = dordpi(u,y,s,Rold)|

%% Description
% This function performs the initial data compression for PI-MOESP subspace
% identification based on measured input-output data [1]. In addition, it
% delivers information usuable for determining the required model order.
% The model structure is the following
%
% $$ x(k+1) = Ax(k) + Bu(k) $$
%
% $$ y(k)   = Cx(k) + Du(k) + v(k) $$
%
% Here, _v(k)_ is zero-mean noise of arbitary color, independent of the
% noise-free input _u(k)_ . Several data batches can be concatenated, as
% shown below. This function acts as a preprocessor to <dmodpi.html |dmodpi|>.

%% Inputs
% |u,y| is the measured input and output data of the system to be
% identified.
%%   
% |s| is the block-size parameter. This scalar should be _>n_.
%%  
% |Rold| is the (optional) data-matrix resulting from a previous call to
% <dordpi.html |dordpi|>.

%% Outputs
% |S| is the first _s_ singular values of the rank-deficient _R32_ matrix
% (see below).
% 
% |R| is a compressed data matrix containing information about the measured
% data, as well as information regarding the system dimensions.

%% Algorithm
% The discrete-time data compression algorithm in [1] is used. The
% following RQ-factorization is made:
%%
% <<dordpi_pic1.jpg>>
%%
% The meaning of the various matrices can be found in the cited article. A
% weighted SVD of the _R32_ matrix is made, and its left singular vectors
% are appended to the |R|-matrix. Its first _s_ singular values are
% returned in |S|.

%% Used By
% This a top-level function that is used directly by the user.

%% See Also
% <dordpo.html |dordpo|>, <dmodpo.html |dmodpo|>,
% <dmodpi.html |dmodpi|>, <dordrs.html |dordrs|>, <dmodrs.html |dmodrs|>

%% References
% [1] M. Verheagen, "Identification of the deterministic part of MIMO
% state space models given in innovations form from input-output data",
% _Automatica_, vol. 30, no. 1, pp. 61-74, 1994.