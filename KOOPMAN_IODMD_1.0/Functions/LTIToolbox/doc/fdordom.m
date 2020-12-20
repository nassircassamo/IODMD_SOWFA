
%% FDORDOM
% Preprocesses frequency-domain data for frequency-domain subspace
% identification of discrete-time state-space models. Delivers an
% order-estimate.

%% Syntax
% |[S,R] = fdordom(H,w,s)|
%%
% |[S,R] = fdordom(H,w,s,Rold)|

%% Description
% This function performs the initial data compression for discrete-time
% subspace identification based on measured frequency reponse function
% (FRF) data. In addition, it delivers information usuable for determining
% the required model order. The model structure is the following: 
% 
% $$ \dot{x}(t) = Ax(t) + Bu(t) $$
%
% $$ y(t) = Cx(t) + Du(t) $$
% 
% Several data batches can be concatenated, as shown below. This function
% acts as a preprocessor to <fdmodom.html |fdmodom|>.

%% Inputs
% |H| is the measured frequency response function (FRF). This should be a
% matrix which follows the convention of MATLAB 6; it should be _l_ x _m_
% x _N_ in which _H(:,:,i)_ contains the complex FRF at the _i_ th complex
% frequency.
%% 
% |w| is the vector of complex frequencies at which the FRF is measured:
% 
% $$ \mathtt{w} = e^{j\omega} $$
% 
%%
% |s| is the block-size parameter. This scalar should be _>n_.
%%
% |Rold| is the (optional) data-matrix resulting from a previous call to
% <fdordom.html |fdordom|>.

%% Outputs
% |S| is the first _s_ singular values of the rank-deficient _R22_ matrix
% (see below).
% 
% |R| is a compressed data matrix containing information about the measured
% data, as well as information regarding the system dimensions.

%% Remarks
% The MEX-implementation may generate the following warning:
% 
% |Cholesky-factorization failed; falling back on QR-factorization.|
% 
% This implies that the fast Cholesky-algorithm failed. The function has
% automatically fallen back onto a slower QR-algorithm. Results from
% |fdordom| can be used without problems if this warning appears.

%% Algorithm
% The discrete-time data compression algorithm in [1] is used. In the
% M-file implementation, the following RQ-factorization is made:
%%
% <<fdordom_pic1.jpg>>
%%
% The meaning of the various matrices can be found in the cited article. In
% the MEX-implementation, the following Cholesky-factorization is attempted
% first:
%%
% <<fdordom_pic2.jpg>>
%%
% If this factorization fails, the algorithm falls back on the above
% RQ-factorization. In all cases, a weighted SVD of the _R22_ matrix is
% made, and its left singular vectors are appended to the |R|-matrix. Its
% first _s_ singular values are returned in |S|.

%% Used By
% This a top-level function that is used directly by the user.

%% Uses Functions
% LAPACK-functions |DPOTRF|, |DGEQRF|, |DGESVD|, |DTRTRS|.
%
% BLAS-functions |DTRMM|.
%
% (All built into the executable)
  
%% See Also
% <fdmodom.html |fdmodom|>, <fcordom.html |fcordom|>

%% References
% [1] T. McKElvey, H. Akcay, and L. Ljung, "Subspace-based multivariable
% system identification from frequency response data", _IEEE Transactions
% on Automatic Control_, vol. 41, pp. 960-979, July 1996.