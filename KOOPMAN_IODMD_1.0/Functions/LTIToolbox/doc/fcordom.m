
%% FCORDOM
% Preprocesses frequency-domain data for frequency-domain subspace
% identification of continuous-time state-space models.

%% Syntax
% |[S,R] = fcordom(H,w,s)|

%% Description
% This function performs the initial data compression for continuous-time
% subspace identification based on measured frequency reponse function
% (FRF) data. In addition, it delivers information usuable for determining
% the required model order. The model structure is the following: 
% 
% $$ \dot{x}(t) = Ax(t) + Bu(t) $$
%
% $$ y(t) = Cx(t) + Du(t) $$
% 
% This function acts as a preprocessor to <fcmodom.html |fcmodom|>. Unlike
% in the discrete-time case, concatenating multiple data batches are not
% supported.

%% Inputs
% |H| is the measured frequency response function (FRF). This should be a
% matrix which follows the convention of MATLAB 6; it should be _l_ x _m_
% x _N_ in which _H(:,:,i)_ contains the complex FRF at the _i_ th complex
% frequency.
%% 
% |w| is the vector of complex frequencies at which the FRF is measured:
% 
% $$ \mathtt{w} = j\omega $$
% 
%%
% |s| is the block-size parameter. This scalar should be _>n_.

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
% <fcordom.html |fcordom|> can be used without problems if this warning
% appears.

%% Algorithm
% The continuous-time data compression algorithm in [1] is used. The same
% factorizations as in the discrete-time function <fdordom.html |fdordom|>
% are used. However, the |W| and |G| matrices are formed by
% Forsythe-recursions to prevents ill-conditioning because the complex
% frequencies are not of unit magnitude [1,2].
% 
% A weighted SVD of the |R22| matrix is made, and its left singular vectors
% are appended to the |R|-matrix. Its first _s_ singular values are
% returned in |S|.

%% Used By
% This a top-level function that is used directly by the user.

%% Uses Functions
% LAPACK-functions |DPOTRF|, |DGEQRF|, |DGESVD|, |DTRTRS|.
%
% BLAS-functions |DTRMM| and |DGEMM|.
%
% (All built into the executable)
  
%% See Also
% <fcmodom.html |fcmodom|>, <fdordom.html |fdordom|>

%% References
% [1] P. van Overschee and B. De Moor, "Continuous-time frequency domain
% subspace system identification", _Signal processing_, vol.52, no. 2, pp.
% 179-194, 1996.
%
% [2] R. Pintelon, "Frequency domain subspace system identfication using
% non-parametric noise models", in _Proceedings of the 40th IEEE Conference
% on Decision and Control_, Orlando, Florida, pp. 3916-3921, Dec 2001.