
%% CHOLICM
% Calculates a Cholesky-factor of the inverse covariance matrix (ICM) of a
% multivariable autoregressive process.

%% Syntax
% |C = cholicm(Af,Ab,Sf,Sb,N)|

%% Description
% The Inverse Covariance Matrix _S_ of a multivariable autoregressive noise
% process according to [1] is calculated. A Cholesky factor _C_ is returned
% such that _C'C_ = _S_
% 
% The noise model contains a causal and an auticausal part, both of which
% describe the actual noise _v(k)_. If _e(k)_ is a Gaussian white
% innovation, the model is given by: 
%
% $$v(k) =
% \mathord{\buildrel{\lower3pt\hbox{$\scriptscriptstyle\rightharpoonup$}}
% \over e}(k) -
% \mathord{\buildrel{\lower3pt\hbox{$\scriptscriptstyle\rightharpoonup$}}
% \over A}_1 v(k-1) - ... -
% \mathord{\buildrel{\lower3pt\hbox{$\scriptscriptstyle\rightharpoonup$}}
% \over A}_d v(k-d)$$
%
% $$v(k) =
% \mathord{\buildrel{\lower3pt\hbox{$\scriptscriptstyle\leftharpoonup$}}
% \over e}(k) -
% \mathord{\buildrel{\lower3pt\hbox{$\scriptscriptstyle\leftharpoonup$}}
% \over A}_1 v(k+1) - ... -
% \mathord{\buildrel{\lower3pt\hbox{$\scriptscriptstyle\leftharpoonup$}}
% \over A}_d v(k+d)$$
%
% Right arrow denote the causal (Forward) components while left arrows
% denote the anti-causal (Backward) ones. 
%
% The multivariable noise model can be obtained using the <destmar.html
% |destmar|> function

%% Inputs
% |Af| is an _l_ x _ld_ matrix containing the causal part of the noise
% process. |Af = [Af1,Af2,...,Afd]|

%%
% |Ab| is an _l_ x _ld_ matrix containing the anti-causal part of the noise
% process. |Ab = [Ab1,Ab2,...,Abd]|

%%
% |Sf| is an _l_ x _l_ matrix describing the covariance |E[ef ef']| 

%%
% |Sb| is an _l_ x _l_ matrix describing the covariance |E[eb eb']| 

%%
% |N| is the number of samples

%% Outputs
% C is a Cholesky factor of the ICM. This matrix is stored in LAPACK/BLAS
% band-storage; its size is _(d_ + _1)l_ x _N_, and the bottom row contains
% the diagonal of _C_. The row above contains a zero, and then the first
% superdiagonal of _C_. The row above contains two zeros, and then the
% second superdiagonal, etc. The top row contains _((d_ + _1)l_ - _1)_
% zeros, and then the _((d_ + _1)l_ - _1)_ th superdiagonal.


%% Limitations
% A covariance matrix of a stationary process is always positive definite.
% However, it is very well possible to specify filter coefficients |Af|,
% |Ab| and covariances |Sf| and |Sb| such that the theoretical ICM
% calculated per [1] is not positive definite. In such cases, no Cholesky
% factor can be calculated, and an identity matrix will be returned along
% with a warning message. The filter should be checked and adjusted in
% these cases.

%% Algorithm
% The upper-triangular block-band of a sparse banded inverse covariance
% matrix according to [1] is filled. A direct sparse Cholesky factorization
% is subsequently performed using MATLAB's internal |chol| function.

%% Used By
% <doptlti.html |doptli|>

%% See Also
% <doptlti.html |doptli|>, <destmar.html |destmar|>

%% References
% [1] B. Davis, _Parameter Estimation in Nonlinear Dynamical Systems with
% Correlated Noise._ PhD thesis, Universite Catholique de Louvain-La-Neuve,
% Belgium, 2001.