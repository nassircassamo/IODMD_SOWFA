
%% FAC2B
% Estimates the _B_ and _D_ matrices in discrete-time and continuous-time
% state-space models from frequency response function (FRF) data.

%% Syntax
% |[B,D] = fac2b(A,C,H,w)|
%%
% |[B,D] = fac2b(A,C,H1,w1,...,Hp,wp)|

%% Description
% This function estimates the _B_ and _D_ matrices corresponding to a
% discrete-time or continuous-time LTI state-space model. The estimate is
% based on the measured frequency response function (FRF) data, and on the
% _A_ and _C_ matrices, which are possibly estimated using <fdmodom.html
% |fdmodom|> or <fcmodom.html |fcmodom|>.  Several data batches can be
% concatenated, though this is possible for discrete-time models only.

%% Inputs
% |A| is the state-space model's _A_ matrix.
%%
% |C| is the state-space model's _C_ matrix.
%%
% |H| is the measured frequency response function (FRF). This should be a
% matrix which follows the convention of MATLAB 6; it should be _l_ x _m_ x
% _N_ in which _H(:,:,i)_ contains the complex FRF at the _i_ th complex
% frequency.
%%          
% |w| is the vector of complex frequencies at which the FRF is measured.
% Although the function can operate using arbitrary complex frequencies,
% the following two choices are rather standard for discrete and continuous
% time models respectively:
%
% $$ \mathtt{w} = e^{j\omega} $$
%
% $$ \mathtt{w} = j\omega $$
%
% For discrete-time models, multiple data batches can be concatenated by
% appending additional |H|, |w| pairs to the parameter list.

%% Outputs
% |B| is the state-space model's _B_ matrix.
%%
% |D| is the state-space model's _D_ matrix.
%%
% |R| is a compressed data matrix that can be used to concatenate another
% data batch in a subsequent call to <fac2bd.html |fac2bd|> (discrete-time
% models only).

%% Algorithm
% Estimating _B_ and _D_ from the frequency response function (FRF) data
% and _A_ and _C_ is a linear regression [1]: 
%%
% <<fac2bd_pic1.jpg>> 
%%
% The regression matrix |Phi| and data matrix |theta| are given by:
%%
% <<fac2b_pic2.jpg>> 
%%
% The function <ltiitr.html |ltifrf|> is used to efficiently fill the
% regression matrix |Phi|.

%% Used By
% This a top-level function that is used directly by the user.

%% Uses Functions
% <ltiitr.html |ltifrf|>

%% See Also
% <fac2b.html |fac2b|>, <fdmodom.html |fdmodom|>, <fcmodom.html
% |fcmodom|>, <ltiitr.html |ltifrf|>

%% References
% [1] T. McKElvey, H. Akcay, and L. Ljung, "Subspace-based multivariable
% system identification from frequency response data", _IEEE Transactions
% on Automatic Control_, vol. 41, pp. 960-979, July 1996.