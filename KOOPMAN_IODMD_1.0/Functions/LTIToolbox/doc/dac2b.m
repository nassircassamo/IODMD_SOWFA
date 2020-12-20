
%% DAC2B
% Estimates the _B_ matrix in discrete-time LTI state-space models from
% input-output measurements.

%% Syntax
% |B = dac2bd(A,C,u,y)|
%%
% |B = dac2b(A,C,u1,y1,...,up,yp)|

%% Description
% This function estimates the _B_ matrix corresponding to a discrete-time
% LTI state-space model. The estimate is based on the measured input-output
% data sequences, and on the _A_ and _C_ matrices, which are possibly
% estimated using <dmodpo.html |dmodpo|>, <dmodpi.html |dmodpi|> or
% <dmodrs.html |dmodrs|>. The _D_ matrix is assumed to be zero. Several data
% batches can be concatenated.

%% Inputs
% |A| is the state-space model's _A_ matrix.
%%
% |C| is the state-space model's _C_ matrix.
%%
% |u,y| is the measured input-output data from the system to be identified.
%
% Multiple data batches can be specified by appending additional |u,y|
% pairs to the parameter list.

%% Outputs
% |B| is the state-space model's _B_ matrix.

%% Algorithm
% Estimating |B| and the initial state |x0| from input-output data and |A|
% and |C| is a linear regression [1]:
%%
% <<dac2b_pic1.jpg>>
%%
% The regression matrix |Phi| and data matrix |theta| are given by:
%%
% <<dac2b_pic2.jpg>>
%%
% The function <ltiitr.html |ltiitr|> is used to efficiently fill the
% regression matrix |Phi|.

%% Used By
% This a top-level function that is used directly by the user.

%% Uses Functions
% <ltiitr.html |ltiitr|>

%% See Also
% <dac2bd.html |dac2bd|>, <dmodpo.html |dmodpo|>, <dmodpi.html |dmodpi|>,
% <dmodrs.html |dmodrs|>, <ltiitr.html |ltiitr|>

%% References
% [1] B. Haverkamp, _Subspace Method Identification, Theory and Practice._
% PhD thesis, Delft University of Technology, Delft, The Netherlands, 2000.