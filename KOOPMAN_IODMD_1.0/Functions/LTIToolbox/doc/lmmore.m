
%% LMMORE
% Performs a More-Hebden Levenberg-Marquardt optimization

%% Syntax
% |x = lmmore('func',xinit,lb,ub,options,arg2,...)|
%%
% |[x,resnorm,residual,exitflag,output,lambda,jacobian] =
% lmmore('func',xinit,lb,ub,options,arg2,...)|

%% Description
% This function is a More-Hebden implementation of the Levenberg-Marquardt
% nonlinear least-squares optimization algorithm. The function is
% interface-compatible with the |lsqnonlin|-function from the MATLAB 6
% Optimization Toolbox.

%% Inputs
% |'func'| is the cost-function that is to be used.
%%
% |xinit| is the parameter-vector's starting point in the non-linear
% optimization.
%% 
% |lb| is the lower-bound on the parameters. This value is _not_ used.
%% 
% |ub| is the upper-bound on the parameters. This value is _not_ used.
%% 
% |options| is a MATLAB 6 compatible |optimset|-structure that contains
% options for the optimization algorithm [1]. In addition, a number of
% extra fields may be present. See the Remarks section below for more
% information.
%% 
% |arg2| will be passed as second argument to the cost-function |'func'|.
% Arguments _3_ to _N_ may be appended after |arg2|.
          
%% Outputs
% |x| is the result of the optimization. The solution |x| is guaranteed to
% have an equal or smaller cost than |xinit|.
% 
% All other parameters are compatible with the MATLAB 6 |lsqnonlin|
% function.

%% Remarks
% The interface to lmmore has been made compatible with the |lsqnonlin|
% optimization function in the MATLAB 6 Optimization Toolbox.
% Note that although a lower and upper bound are given (consistent
% with |lsqnonlin|'s interface), they are _not_ used internally.
%
% This optimization implementation supports overparametrized
% cost-functions. If |options.Manifold| (not part of |optimset|'s normal
% structure) is passed and set to |'on'|, <lmmore.html |lmmore|> expects
% the cost- function to be able to return three arguments: an error-vector
% |EN|, a Jacobian |PsiN| |U2| and a projection matrix |U2|. The columns of
% this matrix |U2| must form an orthonormal basis of the subspace in which
% the cost-function does not change because of over-parametrization.
%
% This optimization implementation supports cost-functions that
% return the _R_-factor of the (projected) Jacobian |PsiN| and the
% error-vector |EN|:
%
% $$[ \Psi_N \; E_N ] = Q R$$
%
% $$[ \Psi_N U_2 \; E_N ] = Q R$$
%
% Cost-functions may use this functionality, e.g. to build up the
% _R_-factor in such a way that less memory is required. In order to
% use this feature with costfunctions that support it, the field
% |options.RFactor| should be set to |'on'|.

%% Algorithm
% This function implements a More-Hebden trust-region based
% Levenberg-Marquardt optimization according to [2,3].
% 
% In addition, this function supports projected gradients according to
% [4,5].

%% Used By
% <doptlti.html |doptlti|>, <foptlti.html |foptlti|>

%% Uses Functions
% <dfunlti.html |dfunlti|>, <ffunlti.html |ffunlti|>

%% See Also
% |lsqnonlin|, |optimset|

%% References
% [1] The MathWorks Inc., Natick, Massachusetts, _Optimization Toolbox 
% User's Guide_, version 2.1 (release 12) ed., Sept 2000.
%
% [2] J. E. Dennis and R. B. Schnabel, _Numerical Methods for Unconstrained
% Optimization and Nonlinear Equations_. New Jersey: Prentice-Hall, 1982.
%
% [3] J. J. More, "The Levenberg-Marquardt algorithm: Implemnetation and
% theory", in _Numerical Analysis (G. A. Watson, ed.), vol. 630 of _Lecture
% Notes in Mathematics_, pp. 106-116, Springer Verlag, 1978.
%
% [4] N. Bergboer, V. Verdult, and M. Verhaegen, "An effcient 
% implementation of maximum likelihood identification of LTI state-space 
% models by local gradient search", in _Proceedings of the 41st IEEE
% Conference on Decision and Control_, Las Vegas, Nevada, Dec. 2002.
%
% [5] L.H. Lee and K. Poolla, "Identification of linear parameter varying
% systems using nonlinear programming", _Journal of Dynamic Systems_,
% Measurement and Control, col. 121, pp. 71-78, Mar 1999.

