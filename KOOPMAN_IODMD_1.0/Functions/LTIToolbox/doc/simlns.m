
%% SIMLNS
% Calculates the left null-space of the basis of similarity transformations.

%% Syntax
% |U2 = simlns(A,B,C,[],[],[])|
%%
% |U2 = simlns(A,B,C,K,fD,fx)|

%% Description
% The function <simlns.html |simlns|> calculates the left null-space of an
% LTI system's similarity map |Mtheta|. In the most general case, when _A_,
% _B_, _C_, _D_, _K_ and _x0_ are part of the parameter vector, this matrix
% is given by [1]:
%%
% <<simlns.jpg>>
%%
% A QR-factorization is used to obtain the left null-space.
%            
% This function is used internally by <dfunlti.html |dfunlti|> and
% <ffunlti.html |ffunlti|> and is not meant for stand-alone use.

%% Inputs
% |A,B,C| are the system matrices describing the LTI State Space system.
%%
% |K| is the (optional) Kalman gain, specify as empty matrix when not
% present.
%%  
% |fD| (optional) specifies whether _D_ is part of the parameter vector,
% specify as empty, _0_ or _1_.
%%                
% |fx| (optional) specifies whether _x0_ is part of the parameter vector,
% specify as empty, _0_ or _1_.
         
%% Outputs
% |U2| is the left null-space of the similarity map.

%% Remarks
% Specifying |fx=1| only causes an _n_ x _n_ identify-matrix to be appended
% to the lower right of the left null-space matrix; in a non-linear
% optimization, applying the left null-space ensures that the state-basis
% does not change. It thus does not have to be projected.

%% Algorithm
% The manifold matrix |Mtheta| is calculated according to [1]. A
% QR-factorization is used subsequently to obtain the left null-space

%% Used By
% <dfunlti.html |dfunlti|>, <ffunlti.html |ffunlti|>

%% References
% [1] L. H. Lee and K. Poolla, "Identification of linear parameter-varying
% systems using nonlinear programming", _Journal of Dynamic Systems,
% Measurement and Control_, vol. 121, pp. 71-78, Mar. 1999.

