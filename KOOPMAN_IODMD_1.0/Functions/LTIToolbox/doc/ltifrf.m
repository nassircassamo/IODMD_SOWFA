
%% LTIFRF
% Calculates an LTI Frequency Response Function

%% Syntax
% |H = ltifrf(A,B,C,[],[],w,outopt)|
%%
% |H = ltifrf(A,B,C,D,[],w,outopt)|
%%
% |H = ltifrf([],[],[],D,[],w,outopt)|
%%
% |H = ltifrf(A,B,C,[],dA,w,outopt)|

%% Description
% <ltifrf.html |ltifrf|> will return the Frequency Response Function (FRF)
% of a linear time-invariant state-space model, evaluated at the complex
% frequencies provided in _w_:
% 
% $$H = C {(\mathtt{w} I_n - A)}^{-1} B + D$$
% 
% This function is used internally by <ffunlti.html |ffunlti|>, <fac2b.html
% |fac2b|> and <fac2bd.html |fac2bd|>. It is not meant for stand-alone use.

%% Inputs
% |A| is the state-space model matrix _A_.
%% 
% |B|	is the state-space model matrix _B_.
%% 
% |C|	is the state-space model matrix _C_.
%% 
% |D|	is the (optional) state-space model matrix _D_.
%% 
% |dA| (optional) calculates the change in FRF given the deviation _dA_ in
% _A_. _D_ and _dA_ are mutually exclusive.
%% 
% |w| is the vector of complex frequencies. For discrete-time systems:
% 
% $$e^{j\omega}$$
% 
% and for continuous-time systems.
% 
% $$ j\omega $$
%% 
% |outopt| controls how _H_ will be returned (see below).
    
%% Outputs
% |H| is the FRF. Usually a 3D-array of size _l_ x _m_ x _N_. However, if
% |outopt| is non-empty and _1_, _H_ will be a vector of size _lmN_ x 1. If
% |outopt| is non-empty and _2_, _H_ will be a matrix of size _l_ x _mN_.

%% Algorithm
% The state-space model is first transformed such that its
% state-transistion matrix _A_ is in upper-Hessenberg form. The matrix
% 
% $${(\mathtt{w} I_n - A)}^{-1} B$$
% 
% is subsequently solved by an efficient upper-Hessenberg solver in SLICOT,
% after which premultiplication by _C_ and addition of _D_ yields the FRF.
% This approach follows [1].
% 
% If a deviation _deltaA_ in _A_ is given, the FRF deviation is given by:
% 
% $$ \delta H = C {(\mathtt{w} I_n - A)}^{-1} \delta A {(\mathtt{w} I_n -
% A)}^{-1} B $$
% 
% Again, the model is transformed so that _A_ has upper-Hessenberg form,
% after which the SLICOT Hessenberg solver is used to obtain
% 
% $${(\omega I_n - A)}^{-1} B$$
% 
% and
% 
% $${(\omega I_n - A)}^{-1} \delta A$$
% 
% Multiplication then yeilds the FRF deviation.

%% Used By
% <ffunlti.html |ffunlti|>, <fac2b.html |dac2b|>, <fac2bd.html |dac2bd|>

%% Uses Functions
% SLICOT-functions |MB02RZ|, |MB02SZ|, |TB05AD|
%%
% LAPACK-functions |DGEHRD| and |DORMHR|
%%
% (All built into the executable)

%% See Also
% <ffunlti.html |ffunlti|>, <fac2b.html |fac2b|>, <fac2bd.html |fac2bd|>

%% References
% [1] A.J. Laub, "Efficient multivariable frequency response calculations",
% _IEEE Transactions on Automatic Control_, vol. 26, pp. 407-408, Apr.
% 1981.


