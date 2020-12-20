
%% DESTMAR
% Fits a multivariable autoregressive model to a time-series.

%% Syntax
% |[Af,Ab,Sf,Sb] = destmar(v,d)|

%% Description
% This function fits a multivariable autoregressive model to a time-series
% _v(k)_. The model-structure is
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
% in which _ef(k)_ and _eb(k)_ are innovation sequences with covariance
% matrices _Sf_ and _Sb_ respectively. The fitting is performed according
% to [1].

%% Inputs
% |v| is the time-series, a _N_ x _l_ matrix for a signal having _N_
% samples and which is _l_-dimensional.
%%
% |d| is the desired order _d_ of the AR model.

%% Outputs
% |Af,Ab| are the coefficient matrices |Af| and |Ab| of the causal and
% anticausal model.
%%
% |Sf,Sb| are the covariance matrices |Sf| and |Sb| of the causal and
% anticausal innovations.

%% Algorithm
% A direct Hankel-matrix based estimation of the AR model is performed
% according to [1].

%% Used By
% This is a top-level function that is used directly by the user.

%% See Also
% <cholicm.m |cholicm|>, <doptlti.m |doptlti|>

%% References
% [1] B. Davis, _Parameter Estimation in Nonlinear Dynamical Systems with
% Correlated Noise._ PhD thesis, Universite Catholique de Louvain-La-Neuve,
% Belgium, 2001.