function [Af,Ab,Sf,Sb] = destmar(v,d) 
%DESTMAR    This function fits a multivariable autoregressive model 
%           to a time-series v(k). The model-structure is: 
% 
%           v(k) = ef(k) - Af_1 v(k-1) - ... - Af_d v(k-d) 
%           v(k) = eb(k) - Ab_1 v(k+1) - ... - Ab_d v(k+d) 
% 
%           in which ef(k) and eb(k) are innovation sequences with 
%           covariance matrices Sf and Sb respectively. 
% 
% Syntax: 
%           [Af,Ab,Sf,Sb] = destmar(v,d); 
% 
% Inputs: 
% v         The time-series, a N x l matrix for a signal having N samples 
%           and which is l-dimensional. 
% d         The order d of the AR model. 
% 
% Outputs: 
% Af,Ab     The coefficient matrices of the causal and anticausal model. 
% Sf,Sb     The covariance matrices of the causal and anticausal 
%           innovations. 
% 
% See Also : cholicm 
 
% Written by Niek Bergboer, December 2001 
% Revised by Niek Bergboer, 2002
% Revised by Ivo Houtzager, 2007
% Copyright (c) 2001-2007, Delft Center of Systems and Control 

% Check number of arguments
if nargin < 2
    error('DESTMAR requires two input arguments.');
end

% Extract data
N = size(v,1);
l = size(v,2);

% Check whether there is at least one signal
if l < 1
    error('DESTMAR requires a signal');
end

% Check whether enough samples are available.
if N < (l+1)*d,
    error('Data sequence is too short; must be at least (l+1)*d samples');
end

% Forward problem (causal)
% Fill GammaF
GammaF = zeros(d*l,N-d);
for i = 1:l
    GammaF((d-1)*l+i:-l:i,:) = hankel(v(1:d,i),v(d:N-1,i)');
end;
gammaF = zeros(l,N-d);

% Fill gammaF
for i = 1:l
    gammaF(i,:) = -v(d+1:N,i)';
end

% Solve regression
Rf = triu(qr([GammaF' gammaF']));
clear GammaF;
clear gammaF;
Rf = Rf(1:(d+1)*l,1:(d+1)*l);
Af = (Rf(1:d*l,1:d*l)\ Rf(1:d*l,d*l+1:(d+1)*l))';

% Estimate covariance matrix
Sf = (Rf(d*l+1:(d+1)*l,d*l+1:(d+1)*l)'*Rf(d*l+1:(d+1)*l,d*l+1:(d+1)*l))/(N-d);

% Backward problem (anti-causal)
% Fill GammaB
GammaB = zeros(d*l,N-d);
for i = 1:l
    GammaB(i:l:(d-1)*l+i,:) = hankel(v(2:d+1,i),v(d+1:N,i)');
end;

% Fill gammaB
gammaB = zeros(l,N-d);
for i = 1:l
    gammaB(i,:) = -v(1:N-d,i)';
end;

% Solve regression
Rb = triu(qr([GammaB' gammaB']));
clear GammaB; clear gammaB;
Rb = Rb(1:(d+1)*l,1:(d+1)*l);
Ab = (Rb(1:d*l,1:d*l)\ Rb(1:d*l,d*l+1:(d+1)*l))';

% Estimate Covariance Matrix
Sb = (Rb(d*l+1:(d+1)*l,d*l+1:(d+1)*l)'*Rb(d*l+1:(d+1)*l,d*l+1:(d+1)*l))/(N-d);








