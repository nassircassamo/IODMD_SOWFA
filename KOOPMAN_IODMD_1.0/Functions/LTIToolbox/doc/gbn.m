
%% GBN
% Produces a generalized pseudo-random binary noise test-signal.

%% Syntax
% |y = gbn(N,ts,A,h,flag)|

%% Description
% This function produces a binary sequence. This kind of testsignal has 
% been described in [1].

%% Inputs
% |N| is the lentgh of the signal [sec].
%%
% |ts| is the settling time of the process [sec].
%%
% |A| is the amplitude of the signal.
%%
% |flag| is;
%
% * 0 if the process is over-damped. 
% * 1 if the process is oscillary (minimum phase).
% * 2 if the process is oscillary (non-minimum phase).
%
%% Outputs
% |y| is a pseudo-random binary noise.

%% Used By
% This is a top-level function that is used directly by the user.

%% See Also
% <prbn.html prbn>

%% References
% [1] H. J. A. F. Tulleken, "Generalized binary noise test-signal concept
% for improved identification-experiment design", _Automatica_, vol. 26,
% no. 1, pp. 37-49, 1990.
