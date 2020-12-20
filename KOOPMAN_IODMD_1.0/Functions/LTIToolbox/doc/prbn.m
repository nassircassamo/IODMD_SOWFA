
%% PRBN
% Produces a pseudo-random binary sequency suitable as identification
% input.

%% Syntax
% |y = prbn(N,rate)|

%% Description
% This function produces a binary sequence, with values _0_ and _1_. The
% chance of switching from level is given by the parameter |rate|. |rate =
% 0| yields a constant value _0_. |rate = 1| gives an signal that changes
% between _0_ and _1_ at every time-instant. Any value in between results
% in a random binary sequence. This kind of testsignal has been described
% in [1].

%% Inputs
% |N| are the number of samples.
%%
% |rate| is the (optional) probability of the signal changing level at each
% time-instant. The default value is _0.5_.
          
%% Outputs
% |y| is a pseudo-random binary noise.

%% Used By
% This is a top-level function that is used directly by the user.

%% References
% [1] H. J. A. F. Tulleken, "Generalized binary noise test-signal concept
% for improved identification-experiment design", _Automatica_, vol. 26,
% no. 1, pp. 37-49, 1990.
