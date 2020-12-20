%% Release Notes

%% Predictor-Based Subspace IDentification Toolbox
%
% V0.6 (February 2015)
% 
% * Fixed a number of bugs
% V0.5 (January 2012)
% 
% * added BPDN regularization for |dordvarx| and |lordvarx|, see example 19b
% * added predictor stabilization option for |lx2abcdk|, see example 19c
% * function |predstab| to assess stability of the predictor form
% * more memory-efficient batch-wise computations in |dordvarx|
% * added function |bodemag| for |idafflpv| object, for LPV Bode magnitude plots
% * correction of bugs in |lordvarx| and |regress|
% * option to use obs. matrix with mu= [1,1..,1] in |lordvarx|
%
% V0.4 (December 2010)
%
% * added functions to estimate the asymptoic variance of the PBSIDopt
% (VARX) estimates
% * added functions for plotting the bode and eigenvalues with error bounds
%
% V0.3 (Oktober 2010)
%
% * regularisation parameter can be now be given by user
% * added the recursive function |rpbsid|
% * added the simulink blocks for real-time recursive identification
% * added examples with the recursive function |rpbsid|
%
% V0.2 (Augustus 2010)
%
% * improved and new functions (|pem|, |pe|, etc) for IDAFFLPV object
% * new functions to convert to IDSS object
% * batch updating of |dordfir|, |dordvarx|, and |dordvarmax| improved
% * new function |idmultisine| to generate multi sine excitation signals
% * add stability option in |dx2abcdk.m|
% * corrected bugs in |lordvarx| and |lx2abcdk|
% * added more documentation
%
% V0.1 (December 2009)
%
% * first release of the PBSID Toolbox
