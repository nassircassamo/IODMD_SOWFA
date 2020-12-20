
%% Release Notes
%
%% LTI System Identification Toolbox
%
% V2.2 (July 2010)
%
% * automated mex file generation using |ltisetup|
% * add option to estimate stable A matrix in |dmodpo| and |dmodpi|
% * add |dmodpi4| and |dmodpo4| where A and C are found via MOESP, B and D via N4SID
%
% V2.1 (July 2009)
%
% * add the work of C.T. Chou, R. Hallouzi, and K. Hinnen
% * new mex-binaries for windows (compiled on Vista/Matlab2008b)
%
% V2.02 (April 2008)
%
% * fixed handling of empty inputs in |dordpo| and |dmodpo|
%
% V2.01 (January 2008)
%
% * fixed bug in |fac2bd|
% * fixed bug calling |dac2bdc|
% * new mex-binaries for windows (compiled on XP/Matlab2007b)
%
% V2.0 (August 2007)
%
% * new help interface
% * removed depreciated SMI v1.0 functions (see below)
% * removed LPV functions (becomes a seperate LPV Toolbox)
% * changed |lpvitr| to |ltiitr| (removed LPV functionality)
% * renamed |example| to |ltidemo|
% * new function |gbn| (based on |prbn|)
% * add check arguments to MEX-files and M-files
% * corrected a lot of M-Lint messages
% * new mex-binaries for windows (compiled on XP/Matlab2006b)
% * new mex-binaries for linux (compiled on FEDORACORE6/Matlab2006b)
%
% Previous versions of the LTI System Identification Toolbox are v1.5, v1.0
% and v0.9 (April 2002). Release notes of these previous versions are not available. 
%
%
%% State space Model Identification (SMI) Toolbox.
% The toolbox software is based partly on the SMI Toolbox v1.0e 
% (August 1999), and forms a superset of most of its functionality, see the
% companion (software manual).