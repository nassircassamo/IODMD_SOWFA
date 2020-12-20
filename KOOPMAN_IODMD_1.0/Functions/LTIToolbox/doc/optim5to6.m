
%% OPTIM5TO6
% Translates a |foptions|-vector into an |optimset|-structure.

%% Syntax
% |options = optim5to6(fopts)|

%% Description
% This function translates a MATLAB 5 optimization options vector ---as
% generated using |foptions| --- into a MATLAB 6 compatible
% |optimset|-structure. Translated fields are:
% 
% * *1* Display
% * *2* TolX
% * *3* TolFun
% * *9* Jacobian
% * *14* MaxFunEval
%

%% Inputs
% |fopts| is a MATLAB 5 compatible |foptions|-vector
         
%% Outputs
% |options| is a MATLAB 6 compatible |optimset|-structure

%% Remarks
% MATLAB 5 uses a default parameter and function tolerance of |1e-4|. This
% is indicated by the second and third element of |fopts|, that are |1e-4|
% in the default case.
% 
% MATLAB 6 uses a default value of |1e-6| for both tolerances, but setting
% the tolerances to |1e-6| if the |fopts| vector contains the default
% values is impossible: there is no way of telling whether the user used
% the default values or that he actually specified |1e-4| as tolerance.
% 
% Consequently, the tolerances are copied verbatim, _and there will thus be
% different results in an optimization when using a default
% |foptions|-vector or a default |optimset|-structure_.

%% Used By
% <mkoptstruc.html |mkoptstruc|>

%% See Also
% |foptions|, |optimset|, <mkoptstruc.html |mkoptstruc|>

