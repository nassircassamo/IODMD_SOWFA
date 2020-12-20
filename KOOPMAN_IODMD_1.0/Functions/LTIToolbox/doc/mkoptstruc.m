
%% MKOPTSTRUC
% Creates a MATLAB 6-compatible |optimset|-structure

%% Syntax
% |optstruc = mkoptstruc|

%% Description
% This function provides a MATLAB 6 |optimset| work-alike. It generates an
% empty |optimset|-structure that can be passed to the high-level |doptlti|
% or |foptlti| function, or to the lower-level |lmmore| function.
% 
% Note that this function only generates a default |optimset|-structure. It
% is not capable to option-merging like the MATLAB 6 |optimset| function.
% It is recommended to use |optimset| if running MATLAB 6.
         
%% Outputs
% |optstruc| is a default |optimset|-structure

%% Used By
% <optim5to6.html |optim5to6|>

%% See Also
% |optimset|

