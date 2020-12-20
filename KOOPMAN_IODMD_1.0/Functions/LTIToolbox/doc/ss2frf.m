
%% SS2FRF
% Calculates an LTI Frequency Response Function

%% Syntax
% |H = ss2frf(A,B,C,D,w)|

%% Description
% This function calculates the frequency-response function of a LTI
% state-space system that has the following model structure:
%
% $$ \dot{x}(t) = Ax(t) + Bu(t) $$
%
% $$ y(t) = Cx(t) + Du(t) $$
%
% $$ x(0) = x0 $$
%
%  or:
%
% $$  x(k+1) = Ax(k) + Bu(k)$$
%
% $$ y(k) = Cx(k) + Du(k) $$
%
%% Inputs
% |A,B,C,D| aer the system matrices describing the state-space system. If
% there is no D-matrix, pass an empty matrix.
%%
% |w| is the vector of complex frequencues at which the frequency response
% function is to be calculated. Typically, this is for continuous-time
% systems
% 
% $$w = e^{j\omega}$$
% 
% and 
% 
% $$w = j\omega $$
%
% for discrete-time systems.
    
%% Outputs
% |H| is the measured frequency response function (FRF). This is a matrix
% that follows the MATLAB 5 FRF Convention: It is a 3D-array of size _l_ x
% _m_ x _N_, in which _H(:,:,j)_ is the FRF at the _j_ th complex
% frequency.

%% Uses Functions
% <ltifrf.html |ltifrf|>

%% See Also
% <ltifrf.html |ltifrf|>



