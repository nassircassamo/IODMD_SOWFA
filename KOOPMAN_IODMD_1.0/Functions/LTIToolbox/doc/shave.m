
%% SHAVE
% Reduces spikes in measured signals.

%% Syntax
% |shave(x)|
%%
% |y = shave(x)|
%%
% |y = shave(x,factor,Wn,lolim,uplim)|

%% Description
% This function is used for reducing spikes in a measured signal. The
% spikes are shaved using the method in [1].
% 
% If no output argument is specified, a figure containing the original
% signal and shaved signal is drawn. The figure also contains the band (see
% "Algorithm" below). Detected spikes are indicated with crosses.

%% Inputs
% |x| is the signal to be shaved.
% 
% |factor| is the (optional) multiplication factor which determines the
% width of the detection band. When the detection is poor, this factor
% should be changed. The default value is |2|.
% 
% |Wn| is the (optional) cut-off frequency of the low-pass filter used for
% trend determination.  It must be in the range |0.0 < Wn < 1.0|, with
% |1.0| corresponding to half the sample rate. Its default value is |0.01|.
% 
% |lolim,uplim| (optional) The signal |x| will be clipped to the band
% |[lolim,up_lim]| before the shaving starts.
         
%% Outputs
% |y| is the shaved signal.

%% Algorithm
% The spike removal algorithm developed in [1] is used. This algorithm can
% be summarized as follows:
% 
% * The trend in the signal |x| is calculated using a fourth-order
% Butterworth filter.
% * The standard deviation of the trend-corrected, clipped signal is
% calculated.
% * The detection band is defined by the signal trend plus and minus a
% certain factor times the standard deviation. All samples outside this
% band are regarded as spikes, and are replaced using linear interpolation.
%

%% Used By
% This is a top-level function that is used directly by the user.

%% References
% [1] A. Backx, _Identification of an Industrial Process: A Markov
% Parameter Approach_. PhD thesis, University of Eindhoven, Eindhoven, The
% Netherlands, 1987.

