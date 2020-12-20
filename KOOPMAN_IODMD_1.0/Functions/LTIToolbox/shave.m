function Shaved_output = shave(Signal,factor,Wn,lo_lim,up_lim)
%SHAVE     This function is used for reducing spikes from a
%          measured signal.
%          The spikes of the signal are 'shaved' as follows:
%          -  From the signal a trend is computed using a
%             fourth-order Butterworth filter.
%          -  The standard deviation of the trend-corrected, clipped
%             signal is computed.
%          -  A detection band is defined by the trend plus and minus a
%             certain factor times the standard deviation. All samples
%             that are outside this band are replaced using linear
%             interpolation.
%          This 'shaving' method has been described in:
%          A. Backx, 'Identification of an Industrial Process: A Markov
%          Parameter Approach', PhD thesis, University of Eindhoven,
%          The Netherlands, November 1987.
% 
%          If no output argument is specified, a figure is drawn which
%          contains the original signal and the `shaved' signal.
%          The band is also plotted.
%          The spikes found are indicated with crosses.
% 
% Syntax:
%           shave(Signal)
%           Shaved_signal=shave(Signal)
%           Shaved_signal=shave(Signal,factor,Wn,lo_lim,up_lim)
% 
% Inputs:
%   Signal         Signal to be 'shaved' (column vector);
%   factor         Multiplication factor which determines the
%                  width of the detection band. When the
%                  detection is poor, you should change this
%                  factor. This argument is optional.
%                  Default is 2.
%   Wn             Cut-off frequency of the low-pass filter used
%                  for trend determination.  It must be in the
%                  range 0.0 < Wn < 1.0, with 1.0 corresponding
%                  to half the sample rate. This argument is
%                  optional. Its default value is 0.01.
%   lo_lim,        If these arguments are present, the signal is
%   up_lim         clipped to a minimum value of lo_lim and a
%                  maximum value of up_lim before the
%                  'shaving' starts.
% 
% Outputs:
%   Shaved_signal  'Shaved' signal.

% Vincent Verdult, August 1997
% Revised by Ivo Houtzager, 2007
% Copyright (c) 1997-2007, Delft Center of Systems and Control 

if nargin == 4
    error('You need to specify two limits.')
end
if nargin < 3
    Wn = 0.01;
end
if nargin < 2
    factor = 2;
end
if factor <= 0
    error('Factor should be larger than zero.')
end
if Wn <= 0 || Wn > 1
    error('Cut-off frequency Wn must be 0 < Wn < 1.')
end
N_samples = length(Signal);
if size(Signal,2) > size(Signal,1)
    Signal = Signal';
end
if size(Signal,2) ~= 1
    error('This function only works for vector input.')
end

Clip_signal = Signal;
if nargin == 5
    if lo_lim >= up_lim
        error('Upper clipping limit should be larger than lower limit')
    end
    for i = 1:N_samples
        if Signal(i) > up_lim
            Clip_signal(i) = up_lim;
        elseif Signal(i)<lo_lim
            Clip_signal(i) = lo_lim;
        end
    end
end

[B,A] = butter(4,Wn);
mean_Signal = mean(Clip_signal);
Clip_signal = [mean_Signal; Clip_signal; mean_Signal;];
Trend = filtfilt(B,A,Clip_signal);
Clip_signal = Clip_signal(2:N_samples+1);
Trend = Trend(2:N_samples+1);
Cor_signal = Clip_signal - Trend;
std_Signal = std(Cor_signal);
Up_bound = Trend+factor*std_Signal;
Lo_bound = Trend-factor*std_Signal;


bad_samples = zeros(1,N_samples);
for i = 1:N_samples
    bad_samples(i) = (Signal(i)>Up_bound(i)|Signal(i)<Lo_bound(i));
end
dif_bad_samples = abs(diff(bad_samples));
times = find(dif_bad_samples>0);


Shaved_signal = Signal;
if bad_samples(1)==1
    len_times = length(times);
    Shaved_signal(1:times(1)) = ones(times(1),1).*Signal(times(1)+1);
    times = times(2:len_times);
end
if bad_samples(N_samples)==1
    len_times = length(times);
    Shaved_signal(times(len_times)+1:N_samples) = ones(N_samples-times(len_times),1).*Signal(times(len_times));
    times = times(1:len_times-1);
end
len_times = length(times);
for i = 1:2:len_times;
    t1 = times(i);
    t2 = times(i+1)+1;
    a = (Signal(t2)-Signal(t1))/(t2-t1);
    b = Signal(t1)-a * t1;
    time_vec = t1+1:t2-1;
    Shaved_signal(time_vec) = a * time_vec+b;
end

if nargout==0
    subplot(2,1,1)
    plot(Signal)
    ylabel('Original Signal')
    hold on
    plot(Up_bound,'m')
    plot(Lo_bound,'m')
    ax3 = min(Lo_bound);
    ax4 = max(Up_bound);
    axis([0,N_samples,ax3,ax4])
    for i = 1:N_samples
        if bad_samples(i)==1
            plot(i,bad_samples(i) * ax4,'gx')
        end
    end
    hold off
    subplot(2,1,2)
    plot(Shaved_signal)
    ylabel('Shaved signal')
    return
end
Shaved_output = Shaved_signal;
