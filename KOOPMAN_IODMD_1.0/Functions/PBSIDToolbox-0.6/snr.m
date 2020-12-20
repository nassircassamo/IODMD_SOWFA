function s = snr(y,y0)
%SNR Signal to noise ratio [dB]
%  S = SNR(Y,Y0) calculates the Signal to Noise Ratio (SNR) of the signals
%  Y and Y0, where Y is the noise corrupted output of the system and Y0 the
%  output of the deterministic transfer only.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010


% check input arguments
if nargin ~= 2
    error('SNR requires two input arguments!')
end

% check dimensions of inputs
if size(y,2) > size(y,1)
    y = y';
end
if size(y0,2) > size(y0,1)
    y0 = y0';
end
N = size(y,1);
if size(y0,1) ~= N
    error('Both signals should have an equal number of samples.');
end
if size(y,2) ~= size(y0,2)
    error('Both signals should have an equal number of components.');
end

% calculate signal to noise ratio
s = 10*log10(var(y)./var(y-y0));