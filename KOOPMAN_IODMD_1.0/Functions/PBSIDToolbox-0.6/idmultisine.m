function [s,t] = idmultisine(N,ampl,ts,nu,F,Fstop,t0,p)
%IDMULTISINE Generation of a Multi Sine Signal
%  [S,T] = IDMULTISINE(N,AMPL,Ts,nu,F,Fstop,T0,P) generates a Multi Sine
%  Signal with the following properties:
% 
%  N                : signal length [s]
%  AMPL (optional)  : amplitude (default: 1)
%  Ts (optional)    : sampling time [s] (default: 1)
%  nu (optional)    : discrete-frequency spacing [Hz] (default: 0.1)
%  F (optional)     : cutoff frequency [Hz] (default: 0.5/pars.ts)
%  Fstop (optional) : PRBS will be band-stop filtered around this frequency 
%                     (default: inf)
%  T0 (optional)    : starting time [s] (default: 0)
%  P (optional)     : number of channels (default: 1)
%
%  See also: IDINPUTS, IDPRBS.

% check number if input arguments
if nargin < 1
    error('IDMULTISINE requires at least one input argument.')
end
if nargin < 8 || isempty(p)
    p = 1;
end
if nargin < 7 || isempty(t0)
    t0 = 0;
end
if nargin < 6 || isempty(Fstop)
    Fstop = inf;
end
if nargin < 4 || isempty(nu)
    nu = 0.1;
end
if nargin < 3 || isempty(ts)
    ts = 1;
end
if nargin < 5 || isempty(F)
    F = 0.5/ts;
end
if nargin < 2 || isempty(ampl)
    ampl = 1;
end

if F==0
    warning('ID:cutoffzero','Selected cutoff frequency is zero, returning zero PRBS signal...')
    s = zeros(1,ceil(N/ts));
    t = t0 + ts*(0:1:length(s)-1)';
    return
end;

Fnyq = 0.5/ts;     % Nyquist frequency
if F > Fnyq
    warning('ID:bandbeyondnyqst','Specified bandwidth [0, %s] Hz is beyond the Nyquist frequency %s Hz',num2str(F),num2str(Fnyq))
    F = Fnyq;
    warning('ID:bandbeyondnyqst','>>> BANDWIDTH MODIFIED TO [0, %s] Hz <<<',num2str(F))
end;

% Determine frequencies
Freq = nu:nu:Fnyq;

% band-stop filtering multisine signal
if Fstop < F
    Il = (Freq < 0.7*Fstop);
    Iu = (Freq > 1.3*Fstop);
    Freq = [Freq(Il) Freq(Iu)];
else
    warning('ID:bandlargercutoff','The frequency to be filtered out, %s Hz, is larger than the cutoff frequency %s Hz',num2str(Fstop), num2str(F))
    warning('ID:bandlargercutoff','>>> Skipping band stop filtering around %s Hz <<<',num2str(Fstop))
end
d = size(Freq,2);

% time stamps for s
t = t0 + ts*(0:1:N-1)';

% Create multisine
%a = 1;
s = zeros(length(t),p);
for k = 1:d
    for m = 1:p
        a = 1/Freq(k);
        s(:,m) = s(:,m) + a*sin((2*pi*(Freq(k)+(m-1)*(nu/p))).*t - ((k*(k-1)*pi)/d).*ones(length(t),1));
    end
end
for m = 1:p
    s(:,m) = (ampl/max(s(:,m))).*s(:,m);
end

end

