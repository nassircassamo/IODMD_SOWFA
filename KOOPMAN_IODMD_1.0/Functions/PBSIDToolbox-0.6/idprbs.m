function [s,t] = idprbs(N,ampl,ts,F,Fstop,t0,p)
%IDPRBS Generation of a Pseudo-Random Binary Signal (PRBS)
%  [S,T] = IDPRBS(N,AMPL,Ts,F,Fstop,T0,P) generatos a Pseudo-Random Binary
%  Signal (PRBS) with the following properties:
% 
%  N                : signal length [s]
%  AMPL (optional)  : amplitude (default: 1)
%  Ts (optional)    : sampling time [s] (default: 1)
%  F (optional)     : cutoff frequency [Hz] (default: 0.5/pars.ts)
%  Fstop (optional) : PRBS will be band-stop filtered around this frequency 
%                     (default: inf)
%  T0 (optional)    : starting time [s] (default: 0)
%  P (optional)     : number of channels (default: 1)
%
%  See also IDINPUTS, IDMULTISINE.

% check number if input arguments
if nargin < 1
    error('ID requires at least one input argument.')
end
if nargin < 7 || isempty(p)
    p = 1;
end
if nargin < 6 || isempty(t0)
    t0 = 0;
end
if nargin < 5 || isempty(Fstop)
    Fstop = inf;
end
if nargin < 3 || isempty(ts)
    ts = 1;
end
if nargin < 4 || isempty(F)
    F = 0.5/ts;
end
if nargin < 2 || isempty(ampl)
    ampl = 1;
end

if F==0
    warning('ID:cutoffzero','Selected cutoff frequency is zero, returning zero PRBS signal...')
    s = zeros(1,ceil(N/ts));
    t = t0 + ts * (0:1:length(s)-1)';
    return
end;

Fnyq = 0.5/ts;     % Nyquist frequency
if F > Fnyq
    warning('ID:bandbeyondnyqst','Specified bandwidth [0, %s] Hz is beyond the Nyquist frequency %s Hz',num2str(F),num2str(Fnyq))
    F = Fnyq;
    warning('ID:bandbeyondnyqst','>>> BANDWIDTH MODIFIED TO [0, %s] Hz <<<',num2str(F))
end;

% Making L=2^N-1
N = ceil(log2(N/ts + 1));
L=2^N-1;

% Bandwidth as fraction of Nyquist frequency
FbyFnyq = F/Fnyq;
Bnd = [0, FbyFnyq];
s = idinput([L, p],'prbs',Bnd,ampl*[-1,1]);

% Low-pass filtering PRBS signal
if FbyFnyq < 1
    Fircoefs = fir1(2000, FbyFnyq, kaiser(2000+1,300));
    %Fircoefs = fir1(100,FbyFnyq,kaiser(100+1,5));
    s = filter(Fircoefs,1,s);
end

% band-stop filtering PRBS signal
if Fstop < F
    Bnd = [0.7, 1.3] * Fstop / Fnyq;
    Bnd(2) = min(1, Bnd(2)); % make sure Fnyq > Bnd*Fnyq
    [B,A] = ellip(4,1,20,Bnd,'stop');
    s = filter(B,A,s);
else
    warning('ID:bandlargercutoff','The frequency to be filtered out, %s Hz, is larger than the cutoff frequency %s Hz',num2str(Fstop), num2str(F))
    warning('ID:bandlargercutoff','>>> Skipping band stop filtering around %s Hz <<<',num2str(Fstop))
end

% time stamps for s
t = t0 + ts * (0:1:length(s)-1)';
