function signal = gbn(N,ts,A,h,flag)
%GBN   Generalized binary noise test-signal. Based on Tulleken.
%
% Syntax:
%           y = gbn(N,ts,A,h,flag)
% Input:
%  N        Length of the signal [s]
%  ts       Process settling time [s]
%  A        Amplitude of the signal [-]
%  h        Process sampling time [s]
%  flag     = 0 if the process is over-damped 
%           = 1 if the process is oscillary (min phase)
%           = 2 if the process is oscillary (non min phase)
%  
% Output:
%  y        Random binary noise.

% Revised by Ivo Houtzager, 2007
% Copyright (c) 2001-2007, Delft Center of Systems and Control  

% Check number of arguments
if nargin ~= 5
    error('GBN requires at five input arguments.');
end

% Help variables
start = 1;
Ns = N/h + 1;     % lenght of the signal [samples]
gbnflag = 0;      % = 1, if the generated test-signal is accepted

% Calculation of the non-switching propability
if flag == 0
    p = 1 - h/ts;
    av = 1*N/ts;    % number of average switches in the signal
elseif flag == 1
    p = 1 - 5*h/ts;
elseif flag == 2
    p = 1 - 3*h/ts;
end

% The loop is repeated until a suitable GBN signal is found
while gbnflag ~= 1
    % Generating the GBN signal
    for i=1:1:Ns;
        help = rand;
        if help > p
            ue(i) = -start*A;
            start = -start;
        else
            ue(i) = start*A;
        end
    end

    % Verifying the designed GBN-signal %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % If the process model is an overdamped system the test-signal
    % should be verified
    %
    % 1) The test signal should not be used, if the number of
    %    switches is too large (more that 120 % of the calculated
    %    average number of switches) or too less (less than 80 %
    %    of the calculated average number of switches)
    % 2) Allow for each 20 expected switches one sequence in the
    %    signal without a switch longer than 2.5*settling time
    % 3) Allow the signal to have 2 switches faster than 0.1 times
    %    the settling time for each 10 expected switches
    %
    if flag == 0
        ue = ue';
        n_s = sum(diff(ue)~=0);
        ind = find([1; (diff(ue)~=0)]);
        for i = 2 : length(ind)
            len(i-1) = ind(i) - ind(i-1);   % in samples
        end
        len = sort(len*h);                 % in seconds
        f_sw = (n_s - rem(n_s,10))/10;
        f_sw = 2*f_sw + 3;
        s_sw = round((n_s - rem(n_s,10))/20);

        if (0.8*av<n_s) && (n_s<1.2*av) && (len(f_sw)>0.1*ts) && (len(n_s-s_sw)<2.5*ts)
            signal = ue;
            gbnflag = 1;
            plot(signal)
        else
            ue = 0;
            %n_s = 0;
            %ind = 0;
            len = 0;
            %f_sw = 0;
        end
    else
        signal = ue';
        gbnflag = 1;
    end
end


