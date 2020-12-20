%% velocitySpectra1D.m
%% FUNCTION TO COMPUTE 1D VELOCITY SPECTRA FROM TIME SERIES.

% Matthew J. Churchfield
% National Renewable Energy Laboratory (NREL)
% 15013 Denver West Parkway, Golden, CO 80501, USA
% 1.303.384.7080
% matt.churchfield@nrel.gov

% April 13, 2015

function [f,XMagSqrAvg,N] = velocitySpectra1D(t,x,nSegment,debugPrint);


% Get the time vector and associated frequency, dt, etc.
N = 2.0*floor(0.5*length(t)/nSegment);
dt = t(2)-t(1);
tmax = (N*nSegment)*dt;
dF = nSegment/tmax;
f = (dF:dF:(N/2)*dF);

% Compute spectra.
XMagSqrAvg = [];
varianceAvg = 0;
for n = 1:nSegment
    
    % apply a Hanning window to make the time series periodic
    w = window(@hann,N);
    meanSignal = mean(x((n-1)*N+1:n*N));
    
    % compute the variance
    variance = var((x((n-1)*N+1:n*N)-meanSignal));
    varianceAvg = varianceAvg + (1.0/nSegment)*variance;
    
    % call the FFT, and get the squared magnitude of it.
    X = fftshift(fft(w.*(x((n-1)*N+1:n*N)-meanSignal)));
    XMagSqr = X.*conj(X);
    
    % compute the area under the spectrum, and scale it so it is equal to 
    % the variance   
    integral = dF*trapz(XMagSqr(N/2+1:end));
    corr = variance/integral;
    XMagSqr = corr*XMagSqr;
    
    % build up an average FFT
    if (n == 1)
        XMagSqrAvg = (1.0/nSegment).*XMagSqr;
    else
        XMagSqrAvg = XMagSqrAvg + (1.0/nSegment).*XMagSqr;
    end
    
    if (debugPrint)
        disp(['Integral: ',num2str(corr*integral)]);
        disp(['Variance: ',num2str(variance)]);
        disp(['Ratio: ',num2str(corr*integral/variance)]);
    end
end