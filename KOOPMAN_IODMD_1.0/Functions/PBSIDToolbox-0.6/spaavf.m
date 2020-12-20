function [G,w,Coh] = spaavf(u,y,r,dt,Nband,Nfft,ZeroPadding,Wname)
%SPAAVF Spectral analysis with frequency averaging
%   [G,W]=SPAAVF(U,Y,Ts,Nband) determines a frequency-domain estimate
%   SYS=FRD(G,W) of the transfer function of the plant. The sample time is
%   given in Ts. Nband is the number of frequency bands to average.
%   Averaging in the frequency domain is used to get a smoother frequency
%   response function. The spectrum is smoothed locally in the region of
%   the target frequencies, as a weighted average of values to the right
%   and left of a target frequency. The variance of the spectrum will
%   decrease as the number of frequencies used in the smoothing increases.
%   As the bandwidth increases, more spectral ordinates are averaged, and
%   hence the resulting estimator becomes smoother, more stable and has
%   smaller variance.
%
%   [G,W]=SPAAVF(U,Y,R,Ts,Nband) determines a frequency-domain estimate
%   SYS=FRD(G,W) of the transfer function of the plant operating in
%   closed-loop. Because the conventional transfer function estimate, will
%   give a biased estimate under closed-loop [2]. An unbiased alternative
%   is to use cross-spectral between the input/output signals with an
%   external excitation signal r [1]. Hence, we define the estimate: 
%             G(exp(j*omega)) = Phi_yr(omega)*inv(Phi_ur(omega))
%
%   [G,W]=SPAAVF(...,Ts,Nband,Nfft) specifies the number of evaluated
%   frequencies. For large data sequences it is wortwhile to choose the
%   value Nfft as function of the power of two. In this case a faster
%   method is used durring the FFT. Nfft <= length(u), unless zeros are
%   added. See also, FFT. 
%
%   [G,W]=SPAAVF(...,Ts,Nband,Nfft,ZeroPading) adds additional zeros to
%   the data sequences. Usefull for increasing the number of evaluated
%   frequencies.
%
%   [G,W]=SPAAVF(...,Ts,Nband,Nfft,ZeroPading,Wname) specifies and aplies
%   an window to the data. Windowing weigths the data, it increases the
%   importance of the data in the middle of the vector and decreases the
%   importance of the data at the end and the beginning, thus reducing the
%   effect of spectral leakage. See also, WINDOW.
%
% References:
%   [1] Akaike, H., Some problems in the application of the cross-spectral 
%       method, In spectral analysis of time series, pp. 81-107, Wiley, 
%       New York, 1967.
%   [2] van den Hof, P., System Identification, Lecture Notes, Delft, 2007.

%  Revision 2: Now also works properly for MIMO cases.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% Check closed-loop
if isequal(length(u),length(r));
    clmode = 1;
else
    clmode = 0;
    if nargin == 7
        Wname = ZeroPadding;
    end
    if nargin == 6
        ZeroPadding = Nfft;
    end
    if nargin == 5
        Nfft = Nband;
    end
    Nband = dt;
    dt = r;
end

% Transpose vectors if needed
if size(u,2) < size(u,1);
    u = u';
end
nu = size(u,1); % number of inputs
if size(y,2) < size(y,1);
    y = y';
end
ny = size(y,1); % number of inputs
if clmode
    if size(r,2) < size(r,1);
        r = r';
    end
    nr = size(r,1); % number of inputs
end

% Apply window and zeros if needed
if nargin == (7+clmode)
   u = u.*(ones(nu,1)*window(Wname,size(u,2))');
   y = y.*(ones(ny,1)*window(Wname,size(u,2))');
   if clmode
       r = r.*(ones(nr,1)*window(Wname,length(u))');
   end
end
if nargin == (6+clmode)
   u = [u zeros(nu,ZeroPadding)]; 
   y = [y zeros(ny,ZeroPadding)];
   if clmode
       r = [r zeros(nr,ZeroPadding)];
   end
end
if nargin < (5+clmode)
    Nfft = [];
end

% Some administration
Fs = 1./dt;                         % Sample frequency
N  = length(u);                     % Number of samples
T  = N*dt;                          % Total time
if isempty(Nfft) 
    f = (0:N-1)'/T;                 % Frequency vector (double sided)
else
    f = (linspace(0,N-1,Nfft))'/T;  % Frequency vector (double sided)
end
Nf = length(f);

% Determine Fourier transforms
U = zeros(nu,1,Nf);
Y = zeros(ny,1,Nf);
for i = 1:nu
    ut = dt*fft(u(i,:),Nfft);
    for k = 1:Nf
        U(i,:,k) = ut(k);
    end
end
for i = 1:ny
    yt = dt*fft(y(i,:),Nfft);
    for k = 1:Nf
        Y(i,:,k) = yt(k);
    end
end
if clmode
    R = zeros(nr,1,Nf);
    for i = 1:nr
        rt = dt*fft(r(i,:),Nfft);
        for k = 1:Nf
            R(i,:,k) = rt(k);
        end
    end
end

% Do closed-loop or open-loop
if clmode
    % Determine spectral densities
    Sur = zeros(nu,nr,Nf);
    Srr = zeros(nr,nr,Nf);
    Syr = zeros(ny,nr,Nf);
    Syy = zeros(ny,ny,Nf);
    for k = 1:Nf
        Sur(:,:,k) = (1/T).*U(:,:,k)*R(:,:,k)';
        Srr(:,:,k) = (1/T).*R(:,:,k)*R(:,:,k)';
        Syr(:,:,k) = (1/T).*Y(:,:,k)*R(:,:,k)';
        Syy(:,:,k) = (1/T).*Y(:,:,k)*Y(:,:,k)';
    end

    % Apply frequency averaging
    Nmod = floor(Nf/Nband);
    mSur = zeros(nu,nr,Nmod);
    mSrr = zeros(nr,nr,Nmod);
    mSyr = zeros(ny,nr,Nmod);
    mSyy = zeros(ny,ny,Nmod);
    mf = freqAvg(f,Nband);
    for i = 1:nr
        for j = 1:nr
            tSur = Sur(i,j,:);
            tSur = freqAvg(tSur(:),Nband);
            mSur(i,j,:) = tSur;
            tSrr = Srr(i,j,:);
            tSrr = freqAvg(tSrr(:),Nband);
            mSrr(i,j,:) = tSrr;
        end
    end
    for i = 1:ny
        for j = 1:ny
            tSyy = Syy(i,j,:);
            tSyy = freqAvg(tSyy(:),Nband);
            mSyy(i,j,:) = tSyy;
        end
        for j = 1:nr
            tSyr = Syr(i,j,:);
            tSyr = freqAvg(tSyr(:),Nband);
            mSyr(i,j,:) = tSyr;
        end
    end

    % Estimate transfer function
    G = zeros(ny,nr,Nmod);
    for k = 1:Nmod
        G(:,:,k) = mSyr(:,:,k)/mSur(:,:,k);
    end

    % Squared coherence Cohuy function between r and y
    if nr == ny && nargout == 3
        Coh = zeros(ny,nr,Nmod);
        for k = 1:Nmod
            Coh(:,:,k) = sqrtm(abs(mSyr(:,:,k))*abs(mSyr(:,:,k))/(mSyy(:,:,k)*mSrr(:,:,k)));
        end
    end
else
    % Determine spectral densities
    Suu = zeros(nu,nu,Nf);
    Syu = zeros(ny,nu,Nf);
    Syy = zeros(ny,ny,Nf);
    for k = 1:Nf
        Suu(:,:,k) = (1/T).*U(:,:,k)*U(:,:,k)';
        Syu(:,:,k) = (1/T).*Y(:,:,k)*U(:,:,k)';
        Syy(:,:,k) = (1/T).*Y(:,:,k)*Y(:,:,k)';
    end
    
    % Apply frequency averaging
    Nmod = floor(Nf/Nband);
    mSuu = zeros(nu,nu,Nmod);
    mSyu = zeros(ny,nu,Nmod);
    mSyy = zeros(ny,ny,Nmod);
    mf = freqAvg(f,Nband);
    for i = 1:nu
        for j = 1:nu
            tSuu = Suu(i,j,:);
            tSuu = freqAvg(tSuu(:),Nband);
            mSuu(i,j,:) = tSuu;
        end
    end
    for i = 1:ny
        for j = 1:ny
            tSyy = Syy(i,j,:);
            tSyy = freqAvg(tSyy(:),Nband);
            mSyy(i,j,:) = tSyy;
        end
        for j = 1:nu
            tSyu = Syu(i,j,:);
            tSyu = freqAvg(tSyu(:),Nband);
            mSyu(i,j,:) = tSyu;
        end
    end

    % Estimate transfer function
    G = zeros(ny,nu,Nmod);
    for k = 1:Nmod
        G(:,:,k) = conj(mSyu(:,:,k)*pinv(mSuu(:,:,k)));
    end

    % Squared coherence Cohuy function between u and y
    if nu == ny && nargout == 3
        Coh = zeros(ny,nu,Nmod);
        for k = 1:Nmod
            Coh(:,:,k) = conj(sqrtm(abs(mSyu(:,:,k))*abs(mSyu(:,:,k))*pinv(mSyy(:,:,k)*mSuu(:,:,k))));
        end
    end
end

% Construct function output
fmax = Fs/2;
fi = find(mf <= fmax);
w = mf(fi).*2*pi;
G = G(:,:,fi);
if nu == ny && nargout == 3
    Coh = Coh(:,:,fi);
end
end

function out = freqAvg(in,nrbands)
%FREQAVG Frequency averaging
%  out=freqAvg(in,nrbands) averages the input 'in' over the number of
%  frequency bands 'nrbands'.

N      = length(in);
Nmod   = floor((N/nrbands));   % number of remaining frequencies after averaging
tmp    = zeros(nrbands,Nmod);  % initialization of temporary matrix for averaging: nrband rows and nmod columns
tmp(:) = in(1:nrbands*Nmod);   % arrange the samples of 'in' in the elements of tmp.
out    = mean(tmp,1)';         % average over columns and make it a vector
end