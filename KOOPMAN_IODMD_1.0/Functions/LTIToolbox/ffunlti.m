function [epsilon,psi,U2] = ffunlti(th,H,w,params,options,timing) 
%FFUNLTI   This function implements the costfuction for the foptlti 
%          frequency domain optimization function. 
%          It is not meant for standalone use. 
% 
% Syntax: 
%          [epsilon]=dfunlin(th,H,params,timing) 
%          [epsilon,psi]=dfunlin(th,H,params,timing) 
%          [epsilon,psi,U2]=dfunlin(th,H,params,timing) 
% 
% Inputs: 
%  th      Parameter vector describing the system 
%  H       The frequency response function of the system to be optimized. 
%          A matrix of size (l x m x N) in which H(:,:,i) contains the 
%          complex FRF at complex frequency number i. 
%  w       Complex frequencies at which the FRF is measured. 
%  params  A structure that contains the dimension parameters of 
%          the system, such as the order, the number of inputs, 
%          whether D, x0 or K is present in the model, etc. 
%  timing  Either 'cont' or 'disc', indicating that the supplied model 
%          is continuous of discrete time. Note that this influences ONLY 
%          the way in which the Output Normal parametrization is built. The 
%          user is still responsible for supplying suitable frequency data 
%          (j.omega or exp(j.omega)). 
% 
% Outputs: 
%  epsilon Output of the costfunction, which is the square of the error 
%          between the actual and predicted vectorized frequency response 
%          function. 
%  psi     Jacobian of epsilon 
%  U2      Left null-space of Manifold matrix for the full 
%          parametrization[1] 
% 
% See also: foptlti 
 
% References: 
% 
% [1]      Vincent Verdult and Michel Verhaegen, 'Idenfitication of 
%          Multivariable LPV State Space Systems by Local Gradient 
%          Search', conference article for the European Control 
%          Conference 2001 
 
% Written by Niek Bergboer, 2001 
% Revised by Ivo Houtzager, 2007
% Copyright (c) 2001-2007, Delft Center of Systems and Control

% Check number of arguments
if nargin < 4
    error('FFUNLTI requires at least four input arguments');
end
if nargin < 6  || isempty(timing)
    timing = 'disc';
end
if nargin < 5 || isempty(options)
    options = mkoptstruc;
    options.RFactor = 'off';
    options.BlockSize = 0;
end

N  =  size(H,3);
l  =  params.l;
n  =  params.n;
m  =  params.m;
fB  =  params.fB;
fD  =  params.fD;
fx  =  params.fx;
fK  =  params.fK;

if strcmp(timing,'disc'),
    TT = 'disc';
elseif strcmp(timing,'cont'),
    TT = 'cont';
else
    error('Timing should be either cont or disc');
end
if ~fB
    error('This optimization does not work without a B matrix');
end;
if fK
    error('Frequency-domain optimization does not work innovation models');
end
if fx
    error('Frequency-domain optimization does not work with an initial state');
end
if size(w,2) > size(w,1),
    w = w.';
end

thn0 = th(1:n*l);
invalid_theta_flag = 0;
invalid_theta_map = zeros(n * l,1);
if strcmp(TT,'disc')
    if strcmp(params.partype,'on')
        for i = 1:n
            si = th(l*(n-i)+1:l*(n-i+1));
            ti = si'*si;
            if ti > 1-eps
                warning('LTI:invalidTheta','invalid theta found in ffunlti')
                th(l*(n-i)+1:l*(n-i+1)) = 0.99*si/sqrt(ti);
                invalid_theta_flag = 1;
                invalid_theta_map(l*(n-i)+1:l*(n-i+1)) = ones(l,1);
            end
        end
    end
end

% Convert parameter-vector into a model
if strcmp(TT,'disc')
    [A,B,C,D] = dth2ss(th,params);
else
    [A,B,C,D] = cth2ss(th,params);
end
if ~fD
    D = zeros(l,m); % To avoid Matrix - [] = error, later on
end

%Compute FRF for current estimate
Hc = ltifrf(A,B,C,D,[],w,1);

% Calculate error-vector
epsilon = H(:)-Hc;
epsilon = [real(epsilon); imag(epsilon)];
if invalid_theta_flag
    epsilon  =  epsilon + sqrt(norm(H(:),2))*norm(th(1:n*l)-thn0);
end
epsilon(~isfinite(epsilon)) = 1e+40;
epsilon(epsilon>1e+40) = 1e+40;
epsilon(epsilon<-1e+40) = -1e+40;

if nargout >= 3,
    if strcmp(params.partype,'fl')
        U2 = simlns(A,B,C,[],fD,[]);
    else
        error('Cannot return a manifold left null-space for current parametrization');
    end
end

if nargout >= 2,
    LargeScale = strcmp(options.LargeScale,'on');
    RFactor = strcmp(options.RFactor,'on');
    if LargeScale,
        BlockSize = options.BlockSize;
    else
        BlockSize = N;
    end; % if LargeScale

    % Number of runs to complete R-factor
    NumberOfRuns = ceil(N/BlockSize);

    % Ns is number of samples in this run
    %Ns = BlockSize;
    Nb = 1;     % Starting at row 1

    % Check whether the block-size is not too small
    if strcmp(params.partype,'fl')
        MinBlockSize = ceil(size(U2,2)/(2 * l * m))+1;
    else
        MinBlockSize = ceil(length(th)/(2 * l * m))+1;
    end
    if BlockSize<MinBlockSize,
        fprintf('OPTIONS.BlockSize must be at least %d for the current problem\n',MinBlockSize);
        error('Invalid OPTIONS.BlockSize');
    end

    % Check whether the block-size is not too large
    if BlockSize>N,
        fprintf('OPTIONS.BlockSize must be <= N (=%d)\n',N);
        error('Invalid OPTIONS.BlockSize');
    end

    % Set R to its initial value (empty)
    R = [];

    for JacRun = 1:NumberOfRuns
        Ns = min(BlockSize,N-Nb+1);  % Set number of rows
        Ne = Nb+Ns-1;    % Set end-row

        % If a Jacobian is desired then compute it
        if strcmp(params.partype,'on')
            psic = zeros(Ns*l*m,length(th)+RFactor);
            psi_col = 0;

            % A/C-related
            for i = 1:n*l
                if strcmp(TT,'disc')
                    [dA,dC] = dth2dac(th,params,i);
                else
                    [dA,dC] = cth2dac(th,params,i);
                end
                psic(:,psi_col+i) = ltifrf(A,B,dC,[],[],w(Nb:Ne),1)+ltifrf(A,B,C,[],dA,w(Nb:Ne),1);
            end
            psi_col = psi_col+n*l;

            % B-related
            for i = 1:n
                dB = zeros(n,m);
                dB(i,1) = 1;
                b_column = ltifrf(A,dB,C,[],[],w(Nb:Ne),1);
                b_column = b_column(1:(Ns*l*m)-(m-1)*l,1);
                for j = 1:m
                    psic((j-1)*l+1:(Ns-1)*l*m+j*l,psi_col+i+(j-1)*n) = b_column;
                end
            end
            psi_col = psi_col+m * n;

            % D-related
            if fD
                for i = 1:l
                    for j = 1:m
                        psic((j-1)*l+i:m*l:(Ns-1)*m*l+(j-1)*l+i,psi_col+i+(j-1)*l) = ones(Ns,1);
                    end
                end
                %psi_col = psi_col+l * m;
            end
            if invalid_theta_flag  % Correction in Jacobian for invalid theta
                for i = 1:n*l
                    if invalid_theta_map(i)
                        psi(:,i) = psi(:,i)+sqrt(norm(H(:),2))*(th(i)-thn0(i))/norm(th(1:n * l)-thn0);
                    end
                end
            end

        elseif strcmp(params.partype,'fl')
            psic = zeros(Ns*l*m,size(U2,2)+RFactor);
            cachedB = zeros(Ns*l*m,n);
            cachedC = zeros(Ns*m,n);
            for i = 1:n
                % Calculate column for this state
                dB = zeros(n,m);
                dB(i,1) = 1;
                cachedB(:,i) = ltifrf(A,dB,C,[],[],w(Nb:Ne),1);
                dC = zeros(1,n);
                dC(1,i) = 1;
                cachedC(:,i) = ltifrf(A,B,dC,[],[],w(Nb:Ne),1);
            end

            % Cut cachedB to right size
            cachedB = cachedB(1:(Ns-1)*l*m+l,:);

            for colnum = 1:size(U2,2)
                % Get dtheta
                dth = U2(:,colnum);

                % Build matrices
                [dA,dB,dC,dD] = dth2ss(dth,params);

                if any(dD(:)>0)
                    % D-related
                    thevec = ltifrf([],[],[],dD,[],w(Nb:Ne),1);
                elseif n~=0
                    % If order>0, perform A,B and C-related calculations.
                    % dA-component: not-reducible
                    thevec = ltifrf(A,B,C,[],dA,w(Nb:Ne),1);

                    % dB-component: use cached copies
                    % per input
                    for j = 1:m
                        thevec((j-1)*l+1:(Ns-1)*l*m+j*l,:) = thevec((j-1)*l+1:(Ns-1)*l*m+j*l,:)+cachedB*dB(:,j);
                    end

                    % dC-component: use cached copies
                    % per output
                    for j = 1:l
                        thevec(j:l:(Ns-1)*l*m+(m-1)*l+j,:) = thevec(j:l:(Ns-1)*l*m+(m-1)*l+j,:)+cachedC*dC(j,:)';
                    end
                end
                % Store result
                psic(1:Ns*l*m,colnum) = thevec;
            end
            
        elseif strcmp(params.partype,'tr')
            % Tri-diagonal parametrization
            psic = zeros(Ns*l*m,length(th)+RFactor);
            psi_col = 0;

            % A-related
            for i = 1:n-1
                dA = zeros(n,n);
                dA(i,i+1) = 1;
                % Upper diagonal
                psic(:,psi_col+i) = ltifrf(A,B,C,[],dA,w(Nb:Ne),1);
            end
            psi_col = psi_col+n-(n>0) * 1;
            for i = 1:n
                dA = zeros(n,n);
                dA(i,i) = 1;
                % Diagonal
                psic(:,psi_col+i) = ltifrf(A,B,C,[],dA,w(Nb:Ne),1);
            end
            psi_col = psi_col+n;
            for i = 1:n-1
                dA = zeros(n,n);
                dA(i+1,i) = 1;
                % Lower diagonal
                psic(:,psi_col+i) = ltifrf(A,B,C,[],dA,w(Nb:Ne),1);
            end
            psi_col = psi_col+n-(n>0)*1;

            % C-related
            for j = 1:n
                dC = zeros(1,n);
                dC(1,j) = 1;
                c_column = ltifrf(A,B,dC,[],[],w(Nb:Ne),1);
                for i = 1:l
                    psic(i:l:Ns*l*m-l+i,psi_col+i+(j-1)*l) = c_column;
                end
            end
            if n~=0
                psi_col = psi_col+n * l;
            end

            % B-related
            for i = 1:n
                dB = zeros(n,m);
                dB(i,1) = 1;
                b_column = ltifrf(A,dB,C,[],[],w(Nb:Ne),1);
                b_column = b_column(1:(Ns*l*m)-(m-1) * l,1);
                for j = 1:m
                    psic((j-1)*l+1:(Ns-1)*l*m+j*l,psi_col+i+(j-1)*n) = b_column;
                end
            end
            psi_col = psi_col+m*n;

            % D-related
            if fD
                for i = 1:l
                    for j = 1:m
                        psic((j-1)*l+i:m*l:(Ns-1)*m*l+(j-1)*l+i,psi_col+i+(j-1)*l) = ones(Ns,1);
                    end;
                end;
                %psi_col = psi_col+l * m;
            end
        else
            error('Unknown parameter type');
        end
        psic(~isfinite(psic)) = 1e+40;
        psic(psic>1e+40) = 1e+40;
        psic(psic<-1e+40) = -1e+40;

        psic = -psic;
        psi = zeros(2 * size(psic,1),size(psic,2));
        psi(size(psic,1)+1:2 * size(psic,1),:) = imag(psic);
        psic = real(psic);
        psi(1:size(psic,1),:) = psic;
        clear psic;

        % If only the R factor of [psi epsilon] is wanted, return it instead
        % of psi.
        if strcmp(options.RFactor,'on'),
            % Add block-rows Nb to Ne
            % Add error-vector
            % First the real part
            psi(1:m*l*Ns,size(psi,2)) = epsilon(m*l*(Nb-1)+1:m*l*Ne,1);

            % And then the imaginary part
            psi(m*l*Ns+1:2*m*l*Ns,size(psi,2)) = epsilon(N*m*l+m*l*(Nb-1)+1:N*m*l+m*l*Ne,1);

            if isempty(R), % First time
                R = qr(psi);
                R = R(1:size(R,2),:);
                R = triu(R);
            else
                R = qr([R(1:size(R,2)-1,:);psi]);
                R = R(1:size(R,2),:);
                R = triu(R);
            end
        end

        % Update the start row for the next run (in the large scale case)
        Nb = Nb+Ns;
    end

    % Copy R to psi
    if RFactor
        psi = R;
    end
end
end


function [dA,dC] = dth2dac(th,params,thnr)
% Returns the matrices:
%
% dA/dtheta\_thnr   dC/dtheta\_thnr
%
% In which thnr is the parameternumber
%
% Method is based on the same Givens-rotation method used throughout
% the SMI-2.0 toolbox
%
% This currently works ONLY for the Output Normal Form and the
% full parametrization
n = params.n;
m = params.m;
l = params.l;
if strcmp(params.partype,'on')
    % The thnr'ed parameter resides in block ceil((n*l-thnr+1)/l)
    blocknr = ceil((n * l-thnr+1)/l);
    Il = eye(l);
    Z = [zeros(l,n);eye(n)];
    for i = 1:n
        if i~=blocknr,
            % Normal case: calculate Givens-rotation as usual
            si = th(l * (n-i)+1:l * (n-i+1));
            ti = si' * si;
            ri = sqrt(1-ti);
            Ti = eye(n+l);
            Ti(i:i+l,i:i+l) = [-si Il-(1-ri)/ti*(si*si'); ri si'];
            Z = Ti' * Z;
        else % i == blocknr
            % This is the block of l Givens-rotations that contains
            % the matrix difference dTi/dtheta
            % Within this block, it is parameter nr. (thnr-1 mod l) + 1
            j = mod(thnr-1,l)+1; % Number j in this block
            si = th(l*(n-i)+1:l*(n-i+1));
            thj = si(j);
            ti = si'*si;
            ri = sqrt(1-ti);
            Ti = zeros(n+l);

            % -si part (upper-left)
            Ti(i+j-1,i) = -1;

            % si' part (lower-right)
            Ti(i+l,i+j) = 1;

            % ri part (lower left)
            Ti(i+l,i) = -thj/ri;

            % Il - (1-ri)/ti * si * si' part (upper-right)
            % First component
            Ti(i:i+l-1,i+1:i+l) = -thj/(ri*ti)*(si*si');

            % Second component
            Ti(i:i+l-1,i+1:i+l) = Ti(i:i+l-1,i+1:i+l)-(1-ri)*(-2*thj/(ti^2))*(si*si');

            % Third component
            Mi = zeros(l,l);
            Mi(j,:) = si'; 
            Mi(:,j) = si; 
            Mi(j,j) = 2*thj;
            Ti(i:i+l-1,i+1:i+l) = Ti(i:i+l-1,i+1:i+l)-((1-ri)/ti)*Mi;

            % Do the transformation
            Z = Ti'*Z;
        end; % if i~=blocknr
    end; % for i=1:n
    % Set the output values
    dA = Z(l+1:l+n,:);
    dC = Z(1:l,:);

elseif strcmp(params.partype,'fl')
    % Full parametrization: one element of dA or dC is 1, the rest 0
    dth = zeros((n+l)*(n+m),1);
    dth(thnr) = 1;
    dM = zeros(n+l,n+m);
    dM(:) = dth;
    dA = dM(1:n,1:n);
    dC = dM(n+1:n+l,1:n);
else % if params.partype == 'on'
    error('Unknown parameter type!');
end
end

function [dA,dC] = cth2dac(th,params,thnr)
n = params.n;
l = params.l;

dth = zeros(n*l,1);
dth(thnr) = 1;
dC = zeros(l,n);
C = zeros(l,n);
dAss = zeros(n,n);
offset = 0;

% Build C, dC and dAss (following Haverkamp)
for j = 1:l
    C(j:l,j) = th(offset+1:offset+l-j+1);
    dC(j:l,j) = dth(offset+1:offset+l-j+1);
    offset = offset+l-j+1;
end
for j = 1:min(l,n-1)
    dAss = dAss+diag(dth(offset+1:offset+n-j),j)-diag(dth(offset+1:offset+n-j),-j);
    offset = offset+n-j;
end;

% Build dA/dtheta (dC/dtheta has already been built)
dA = -0.5*C'*dC-0.5*dC'*C+dAss;

end


