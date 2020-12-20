function [epsilon,psi,U2] = dfunlti(th,u,y,params,options,OptType,sigman,filtera,CorrD) 
%DFUNLTI   This function implements the cost-fuction for doptlti 
%          It is not meant for standalone use. 
% 
% Syntax: 
%          [epsilon]=dfunlti(th,u,y,params) 
%          [epsilon,psi]=dfunlti(th,u,y,params) 
%          [epsilon,psi,U2]=dfunlti(th,u,y,params) 
% 
% Inputs: 
%  th      Parameter vector describing the system 
%  u,y     The input and output data of the system to be optimized. 
%  params  A structure that contains the dimension parameters of 
%          the system, such as the order, the number of inputs, 
%          whether D, x0 or K is present in the model. 
%  options (optional) An optimset compatible options-structure. The 
%          fields options.RFactor, options.LargeScale, options.Manifold 
%          and options.BlockSize should have been added by \verb$doptlti$. 
%  OptType (optional) Indicates what kind of weighted least squares or 
%          maximum likelihood optimization is being performed: 
%          'no_mle' implies a nonlinear (weighted) least squares 
%          optimization. 
%          'uncorr' implies a maximum likelihood optimization without 
%          correlation among the output perturbances. 
%          'flcorr' implies a maximum likelihood optimization with 
%          correlated output perturbances. 
%  sigman  (optional) If OptType is \verb$'no_mle'$, this can be a vector of 
%          size 1 x l that indicates the standard deviation of the perturbance 
%          of each of the outputs (i.e. inverse weights for the outputs). 
%          If OptType is \verb$'uncorr'$, this should be a vector of size 
%          1 x l that indicates the standard deviation of the white noise 
%          innovations for the output perturbance AR model. 
%          If OptType is 'flcorr', this should be a Cholesky factor of the 
%          AR process' inverse covariance matrix, as obtained by cholicm. 
%  filtera (optional) The a-polynomial of an AR noise model. The first 
%          element should be 1, and the other elements should 
%          be d filter coefficients. In the multi-output case 
%          filtera should be a matrix having max(d_i)+1 rows and l 
%          columns. If a certain output noise model has a lower order 
%          then pad the coefficient vector with NaNs. 
%  CorrD   A d x d x l array that describes the D2 correction matrix 
%          for every output. If for output j a filter 
%          order k < d is used, then only CorrD(1:k,1:k,j) is relevant. 
%          No detail will be provided here. 
% 
% Outputs: 
%  epsilon Output of the costfunction, which is the square of the error 
%          between the output and the estimated output. 
%  psi     Jacobian of epsilon 
%  U2      Left null-space of Manifold matrix for the Full 
%          parametrization[1] 
% 
% See also: doptlti, ffunlti 
 
% References 
% 
% [1]      Vincent Verdult and Michel Verhaegen, 'Identification of 
%          Multivariable LPV State Space Systems by Local Gradient 
%          Search', conference article for the European Control 
%          Conference 2001 
% 
% Written by Niek Bergboer, December 2001 
% Revised by Ivo Houtzager, 2007
% Copyright (c) 2001-2007, Delft Center of Systems and Control 

% Assign default values to unspecified parameters
if nargin < 9
    CorrD = [];
end
if nargin < 8
    filtera = [];
end
if nargin < 7
    sigman = [];
end
if nargin < 6
    % Default optimization type
    OptType = 'no_mle';
end
if nargin < 5
    % Default options if none specified
    options = mkoptstruc;
    options.RFactor = 'off';
    options.BlockSize = 0;
    options.OEMStable = 'off';
end
if nargin < 4
    error('DFUNLTI requires at least four input arguments.');
end

N = size(u,1);
l = params.l;
n = params.n;
m = params.m;
fB = params.fB;
fD = params.fD;
fx = params.fx;
fK = params.fK;

if ~fB
    error('This optimization does not work without a B matrix');
end

if strcmp(OptType,'no_mle')
    if ~isempty(sigman)
        if ~all(size(sigman)==[1 l])
            error('The number of elements in sigman should equal the number of outputs');
        end
    end

    if ~isempty(filtera)
        error('filtera must not be specified in a (weighted) least squares optimization');
    end
    if ~isempty(CorrD)
        error('CorrD must not be specified in a (weighted) least squares optimization');
    end
    maxfiltorder = 0;

elseif strcmp(OptType,'uncorr')
    % An uncorrelated maximum likelihood optimization
    if ~isempty(sigman)
        if ~all(size(sigman)==[1 l])
            error('The number of elements in sigman should equal the number of outputs');
        end
    end
    if isempty(filtera) || isempty(CorrD)
        error('filtera and CorrD must be specified in an uncorrelated maximum likelihood optimization');
    end

    % Parse the contents of filtera
    if size(filtera,2) ~= l
        error('filtera and y should have the same number of columns');
    end
    maxfiltorder = size(filtera,1)-1;
    if maxfiltorder==0;
        warning('LTI:maxFilterOrder','Filter order should at least be one. Continuing without filter.');
        filtera = [];
    else

        % build vector with filter orders
        filtorder = zeros(l,1);
        for i = 1:l
            if isempty(find(isnan(filtera(:,i)),1)),
                % Maximum order
                filtorder(i,1) = maxfiltorder;
            else
                % Lesser order
                filtorder(i,1) = find(isnan(filtera(:,i)),1)-2;
                if filtorder(i,1)<1,
                    filtorder(i,1) = 0;
                end
            end
        end
    end

    % Check the size of CorrD
    if size(CorrD,1)~=maxfiltorder || size(CorrD,2)~=maxfiltorder || size(CorrD,3)~=l
        error('CorrD must be of size d x d x l');
    end

elseif strcmp(OptType,'flcorr')
    if isempty(sigman)
        error('cicm must be specified in a correlated maximum likelihood optimization');
    end
    cicm = sigman;

    sigman = [];
    if size(cicm,2)~=N*l
        error('cicm must have N*l columns');
    end

    % Convert cicm from LAPACK/BLAS band storage format into sparse format
    maxfiltorder = size(cicm,1)-1; % Required for filtering later
    cicm = spdiags(cicm',size(cicm,1)-1:-1:0,N*l,N*l);

else
    error('Unknown optimization type');
end

thn0 = th(1:n*l);
invalid_theta_flag = 0;
invalid_theta_map = zeros(n*l,1);
if strcmp(params.partype,'on')
    for i = 1:n
        si = th(l*(n-i)+1:l*(n-i+1));
        ti = si'*si;
        if ti>1-eps
            warning('LTI:invalidTheta','invalid theta found in dfunlti')
            th(l*(n-i)+1:l*(n-i+1)) = 0.99*si/sqrt(ti);
            invalid_theta_flag = 1;
            invalid_theta_map(l*(n-i)+1:l*(n-i+1)) = ones(l,1);
        end
    end
end

[A,B,C,D,x0,K] = dth2ss(th,params);
if strcmp(options.OEMStable,'on')
    K = options.K;
end

if ~fD
    D = zeros(l,m); % To avoid Matrix - [] = error, later on
end
if ~fx
    x0 = zeros(n,1);
end

% Calculate the estimated state trajectory
if fK || strcmp(options.OEMStable,'on')
    if strcmp(params.partype,'fl')
        xe = ltiitr(A-K*C,[B-K*D,K],[u,y],[],x0);
    else
        xe = ltiitr(A,[B,K],[u,y],[],x0);
    end;
else
    xe = ltiitr(A,B,u,[],x0);
end

% Calculate the estimated output
ye = C*xe'+D*u';

% Calculate the error-vector
epsilon = y'-ye;
epsilon = epsilon(:);

if invalid_theta_flag
    epsilon  =  epsilon+mean(std(y))*norm(th(1:n*l)-thn0);
end

epsilon(~isfinite(epsilon)) = 1e+40;
epsilon(epsilon>1e+40) = 1e+40;
epsilon(epsilon<-1e+40) = -1e+40;

if strcmp(OptType,'uncorr')
    % Make a maximum likelihood correction
    for i = 1:l
        % Adjust data for output i
        if filtorder(i)>0,
            d = filtorder(i);

            % Store filtered vector (flipped)
            vectemp = filter(filtera(1:d+1,i),1,epsilon((N-1) * l+i:-l:i,1));

            % Edit top part of error for output i
            epsilon(i:l:(N-d-1) * l+i,1) = vectemp(N:-1:d+1,1);

            % Edit bottom part
            epsilon((N-d) * l+i:l:(N-1) * l+i,1) = vectemp(d:-1:1,1)+CorrD(1:d,1:d,i) * epsilon((N-d) * l+i:l:(N-1) * l+i,1);
        end
    end
elseif strcmp(OptType,'flcorr')
    % Make a correlated maximum likelihood correction
    epsilon = cicm*epsilon;
end

if ~isempty(sigman)
    for i = 1:l
        epsilon(i:l:(N-1)*l+i,1) = epsilon(i:l:(N-1)*l+i,1)/sigman(i);
    end
end

if nargout >= 3
    if strcmp(params.partype,'fl')
        if fK
            U2 = simlns(A,B,C,K,fD,fx);
        else    
            U2 = simlns(A,B,C,[],fD,fx);
        end
    else
        error('Cannot return a manifold left null-space for current parametrization');
    end
end

if nargout >= 2
    LargeScale = strcmp(options.LargeScale,'on');
    RFactor = strcmp(options.RFactor,'on');
    if LargeScale
        BlockSize = options.BlockSize;
    else
        BlockSize = N;
    end % if LargeScale

    % Number of runs to complete R-factor
    NumberOfRuns = ceil((N-maxfiltorder)/(BlockSize-maxfiltorder));

    % Ns is number of samples in this run
    %Ns = BlockSize;
    Nb = 1;     % Starting at row 1
    if strcmp(params.partype,'fl')
        XInit = zeros(n,size(U2,2));
    else
        XInit = zeros(n,length(th));
    end; % if params.partyle=='fl'

    % Check whether the block-size is not too small
    MinBlockSize = ceil((size(XInit,2))/l)+maxfiltorder+1;
    if BlockSize<MinBlockSize,
        fprintf('OPTIONS.BlockSize must be at least %d for the current problem\n',MinBlockSize);
        error('Invalid OPTIONS.BlockSize');
    end;

    % Check whether the block-size is not too large
    if BlockSize>N,
        fprintf('OPTIONS.BlockSize must be <= N (=%d)\n',N);
        error('Invalid OPTIONS.BlockSize');
    end;

    % Set R to its initial value (empty)
    R = [];

    for JacRun = 1:NumberOfRuns
        Ns = min(BlockSize,N-Nb+1);  % Set number of rows
        Ne = Nb+Ns-1;    % Set end-row

        if strcmp(params.partype,'on')
            psi = zeros(Ns*l,length(th)+RFactor);
            psi_col = 0;
            % First compute A/C-related columns
            for i = 1:n*l
                [dA,dC] = dth2dac(th,params,i);
                X = ltiitr(A,dA,xe(Nb:Ne,:),[],XInit(:,psi_col+i));
                dye = -C*X' - dC*xe(Nb:Ne,:)';
                psi(:,psi_col+i) = dye(:);
                % Save state and shift forward

                XInit(:,psi_col+i) = X(Ns-maxfiltorder,:)';
                XInit(:,psi_col+i) = A*XInit(:,psi_col+i) + dA*xe(Ne-maxfiltorder,:)';
            end
            psi_col = psi_col + n*l;

            % B-related columns
            for i = 1:m*n
                db = zeros(n*m,1);
                db(i) = 1;
                dB = zeros(n,m);
                dB(:) = db;
                X = ltiitr(A,dB,u(Nb:Ne,:),[],XInit(:,psi_col+i));
                dye = -C*X';
                psi(:,psi_col+i) = dye(:);

                % Save state and shift forward
                XInit(:,psi_col+i) = X(Ns-maxfiltorder,:)';
                XInit(:,psi_col+i) = A*XInit(:,psi_col+i) + dB*u(Ne-maxfiltorder,:)';
            end;
            psi_col = psi_col + m*n;

            if fD
                % D-related columns
                for i = 1:l
                    for j = 1:m
                        % Copy j^th column of u to the Jacobian-column
                        % corresponding to the j^th column of D. The rows
                        % in the Jacobian correspond to the i^th output
                        % (i^th row of D)

                        psi(i:l:(Ns-1)*l+i,psi_col+i+(j-1)*l) = -u(Nb:Ne,j);
                    end
                end
                psi_col = psi_col+l * m;
            end
            if fx
                % x0-related columns
                for i = 1:n
                    dx0 = zeros(n,1);
                    dx0(i) = 1;
                    if isempty(R),
                        X = ltiitr(A,zeros(n,m),zeros(Ns,m),[],dx0);
                    else
                        X = ltiitr(A,zeros(n,m),zeros(Ns,m),[],XInit(:,psi_col+i));
                    end
                    dye = -C*X';
                    psi(:,psi_col+i) = dye(:);

                    % Save state and shift forward
                    XInit(:,psi_col+i) = X(Ns-maxfiltorder,:)';
                    XInit(:,psi_col+i) = A * XInit(:,psi_col+i);
                end
                psi_col = psi_col+n;

            end
            if fK
                % K-related columns
                for i = 1:n*l
                    dk = zeros(n*l,1);
                    dk(i) = 1;
                    dK = zeros(n,l);
                    dK(:) = dk;
                    X = ltiitr(A,dK,y(Nb:Ne,:),[],XInit(:,psi_col+i));
                    dye = -C*X';
                    psi(:,psi_col+i) = dye(:);

                    % Save state and shift forward
                    XInit(:,psi_col+i) = X(Ns-maxfiltorder,:)';
                    XInit(:,psi_col+i) = A*XInit(:,psi_col+i)+dK*y(Ne-maxfiltorder,:)';
                end
                %psi_col = psi_col + n*l;
            end
            if invalid_theta_flag  % Correction in Jacobian for invalid theta
                for i = 1:n*l
                    if invalid_theta_map(i)
                        psi(:,i) = psi(:,i)+mean(std(y))*(th(i)-thn0(i))/norm(th(1:n * l)-thn0);
                    end
                end
            end

        elseif strcmp(params.partype,'fl')
            % Reserve a column in psi for epsilon if RFactor == 1
            psi = zeros(Ns*l,size(U2,2)+RFactor);
            for i = 1:size(U2,2)

                % Select column and set theta-deviation
                dth = U2(:,i);

                % Use dth2ss to extract delta matrices
                [dA,dB,dC,dD,dx0,dK] = dth2ss(dth,params);
                if strcmp(options.OEMStable,'on')
                    dK = zeros(n,l);
                end

                % Prevent error in PEM when D is not estimated
                if ~fD
                    fD = zeros(l,m);
                end

                % Determine whether the change is in D or x0; if it is, they
                % are the only quantities that change.

                if any(abs(dD(:))>1e-12),
                    % Only D changes, determine at which position
                    [Di,Dj] = find(dD>0); % Only one element is non-zero

                    % Copy input Dj the output Di
                    % OEM contribution
                    psi(Di:l:(Ns-1) * l+Di,i) = -u(Nb:Ne,Dj);

                    if fK || strcmp(options.OEMStable,'on')
                        % PEM (w/ Kalman Gain)
                        X = ltiitr(A-K*C,-K*dD,u(Nb:Ne,:),[],XInit(:,i));
                        thevec = -C * X';
                        psi(:,i) = psi(:,i)+thevec(:);

                        % Shift state forward (for large-scale problems)
                        XInit(:,i) = X(Ns-maxfiltorder,:)';
                        XInit(:,i) = (A-K*C)*XInit(:,i)-K*dD*u(Ne-maxfiltorder,:)';
                    end

                elseif any(dx0>0),
                    % Only x0 changes, determine at which position
                    %xi = find(dx0>0);  % Only one element is non-zero
                    % Simulate transient reponse from initial state

                    if fK || strcmp(options.OEMStable,'on'),
                        % Use a Kalman gain
                        if isempty(R), % First time
                            X = ltiitr(A-K*C,zeros(n,m),zeros(Ns,m),[],dx0);
                        else
                            X = ltiitr(A-K*C,zeros(n,m),zeros(Ns,m),[],XInit(:,i));
                        end;

                        % Shift state forward
                        XInit(:,i) = X(Ns-maxfiltorder,:)';
                        XInit(:,i) = (A-K*C) * XInit(:,i);
                    else % if fK
                        % The normal output-error case
                        if isempty(R), % First time
                            X = ltiitr(A,zeros(n,m),zeros(Ns,m),[],dx0);
                        else
                            X = ltiitr(A,zeros(n,m),zeros(Ns,m),[],XInit(:,i));
                        end;

                        % Shift state forward
                        XInit(:,i) = X(Ns-maxfiltorder,:)';
                        XInit(:,i) = A*XInit(:,i);
                    end; % if fK
                    thevec = -C*X';
                    psi(:,i) = thevec(:);
                    
                else
                    % D and x0 do not change, so all the other matrices do. 
                    % Calculate the full Jacobian column.
                    %thevec = zeros(N * l,1);

                    % Calculate the Jacobian column
                    if fK || strcmp(options.OEMStable,'on')
                        % PEM (w/ Kalman Gain)
                        dA_KC = dA - dK*C - K*dC;
                        dB_KD = dB - dK*D; %-K*dD;
                        X = ltiitr(A - K*C,[dA_KC dB_KD dK],[xe(Nb:Ne,:) u(Nb:Ne,:) y(Nb:Ne,:)],[],XInit(:,i));
                        thevec = -C*X' - dC*xe(Nb:Ne,:)';

                        % Shift state forward (for large-scale problems)
                        XInit(:,i) = X(Ns-maxfiltorder,:)';
                        XInit(:,i) = (A-K*C)*XInit(:,i)+[dA_KC dB_KD dK]*[xe(Ne-maxfiltorder,:) u(Ne-maxfiltorder,:) y(Ne-maxfiltorder,:)]';
                    else
                        % OEM (w/o Kalman Gain)
                        X = ltiitr(A,[dA dB],[xe(Nb:Ne,:) u(Nb:Ne,:)],[],XInit(:,i));
                        thevec = -C*X' - dC*xe(Nb:Ne,:)';

                        % Shift state forward (for large-scale problems)
                        XInit(:,i) = X(Ns-maxfiltorder,:)';
                        XInit(:,i) = A*XInit(:,i) + [dA dB]*[xe(Ne-maxfiltorder,:) u(Ne-maxfiltorder,:)]';
                    end
                    psi(:,i) = thevec(:);
                end
            end


        elseif strcmp(params.partype,'tr')
            % Tri-diagonal parametrization
            psi = zeros(Ns * l,length(th)+RFactor);
            psi_col = 0;
            
            % A-related columns
            for i = 1:n-1
                dA = zeros(n,n);
                dA(i,i+1) = 1; % Upper diagonal
                X = ltiitr(A,dA,xe(Nb:Ne,:),[],XInit(:,psi_col+i));
                dye = -C*X';
                psi(:,psi_col+i) = dye(:);
                % Save state and shift forward

                XInit(:,psi_col+i) = X(Ns-maxfiltorder,:)';
                XInit(:,psi_col+i) = A*XInit(:,psi_col+i) + dA*xe(Ne-maxfiltorder,:)';
            end
            
            psi_col = psi_col+n - (n>0) * 1;
            
            for i = 1:n
                dA = zeros(n,n);
                dA(i,i) = 1; % Diagonal
                X = ltiitr(A,dA,xe(Nb:Ne,:),[],XInit(:,psi_col+i));
                dye = -C * X';
                psi(:,psi_col+i) = dye(:);
                
                % Save state and shift forward
                XInit(:,psi_col+i) = X(Ns-maxfiltorder,:)';
                XInit(:,psi_col+i) = A*XInit(:,psi_col+i) + dA*xe(Ne-maxfiltorder,:)';
            end
            psi_col = psi_col + n;
            
            for i = 1:n-1
                dA = zeros(n,n);
                dA(i+1,i) = 1; % Lower diagonal
                X = ltiitr(A,dA,xe(Nb:Ne,:),[],XInit(:,psi_col+i));
                dye = -C*X';
                psi(:,psi_col+i) = dye(:);
                
                % Save state and shift forward
                XInit(:,psi_col+i) = X(Ns-maxfiltorder,:)';
                XInit(:,psi_col+i) = A * XInit(:,psi_col+i)+dA * xe(Ne-maxfiltorder,:)';
            end
            psi_col = psi_col+n - (n>0)*1; % Also accomodate 0th order systems.

            % C-related columns
            for i = 1:l
                for j = 1:n
                    % Copy j^th column of xe to the Jacobian-column
                    % corresponding to the j^th column of C. The rows
                    % in the Jacobian correspond to the i^th output
                    % (i^th row of C)
                    psi(i:l:(Ns-1)*l+i,psi_col+i+l*j-1) = -xe(Nb:Ne,j);
                end
            end
            psi_col = psi_col + l*n;

            % B-related columns
            for i = 1:m*n
                db = zeros(n * m,1);
                db(i) = 1;
                dB = zeros(n,m);
                dB(:) = db;
                X = ltiitr(A,dB,u(Nb:Ne,:),[],XInit(:,psi_col+i));
                dye = -C*X';
                psi(:,psi_col+i) = dye(:);

                % Save state and shift forward
                XInit(:,psi_col+i) = X(Ns-maxfiltorder,:)';
                XInit(:,psi_col+i) = A*XInit(:,psi_col+i) + dB*u(Ne-maxfiltorder,:)';
            end
            psi_col = psi_col + m*n;

            if fD
                % D-related columns
                for i = 1:l
                    for j = 1:m
                        % Copy j^th column of u to the Jacobian-column
                        % corresponding to the j^th column of D. The rows
                        % in the Jacobian correspond to the i^th output
                        % (i^th row of D)
                        psi(i:l:(Ns-1)*l+i,psi_col+i+(j-1)*l) = -u(Nb:Ne,j);
                    end
                end
                psi_col = psi_col+l * m;
            end

            if fx
                % x0-related columns
                for i = 1:n
                    dx0 = zeros(n,1);
                    dx0(i) = 1;
                    if isempty(R)
                        X = ltiitr(A,zeros(n,m),zeros(Ns,m),[],dx0);
                    else
                        X = ltiitr(A,zeros(n,m),zeros(Ns,m),[],XInit(:,psi_col+i));
                    end;
                    dye = -C*X';
                    psi(:,psi_col+i) = dye(:);

                    % Save state and shift forward
                    XInit(:,psi_col+i) = X(Ns-maxfiltorder,:)';
                    XInit(:,psi_col+i) = A*XInit(:,psi_col+i);
                end;
                psi_col = psi_col+n;
            end
            
            if fK
                % K-related columns
                for i = 1:n*l
                    dk = zeros(n*l,1);
                    dk(i) = 1;
                    dK = zeros(n,l);
                    dK(:) = dk;
                    X = ltiitr(A,dK,y(Nb:Ne,:),[],XInit(:,psi_col+i));
                    dye = -C*X';
                    psi(:,psi_col+i) = dye(:);

                    % Save state and shift forward
                    XInit(:,psi_col+i) = X(Ns-maxfiltorder,:)';
                    XInit(:,psi_col+i) = A * XInit(:,psi_col+i)+dK * y(Ne-maxfiltorder,:)';
                end
                %psi_col = psi_col + n*l;
            end
        else
            error('Unknown parameter type');
        end

        psi(~isfinite(psi)) = 1e+40;
        psi(psi>1e+40) = 1e+40;
        psi(psi<-1e+40) = -1e+40;

        if strcmp(OptType,'uncorr')
            % Make a maximum likelihood correction for the Jacobian
            for i = 1:size(psi,2)-RFactor
                % Adjust column i
                for j = 1:l
                    % Adjust rows corresponding to output j
                    if filtorder(j)>0
                        d = filtorder(j);

                        % Create filtered vector (flipped)
                        vectemp = filter(filtera(1:d+1,j),1,psi((Ns-1)*l+j:-l:j,i));

                        % Depending on whether Ne == N, do or do not correct the bottom part
                        % of the psi-column

                        if Ne==N
                            % Edit top part of psi-column
                            psi(j:l:(Ns-d-1)*l+j,i) = vectemp(Ns:-1:d+1,1);

                            % Edit bottom part
                            psi((Ns-d) * l+j:l:(Ns-1)*l+j,i) = vectemp(d:-1:1,1)+CorrD(1:d,1:d,j)*psi((Ns-d)*l+j:l:(Ns-1)*l+j,i);
                        else
                            % Edit the psi-column
                            psi(j:l:(Ns-1)*l+j,i) = vectemp(Ns:-1:1,1);
                        end
                    end
                end
            end
        elseif strcmp(OptType,'flcorr')
            % Make a correlated maximum likelihood correction
            psi = cicm((Nb-1)*l+1:Ne*l,(Nb-1)*l+1:Ne*l)*psi;
        end;

        if ~isempty(sigman),
            % Make a weighting correction
            for i = 1:l
                psi(i:l:(Ns-1)*l+i,:) = psi(i:l:(Ns-1)*l+i,:)/sigman(i);
            end
        end

        % If only the R factor of [psi epsilon] is wanted, return it instead
        % of psi.
        if strcmp(options.RFactor,'on')
            % Add block-rows Nb to Ne-maxfiltorder
            % Add error-vector
            psi(:,size(psi,2)) = epsilon(l * (Nb-1)+1:l * Ne,1);

            % Truncate psi (optionally) if filtering takes place
            if maxfiltorder>0 && (Ne~=N)
                psi = psi(1:l * (Ns-maxfiltorder),:);
            end
            if isempty(R) % First time
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
        Nb = Nb+Ns-maxfiltorder;
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
% dA/dtheta_j   dC/dtheta_j
n = params.n;
l = params.l;
if strcmp(params.partype,'on')
    % The thnr'ed parameter resides in block ceil((n*l-thnr+1)/l)
    blocknr = ceil((n * l-thnr+1)/l);
    Il = eye(l);
    Z = [zeros(l,n);eye(n)];
    for i = 1:n
        if i~=blocknr,

            % Normal case: calculate rotation as usual
            si = th(l*(n-i)+1:l*(n-i+1));
            ti = si'*si;
            ri = sqrt(1-ti);
            Ti = eye(n+l);
            Ti(i:i+l,i:i+l) = [-si,Il-(1-ri)/ti*si* i'; ri,si'];
            Z = Ti'*Z;
            
        else % i == blocknr
            % This is the block that contains theta_j
            % Within this block, it is parameter nr. (thnr-1 mod l) + 1
            j = mod(thnr-1,l)+1;
            si = th(l*(n-i)+1:l*(n-i+1));
            thj = si(j);
            ti = si' * si;
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
        end
    end
    
    % Set the output values
    dA = Z(l+1:l+n,:);
    dC = Z(1:l,:);
else
    error('Unknown parameter type!');
end
end





