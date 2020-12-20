function [A,B,C,D,x0,K,output] = doptlti(u,y,A,B,C,D,x0,K,partype,options,sigman,filtera) 
%DOPTLTI   Performs a Least Squares optimization of a discrete 
%          time linear state space system system with model 
%          structure: 
% 
%              x(k+1) = Ax(k) + Bu(k) + Ke(k) 
%              y(k)   = Cx(k) + Du(k) +  e(k) 
% 
%          First, the state space matrices are parameterized. 
%          The Output Normal parametrization, the Tridiagonal 
%          parametrization and the Full parametrization can be 
%          used. 
%          The parameterized model is optimized using the supplied 
%          lmmore Levenberg-Marquardt function. The matrices A,B, 
%          and C are always returned. If needed, D, the initial 
%          state and a Kalman gain can also be optimized. 
% 
% Syntax: 
%          [A,B,C,D]=doptlti(u,y,A,B,C,D) 
%          [A,B,C,D,x0,K,options] = doptlti(u,y,A,B,C,D,x0,K,partype,options) 
% 
% Input: 
%   u,y        The input and output data of the system to be 
%              optimized. 
%   A,B,C,D    Initial estimates of the system matrices A, B, C and D. 
%   partype    This parameter specifies the type of parameterization 
%              that is used to parameterize the state space model. 
%              Three types of parameterization are supported: 
%              'on'= Output Normal, 'tr'=TRidiagonal and 
%              'fl'= FuLl. Default is 'fl'. 
%   options    Input parameters that are passed on directy to the 
%              optimization function. The options are compatible 
%              with those of the Optimization Toolbox. 
%              An extra field options.Manifold may be set to 'on' 
%              or 'off' if the full parametrization is used. The 
%              Manifold field indicates whether the search 
%              direction should be confined to directions in which 
%              the cost-function changes. 
%              If options.Manifold is not set, doptlti will set it 
%              to 'off' for the Output Normal and Tridiagonal 
%              parametrizations, and to 'on' for the Full 
%              parametrization. 
%              See foptions or optimset for more information 
%   sigman     The function of this parameters depends on its format: 
%              1. If sigman is a 1 x l vector, the output errors will 
%                 be weighted by the inverse of these factors. In a 
%                 weighted least squares estimation, sigman should 
%                 contain the standard deviations of the noise on each 
%                 of the outputs. 
%                 In a maximum likelihood estimation which assumes no 
%                 correlation between the noise on different outputs [1], 
%                 sigman should contain the standard deviations of the 
%                 white noise which, when put through the AR filter 
%                 specified by 'filtera', generates the output measure- 
%                 ment noise. 
%              2. If sigman is a l x 2l matrix, a maximum likelihood 
%                 estimation which does support correlation between the 
%                 output noises will be carried out. The 'filtera' 
%                 parameter MUST be specified in this case. 
%                 sigman should be [Sf Sb], in which Sf is the covariance 
%                 matrix of the multivariable white noise sequence that 
%                 is put through the causal filter Af (see 'filtera'). 
%                 Sb is the covariance matrix of the white noise sequence 
%                 that will be put through the anticausal filter Ab [2]. 
% 
%   filtera    The specification of the AR noise model. This should be 
%              either a matrix of size order x l, or a matrix of size 
%              2l x l*order. 
%              In the first output case filtera should be a matrix having 
%              max(d_i)+1 rows and l columns. If a certain output noise 
%              model has a lower order then pad the coefficient vector with 
%              NaNs [1]. 
%              In the second case, filtera should be [Af;Ab] in which Af 
%              specifies the causal AR filter, and Ab specifies the anti- 
%              causal AR filter [2]. 
% 
% Output: 
%   A,B,C,D    System matrices of the optimized linear model. 
%              If the D matrix is not estimated, it will be empty. 
%   x0         Estimate of the initial state. If the x0 matrix is 
%              not estimated, it will be returned empty. 
%   K          Estimate of the Kalman gain. If the K matrix is 
%              not estimated, it will be returned empty. 
%   options    Output parameters from the Optimization Toolbox. 
%              See foptions or optimset. 
% 
% Large-sample size optimizations: 
%              For a very large number of samples, memory-usage can 
%              be reduced by building the Jacobian Nb block-rows 
%              at a time [3]. To this end, the following options 
%              should be set: 
% 
%              options.RFactor='on' 
%              options.LargeScale='on' 
%              options.BlockSize=Nb 
% 
% Unstable systems: 
%              When the initial system model in unstable, one can 
%              specify a fixed Kalman K that causes the initial 
%              (A-KC) to be stable. K can be passed as a parameter 
%              as usual, and options.OEMStable='on' should be passed. 
%              This way, the K-factor is fixed, which speeds up 
%              convergence compared to a pure Prediction Error 
%              optimization. 
% 
% Limitations: 
%              The Levenberg-Marquardt optimization used by this 
%              function is a nonlinear gradient search. This implies 
%              that there is the inherent risk of ending up in a local 
%              minimum. Therefore, if the optimization does not converge, 
%              or converges to a local minimum, another initial estimate 
%              should be tried. 
% 
% Remarks: 
%              This optimization function has been targeted at 
%              MATLAB version 6 or higher. However, the function 
%              will run on MATLAB version 5 using a compatibility 
%              kludge. 
%              This kludge implies that the options input parameter 
%              can either be a MATLAB 6 'optimset'-structure, or 
%              a MATLAB 5 compatible 'foptions'-vector. However, 
%              the latter is discouraged since it does not allow 
%              the Manifold-field to be set. 
% 
% See also:    lmmore, dfunlti, optimset/foptions 
% 
% References: 
%  [1]      Benoit David and Georges Bastin, 'An estimator of the inverse 
%           covariance matrix and its application to ML parameter estimation 
%           in dynamical systems', Automatica, volume 37, number 1, 2001, 
%           pages 99-106. 
%  [2]      Benoit David, 'Parameter Estimation in Nonlinear Dynamical 
%           Systems with Correlated Noise', Ph.D. thesis, Universite 
%           Catholique de Louvain, November 2001. 
%  [3]      Niek Bergboer, Vincent Verdult and Michel Verhaegen, 
%           'An Efficient Implementation of Maximum Likelihood Identification 
%           of LTI State-Space Models by Local Gradient Search', submitted 
%           to the Conference on Decision and Control 2002. 
 
% This function is closely based on dslslin.m by 
% B. Haverkamp and V. Verdult, August 1998 
%
% Revised by Niek Bergboer, 2001 
% Revised by Ivo Houtzager, 2007
% Copyright (c) 1994-2007, Delft Center of Systems and Control
 
% Check number of arguments
if nargin < 5
    error('DOPTLTI requires at least five input arguments')
end

% Check the number of input parameters and set default values
if nargin<12,
    filtera = [];
end;
if nargin<11,
    sigman = [];
end;
if nargin<10 || isempty(options),
    options = mkoptstruc;
end
if nargin<9 || isempty(partype),
    partype = 'fl';
end
if nargin<8
    K = [];
end
if nargin<7,
    x0 = [];
end
if nargin<6,
    D = [];
end

if strcmp(partype,'on')
    params.partype = partype;
elseif strcmp(partype,'tr')
    params.partype = partype;
elseif strcmp(partype,'fl')
    params.partype = partype;
else
    error('You specified an unknown type of parameterization')
end

% Check the optimization options: these may be either empty,
% foptions, or optimset compatible.
if ~isstruct(options) && ~isempty(options),
    % options is not a struct: check whether it's a MATLAB 5
    % foptions vector
    if ~(all(size(options) == [1 18])),
        % No, wrong size: bail out
        error('options is not an optimset- or foption-compatible structure');
    else
        % Size of options-vector consistent with foptions
        warning('LTI:matlab5','Using MATLAB 5 compatibility kludge');

        % Translate information from old options-vector to
        % new optimset structure
        options = optim5to6(options);
    end
end

if isempty(options)
    % Set default struct
    options = mkoptstruc;
end

% At this point, an optimset-structure is available
% Set options.Manifold to a default value is not specified
if isfield(options,'Manifold'),
    if isempty(options.Manifold),
        if strcmp(params.partype,'fl')
            options.Manifold = 'on';
        else
            options.Manifold = 'off';
        end
    end
else
    if strcmp(params.partype,'fl')
        options.Manifold = 'on';
    else
        options.Manifold = 'off';
    end
end

% Enforce Manifold='on' for the Full parametrization
if strcmp(params.partype,'fl') && ~strcmp(options.Manifold,'on')
    error('Manifold must be on for Full parametrization');
end

% Check whether the options.RFactor, options.BlockSize and
% options.OEMStable fields exist. And set them.
if ~isfield(options,'RFactor')
    options.RFactor = [];
end
if ~isfield(options,'BlockSize')
    options.BlockSize = 0;
end
if ~isfield(options,'OEMStable')
    options.OEMStable = [];
end

% If the blocksize has been set to something >0 (and thus other than
% the default), set options.LargeScale to 'on'
if options.BlockSize>0 && isempty(options.LargeScale)
    warning('LTI:resetOptions','Resetting options.LargeScale from default to ''on''');
    options.LargeScale = 'on';
end

% Check whether the large-scale optimization has been selected,
% and enforce RFactor='on'; in that case.
% The overparametrized (and unprojected) Tridiagonal parametrization

% does not work together with the LargeScale framework.
if strcmp(options.LargeScale,'on')
    if isempty(options.RFactor),
        warning('LTI:resetOptions','Resetting options.RFactor from default to ''on''');
        options.RFactor = 'on';
    end
    if ~strcmp(options.RFactor,'on')
        error('Large-scale optimizations require OPTIONS.RFactor=''on''');
    end
    if strcmp(params.partype,'tr')
        error('Large-scale optimizations are incompatible with the tridiagonal parametrization');
    end
end

if size(y,2)>size(y,1)
    y  =  y';
end
if size(u,2)>size(u,1)
    u  =  u';
end

m = size(u,2);
l = size(y,2);
N = size(u,1);
params.N = N;
params.m = m;
params.l = l;

if l == 0
    error('We need an output')
end
if m == 0
    error('We need an input')
end

% Set flags for various matrices
fB = ~all(size(B) == 0);
fD = ~all(size(D) == 0);
fx = ~isempty(x0);
params.fB  =  fB;
params.fD  =  fD;
params.fx  =  fx;

% K is a special case: if options.OEMStable='on', then fK is
% zero even though K is passed as a parameter. In that case,
% we use the method by Chou and Verhaegen.
if ~all(size(K) == 0),
    if strcmp(options.OEMStable,'on'),
        fK = 0;
        % Note: options.K will be set later, when the similarity
        % transformation is known
    else
        fK = 1;
    end;
else
    if strcmp(options.OEMStable,'on'),
        error('K MUST be passed if options.OEMStable is on');
    else
        fK = 0;
    end;
end;
params.fK = fK;

% Check sizes and properties of ABCD
n = size(A,1);
params.n = n;
if size(A,2) ~= n
    error('A must be square');
end
if size(C,2) ~= n
    error('The size of C does not correspond to the number of states');
end
if size(C,1) ~= l
    error('The size of C does not correspond to the number of outputs')
end
if ~fB
    error('This optimization does not work without a B matrix');
end;
if size(B,2)~=m
    error('The size of B does not correspond to the number of inputs')
end;
if fD
    if ~all(size(D)==[l m])
        error('The size of D is incorrect');
    end
end
if fx
    if ~all(size(x0)==[n 1])
        error('The size of x0 is incorrect');
    end
end
if fK || strcmp(options.OEMStable,'on')
    if ~all(size(K)==[n l]),
        error('The size of K is incorrect');
    end
end

% Check stability and observability of the initial model
if ~fK
    % Output error optimization
    if strcmp(options.OEMStable,'on'),
        % However, in this case K is used to stabilize A-K*C
        if max(abs(eig(A-K * C)))>1
            error('Initial stabilized A-K*C must be stable.')
        end
        Ac = A - K*C;
        Cc = C;
    else
        % The normal OEM case
        if max(abs(eig(A)))>1
            error('Initial A must be stable.')
        end
        Ac = A; 
        Cc = C;
    end
else
    % The PEM case
    if (max(abs(eig(A-K * C)))>1)
        error('Initial A-K*C must be stable.')
    end
    Ac = A-K*C;
    Cc = C;
end

% Determine observability matrix
ObMat = zeros(n*l,n);
ObMat(1:l,:) = Cc;
for blockrow = 2:n
    ObMat((blockrow-1)*n+1:blockrow*n,:) = ObMat((blockrow-2)*n+1:(blockrow-1)*n,:)*Ac;
end
if rank(ObMat)<n
    error('Initial model must be observable');
end

if ~isempty(sigman),
    if all(size(sigman)==[1 l])
        if ~isempty(filtera)
            OptType = 'uncorr';
        else
            OptType = 'no_mle';
        end
    elseif all(size(sigman)==[l 2*l])
        OptType = 'flcorr';
    else
        error('The sigman parameter has an illegal size');
    end
else
    if ~isempty(filtera)
        OptType = 'uncorr';
    else
        OptType = 'no_mle';
    end
end

% Check conflict with PEM or OEMStable
if (strcmp(OptType,'uncorr') || strcmp(OptType,'flcorr')) && (fK || strcmp(options.OEMStable,'on'))
    error('PEM or stabilized OEM, and maximum likelihood optimization, are mutually exclusive.');
end

% Check the AR noise filter parameters for the Maximum Likelihood
% optimization.
if ~isempty(filtera),
    % First determine if an uncorrelated ML estimation according to [1],
    % or a correlated according to [2] is to be performed.

    if strcmp(OptType,'no_mle') || strcmp(OptType,'uncorr')
        % Either weighted least squares or an uncorrelated ML estimation
        % according to [1].
        if size(filtera,2)~=l
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
                        warning('LTI:noiseFilter','Noise filter for output %i is not used; zero order.',i);
                    end
                end
            end

            % Check stability
            for i = 1:l
                if max(abs(roots(filtera(1:filtorder(i)+1,i)))) > 1 - 1e-9,
                    error('Noise filter for output %i is insufficiently stable',i);
                end
            end
        end
    else
        if size(filtera,1)~=2*l,
            error('filtera should be of size 2*l x order*l');
        end
        cicm = cholicm(filtera(1:l,:),filtera(l+1:2*l,:),sigman(:,1:l),sigman(:,l+1:2*l),N);
    end
else
    % Filter is empty
    if strcmp(OptType,'flcorr')
        % A filter MUST be specified in the correlated ML estimation case.
        error('filtera must be specified in a correlated maximum likelihood estimation');
    end
end

if ~isempty(filtera) && strcmp(OptType,'uncorr') && max(filtorder)>0
    % Create 3D array that holds the l correction matrices
    CorrD = zeros(max(filtorder),max(filtorder),l);

    % Loop over the outputs
    for i = 1:l
        % Only correct when output i has a noise filter
        if filtorder(i)>0,
            % Calculate a correction matrix: see internal document
            d = filtorder(i);

            % Lower right corner of M
            M2 = toeplitz([filtera(1,i);zeros(d-1,1)],filtera(1:d,i));

            % Calculate upper left corner of V'*V
            V = toeplitz(filtera(d+1:-1:2,i),[filtera(d+1,i),zeros(1,d-1)]);
            VTV = V' * V;

            % Calculate M2+D2
            MpD = chol(M2' * M2-VTV(d:-1:1,d:-1:1));

            % Calculate D2 and store in correction array
            CorrD(1:d,1:d,i) = MpD-M2;
        end
    end
else
    CorrD = [];
end

% PEM Conversion of the initial model
if fK || strcmp(options.OEMStable,'on'),
    if strcmp(params.partype,'on') || strcmp(params.partype,'tr')
        A = A-K*C;
        if fD,
            B = B-K*D;
        end
    end
    % from now on, A and B are the innovation-form A and B
    % for the Output Normal and Tridiagonal parametrizations
end

if fK
    [th,params,T] = dss2th(A,B,C,D,x0,K,params.partype);
else
    [th,params,T] = dss2th(A,B,C,D,x0,[],params.partype);
end

if strcmp(options.OEMStable,'on'),
    % Transform K according to the same similarity transformation
    options.K = T \ K;
end;

if strcmp(OptType,'uncorr') || strcmp(OptType,'no_mle')
    % Either a (Weighted) Least Squares optimization or an uncorrelated
    % maximum likelihood optimization.
    [th,resnorm,residual,exitflag,output] = lmmore('dfunlti',th,[],[],options,u,y,params,options,OptType,sigman,filtera,CorrD);
else
    % A correlated maximum likelihood optimization.
    [th,resnorm,residual,exitflag,output] = lmmore('dfunlti',th,[],[],options,u,y,params,options,OptType,cicm);
end

[A,B,C,D,x0,K] = dth2ss(th,params);

% Set possibly transformed K
if strcmp(options.OEMStable,'on')
    K = options.K;
end

if isempty(A) && (n~=0)
    warning('LTI:notStable','Optimization did not result in stable model')
    B = [];
    C = [];
    D = [];
    x0 = [];
    K = [];
else
    if (fK || strcmp(options.OEMStable,'on')) && (strcmp(partype,'on') || strcmp(partype,'tr'))
        % Convert model back into normal form
        % Once again: not for the Full parametrization
        A = A + K*C;
        if fD
            B = B + K*D;
        end
        % now A and B are the real A and B again
    end
end





