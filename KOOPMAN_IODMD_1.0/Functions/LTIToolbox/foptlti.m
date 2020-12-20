function [A,B,C,D,output] = foptlti(H,w,A,B,C,D,partype,options,timing) 
%FOPTLTI   Performs a Least Squares optimization of a discrete 
%          or continuous time linear state space system based on 
%          frequency reponse data. 
%          The model structure is the following: 
% 
%              x(k+1) = Ax(k) + Bu(k) 
%              y(k)   = Cx(k) + Du(k) 
% 
%          First, the state space matrices are parameterized. 
%          The Output Normal parametrization, the tridiagonal 
%          parametrization and the full parametrization can be 
%          used. 
%          The parameterized model is optimized using the supplied 
%          lmmore Levenberg-Marquardt function. The matrices A,B, 
%          C and D are always returned. 
% 
% Syntax: 
%          [A,B,C,D,options]=foptlti(H,w,A,B,C,D,model,partype,options) 
% Input: 
%   H          The measured frequency response function (FRF). This should 
%              be a matrix which follows the MATLAB 5 FRF storage convention: 
%              it should be a 3D-array of size l x m x N, in which H(:,:,j) 
%              is the FRF at the j-th complex frequency. 
% 
%   w          Vector is complex frequencies at which the FRF is measured. 
%              Although the function can operate using arbitrary complex 
%              frequencies, the following two choices are rather standard: 
% 
%              w = exp(j*omega)  For discrete-time systems 
%              w = j * omega     For continuous-time systems 
% 
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
%              If options.Manifold is not set, foptlti will set it 
%              to 'off' for the Output Normal and Tridiagonal 
%              parametrizations, and to 'on' for the Full 
%              parametrization. 
%              See foptions or optimset for more information 
%   timing     Must be either 'cont' or 'disc' to specify that 
%              the model if continuous or discrete time. Note 
%              that this changes ONLY the stability check and 
%              the Output Normal Form parametrization. It is 
%              still up to the user to supply suitable frequency 
%              data. 
% 
% Output: 
%   A,B,C,D    System matrices of the optimized linear model. 
%              If the D matrix is not estimated, it will be empty. 
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
% See also:    lmmore, ffunlti, optimset/foptions 
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
 
% N.H. Bergboer, October 2001 
% Revised by Ivo Houtzager, 2007
% Copyright (c) 2001-2007, Delft Center of Systems and Control

% Check number of arguments
if nargin < 5
    error('FOPTLI requires at least five input arguments')
end
if nargin < 9 || isempty(timing)
    timing = 'disc';
end
if nargin < 8 || isempty(options)
    options  =  mkoptstruc;
end
if nargin < 7 || isempty(partype)
    partype  =  'fl';
end
if nargin<6
    D = [];
end
if size(w,2)>size(w,1)
    w = w.';
end
if strcmp(timing,'disc')
    TT = 'disc';
elseif strcmp(timing,'cont')
    TT = 'cont';
else
    error('Timing should be either cont or disc');
end
N = size(w,1);
if ndims(H)~=3
    error('H should be l x m x N');
end
if size(H,3)~=N
    error('H should be l x m x N');
end

l  =  size(H,1);
m  =  size(H,2);
params.N  =  N;
params.m  =  m;
params.l  =  l;

if l == 0
    error('We need an output')
end
if m == 0
    error('We need an input')
end

fB = ~all(size(B) == 0);
fD = ~all(size(D) == 0);
params.fB = fB;
params.fD = fD;
params.fK = 0;
params.fx = 0;

% Check sizes and properties of ABCD
n = size(A,1);
params.n = n;
if size(A,2) ~= n
    error('A must be square');
end;
if size(C,2) ~= n
    error('The size of C does not correspond to the number of states');
end
if size(C,1) ~= l
    error('The size of C does not correspond to the number of outputs')
end
if ~fB
    error('This optimization does not work without a B matrix');
end
if size(B,2) ~= m
    error('The size of B does not correspond to the number of inputs')
end
if fD
    if ~all(size(D)==[l m])
        error('The size of D is incorrect');
    end
end

if  strcmp(partype,'on')
    params.partype  =  partype;
elseif  strcmp(partype,'tr')
    params.partype  =  partype;
elseif  strcmp(partype,'fl')
    params.partype  =  partype;
else
    error('You specified an unknown type of parameterization')
end

if strcmp(TT,'disc')
    if (max(abs(eig(A)))>(1-eps))
        error('Initial A must be stable.')
    end
else
    if (max(real(eig(A)))>-eps)
        error('Initial A must be stable.')
    end
end

if strcmp(TT,'disc')
    Ac = A;
    Cc = C;
    % Determine observability matrix
    ObMat = zeros(n * l,n);
    ObMat(1:l,:) = Cc;
    for blockrow = 2:n
        ObMat((blockrow-1)*n+1:blockrow*n,:) = ObMat((blockrow-2)*n+1:(blockrow-1)*n,:)*Ac;
    end
    if rank(ObMat)<n
        error('Initial model must be observable');
    end
end

% Check the optimization options: these may be either empty,
% foptions, or optimset compatible.
if ~isstruct(options) && ~isempty(options)
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
if isfield(options,'Manifold')
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
end;

% Check whether the large-scale optimization has been selected,
% and enforce RFactor='on'; in that case.
% The overparametrized (and unprojected) Tridiagonal parametrization
% does not work together with the LargeScale framework.
if strcmp(options.LargeScale,'on'),
    if isempty(options.RFactor),
        warning('LTI:resetOptions','Resetting options.RFactor from default to ''on''');
        options.RFactor = 'on';
    end;
    if ~strcmp(options.RFactor,'on'),
        error('Large-scale optimizations require OPTIONS.RFactor=''on''');
    end;
    if strcmp(params.partype,'tr')
        error('Large-scale optimizations are incompatible with the tridiagonal parametrization');
    end;
end

if strcmp(TT,'disc')
    th = dss2th(A,B,C,D,[],[],params.partype);
else
    th = css2th(A,B,C,D,[],[],params.partype);
end

[th,resnorm,residual,exitflag,output] = lmmore('ffunlti',th,[],[],options,H,w,params,options,TT);

if strcmp(TT,'disc')
    [A,B,C,D] = dth2ss(th,params);
else
    [A,B,C,D] = cth2ss(th,params);
end;

if isempty(A) && n~=0,
    warning('LTI:notStable','Optimization did not result in stable model')
    B = [];
    C = [];
    D = [];
end
