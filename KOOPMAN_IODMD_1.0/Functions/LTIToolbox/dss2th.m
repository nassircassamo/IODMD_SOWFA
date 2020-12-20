function [theta,params,T] = dss2th(varargin)
%DSS2TH     This function converts a discrete time state space
%           model to a parameter vector that describes the model
%           Model structure:
%                x(k+1) = Ax(k) + Bu(k) + Ke(k)
%                y(k)   = Cx(k) + Du(k) + e(k)
% Syntax:
%           [theta,params,T] = dss2th(A,C,partype);
%           [theta,params,T] = dss2th(A,B,C,partype);
%           [theta,params,T] = dss2th(A,B,C,D,partype);
%           [theta,params,T] = dss2th(A,B,C,D,x0,partype);
%           [theta,params,T] = dss2th(A,B,C,D,x0,K,partype);
%
% Input:
% A,B,C,D   System matrices describing the state space system.
%           The B and D matrices are optional and can be left out
%           or given as an empty matrix to indicate it is not part
%           of the model.
% x0        Initial state, This is optional.
% K         Kalman gain. Also this matrix is optional.
%
% Rules for input parameters:
%
%           The final parameter should always be the parametrization
%           type. The order for the parameters prior to partype
%           is A,B,C,D,x0,K. The only exception is A,C, when only
%           those are to be parametrized. This kludge exists to
%           prodive compatibility with the dslslin function.
%
%           All parameters after A,B,C and before partype are optional.
%           If the last one is not to be parametrized it can be omitted.
%           If any other is not to be parametrized, an empty matrix
%           should be passed.
%
%           (A,B,C,partype) thus is equivalent to (A,B,C,[],[],[],partype)
%           However, (A,B,C,[],x0,partype) cannot be abbreviated.
%
% partype   This parameter specifies the type of parameterization
%           that is used to parameterize the state space model.
%           Three types of parameterization are supported:
%           'on'= Output Normal, 'tr'=TRidiagonal and
%           'fl'= FuLl.
%
% Output:
% theta     Parameters vector describing the system.
% params    A structure that contains the dimension parameters of
%           the system, such as the order, the number of inputs
%           whether D, x0 or K is present, etc.
% T         Transformation matrix between the input state space system
%           and the state space system in the form described by theta.
%
% See Also : dth2ss

% This function is loosely based on dss2th.m in the SMI-2.0
% toolbox. Johan Bruls 1996, Bert Haverkamp 2000
%
% Revised by Niek Bergboer, 2001
% Revised by Ivo Houtzager, 2007
% Copyright (c) 1996-2007, Delft Center of Systems and Control

% Issue help if no parameters are given
if nargin < 3
    error('DSS2TH requires more then three input arguments.')
end

% Set binary mask for matrices to be parametrized
partype = varargin{nargin};
fB = 0;
fD = 0;
fx = 0;
fK = 0;
if nargin==3,
    A  =  varargin{1};
    C  =  varargin{2};
    m = 0;
    
    % Check input parameters
    [ma,na] = size(A);
    nc = size(C,2);
    if ma ~= na,
        error('A matrix must be square.')
    end
    if na ~= nc,
        error('C matrix must have as many columns as A.')
    end
else
    % A,B,C are given, always parametrize them
    fB = 1;
    A = varargin{1};
    B = varargin{2};
    C = varargin{3};
    m = size(B,2);
    
    % Check input parameters
    [ma,na] = size(A);
    [mb,nb] = size(B);
    [mc,nc] = size(C);
    if ma ~= na
        error('A matrix must be square.')
    end
    if mb ~= ma
        error('B matrix must have as many rows as A.')
    end
    if na ~= nc
        error('C matrix must have as many columns as A.')
    end

    if nargin >= 5
        % Four matrices given
        D = varargin{4};
        fD = ~isempty(D);
        if fD == 1
            % Check input parameters
            [md,nd] = size(D);
            if md ~= mc
                error('D matrix must have as many rows as C.')
            end
            if nb ~= nd
                error('D matrix must have as many columns as B.')
            end
        end
    end
    if nargin >= 6
        % Four matrices and initial state given
        x0 = varargin{5};
        fx = ~isempty(x0);
        if fx == 1
            % Check input parameters
            [mx,nx] = size(x0);
            if mx ~= ma
                error('X0 vector must have as many rows as A.')
            end
            if nx ~= 1
                error('X0 vector must consists of 1 column.')
            end
        end
    end
    if nargin == 7
        % Four matrices, initial state and Kalman gain given
        K = varargin{6};
        fK = ~isempty(K);
        if fK == 1
            % Check input parameters
            [mk,nk] = size(K);
            if mk ~= ma
                error('K matrix must have as many rows as A.')
            end
            if nk ~= mc
                error('The number of columns of K must be equal to the number of rows of C.')
            end
        end
    end
    if nargin > 7
        error('Too many input arguments');
    end
end

[l,n] = size(C);
if length(partype)>2
    partype = partype(1:2);
end
if strcmp(partype,'on')
%    nn = n*l;
elseif strcmp(partype,'tr')
%    nn = n*l + 3*n - 2;
elseif strcmp(partype,'fl')
%    nn = n*(n + l);
else
    error('You specified an unknown type of parameterization')
end
%nl = fB*m*n + fD*m*l + fx*n + fK*n*l;

params = struct('n',n,'m',m,'l',l,'fB',fB,'fD',fD,'fx',fx,'fK',fK,'partype',partype);
if ~fB
    B = [];
end
if ~fK
    K = [];
end
if ~fx
    x0 = [];
end
if ~fD
    D = [];
end

if strcmp(partype,'on') && any(max(abs(eig(A)))>=(1-eps))
    error('A Pole was found on the imaginary axis. The output normal form is not possible for unstable or marginally stable systems.');
end

if ~strcmp(params.partype,'fl')
    [thn,T] = dac2thn(A,C,params);
    if nargin > 3
        if params.fB
            B = T\B;
        end
        if params.fK
            K = T\K;
        end
        if params.fx
            x0 = T\x0;
        end
        thl = bdxk2thl(B,D,x0,K,params);
    else
        thl = [];
    end
    if n ~= 0
        theta = [thn; thl];
    else
        % Special case for static systems
        theta = thl;
    end
else
    if fD
        M = [A B; C D];
        theta = M(:);
    else
        M = [A; C];
        theta = [M(:); B(:)];
    end

    % Append K below, if desired
    if fK
        theta = [theta; K(:)];
    end

    % And x0, if selected
    if fx
        theta = [theta;x0];
    end
    T = eye(n);
end
end

function [thn,T] = dac2thn(A,C,params)
partype = params.partype;
n = params.n;
l = params.l;
Il = eye(l);
if strcmp(partype,'on')
    % output normal
    [A,C,T] = dac2on(A,C,params);
    thn = zeros(n*l,1);
    Z = [C; A];
    for i = 1:n
        si = Z(n-i+2:n+l-i+1,n-i+1);
        ri = Z(n-i+1,n-i+1);
        ti = si'*si;
        Ti = eye(n+l,n+l);
        if ti > eps
            Ti(n-i+1:n+l-i+1,n-i+1:n+l-i+1) = [-si Il-(1-ri)/ti*(si*si'); ri si'];
        else
            Ti(n-i+1:n+l-i+1,n-i+1:n+l-i+1) = [-si Il; ri transpose(si)];
        end
        Z = Ti*Z;
        thn(l*(i-1)+1:l*i) = si;
    end
elseif strcmp(partype,'tr')
    % Tridiagonal
    [A,C,T] = ac2tr(A,C);
    if n > 1
        thn = [diag(A,1); diag(A); diag(A,-1); C(:)];
    else
        % First-order and static systems.
        thn = [A; C];
    end
end
end

function thl = bdxk2thl(B,D,x0,K,params)
fD = params.fD;
fx = params.fx;
fK = params.fK;
thl = B(:);
if fD
    thl = [thl; D(:)];
end
if fx
    thl = [thl; x0];
end
if fK
    thl = [thl; K(:)];
end
end

function [A,C,T] = dac2on(A,C,params)
n = params.n;
%l = params.l;

% Make system output normal
Wo = dgram(A',C');
Q1 = chol(Wo);
A = Q1*A/Q1;
C = C/Q1;
%Z = [C; A];

% Make observability matrix positive lower triangular
Gam = obsv(A,C);
[Q2,R] = qr(Gam');

%Gam*Q2
Q3 = diag(sign(diag(R(1:n,1:n))));
Q4 = Q2*Q3;
A = Q4\A*Q4;
C = C*Q4;
T = pinv(Q1)*Q4;
end

function [A,C,T] = ac2tr(A,C)
[V,Lambda] = eig(A);
Lambda = diag(Lambda);

% First, separate the real and complex eigenvalues.
n = length(Lambda);
flag = zeros(n,1);
for kl = 1:n
    flag(kl) = isreal(Lambda(kl));
end
RealIndices = find(flag);
ComplexIndices = find(1-flag);
nc = length(ComplexIndices);
nr = length(RealIndices);
LambdaReal = Lambda(RealIndices);
LambdaComplex = Lambda(ComplexIndices);
Lambda = [LambdaComplex; LambdaReal];
VReal = V(:,RealIndices);
VComplex = V(:,ComplexIndices);
V = [VComplex VReal];

% Now, convert pairs of eigenvalues into 2by2 blocks -- except for the last
% (real) eigenvalue, when n is odd.
T = eye(n);
for kl = 1:2:nc
    lambda = Lambda(kl);
    a = real(lambda);
    b = imag(lambda);
    Mk = [a b;-b a];
    [Vk,Dk] = eig(Mk);
    Tk = inv(Vk);
    T(kl:kl+1,kl:kl+1) = Tk;
end
for kl = nc+1:2:nc+nr-1
    lambda1 = Lambda(kl);
    lambda2 = Lambda(kl+1);
    a = (lambda1 + lambda2)/2;
    b = (lambda1 - lambda2)/2;
    Mk = [a b; b a];
    [Vk,Dk] = eig(Mk);
    Tk = inv(Vk);
    T(kl:kl+1,kl:kl+1) = Tk;
end
T = V*T;
A = T\A*T;
A = tril(triu(A,-1),1);
C = C*T;
end

