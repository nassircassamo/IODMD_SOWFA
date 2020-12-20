function [A,B,C,D,x0,K] = cth2ss(theta,params,T) 
%CTH2SS     This function converts a parameter vector 
%           that describes a continuous time state space model 
%           in output normal form to the state space matrices 
%           of that model. 
%           Model structure 
%              . 
%              x(t) = Ax(t) + Bu(t)+ K e(t) 
%              y(t) = Cx(t) + Du(t) + e(t) 
%              x(0) = x0 
% 
% Syntax: 
%           [A,C] = cth2ss(theta,params) 
%           [A,B,C] = cth2ss(theta,params) 
%           [A,B,C,D] = cth2ss(theta,params) 
%           [A,B,C,D,x0] = cth2ss(theta,params) 
%           [A,B,C,D,x0,K] = cth2ss(theta,params) 
% Input: 
%  theta    Parameter vector describing the system. 
%  params   A structure that contains the dimension parameters of 
%           the system, such as the order, the number of inputs, 
%           whether D, x0 or K are present, etc. 
%  T        Transformation matrix to be applied to the state space 
%           system that is constructed from theta. This transformation 
%           might come from the function dss2th. 
% 
% Output: 
%   A,B,C,D System matrices describing the state-space 
%           system in output normal form. If theta does 
%           not contain  parameters for D, this matrix 
%           will be returned as an empty matrix. 
%   x0      Initial condition. If theta does not contain 
%           parameters for x0, this vector will be returned 
%           as an empty matrix. 
%   K       Kalman gain. If not parametrized, it will be 
%           returned empty. 
% 
% See also css2th 
 
% This function is closely based on cth2ss.m from the SMI-2.0 
% toolbox by Johan Bruls 1996, Bert Haverkamp 2000 
% 
% Revised by Niek Bergboer, 2001 
% Revised by Ivo Houtzager, 2007
% Copyright (c) 1996-2007, Delft Center of Systems and Control 

% Issue help if no parameters are given
if nargin < 1
    error('CTH2SS requires more then one input arguments.')
end

% Assign default values to unspecified parameters
if nargin < 3
    T  =  [];
end

if isstruct(params)
    n = params.n;
    l = params.l;
    m = params.m;
    if nargout==2
        params.fB = 0;
        params.fD = 0;
        params.fx = 0;
        params.fK = 0;
    end
    fB = params.fB;
    fD = params.fD;
    fx = params.fx;
    fK = params.fK;
    partype = params.partype;
end

if strcmp(partype,'on');
    nn = n*l;
elseif strcmp(partype,'tr');
    if n ~= 0
        nn = n*l + 3*n - 2;
    else
        % Special case for static systems
        nn  =  0;
    end
elseif strcmp(partype,'fl');
    nn = n*(n + l);
else
    error('You specified an unknown type of parameterization in params.partype')
end
nt = length(theta);
nl = fx*n + fB*n*m + fD*l*m + fK*n*l;

if nt ~= nn + nl
    error('The length of theta and the values in params do not correspond.');
end

if strcmp(params.partype,'fl')
    thn = theta(1:nn);
    thl = theta(nn+1:nt);

    if n ~= 0
        [A,C] = cthn2ac(thn,params);
    else
        A = [];
        C = zeros(l,0);
    end
    if nargout <= 2
        B = C;
    else
        [B,D,x0,K] = thl2bdxk(thl,params);
    end
    if ~isempty(T)
        A = T*A/T;
        C = C/T;
        if fB
            B = T*B;
        end
        if fx
            x0 = T*x0;
        end
        if fK
            K = T*K;
        end
    end
else
    if fD
        M = zeros(n+l,n+m);
        M(:) = theta(1:(n+l)*(n+m));
        A = M(1:n,1:n);
        B = M(1:n,n+1:n+m);
        C = M(n+1:n+l,1:n);
        D = M(n+1:n+l,n+1:n+m);

        % thpos indicates from which position x0 can be extracted (if wanted)
        thpos = (n+l)*(n+m);
    else
        M = zeros(n+l,n);
        M(:) = theta(1:(n+l)*n);
        A = M(1:n,1:n);
        C = M(n+1:n+l,1:n);
        B = zeros(n,m);
        B(:) = theta((n+l)*n+1:(n+l)*n + n*m);
        thpos = (n+l)*n + n*m;
        D = [];
        if ~fB
            B = [];
        end
    end
    if fK
        K = zeros(n,l);
        K(:) = theta(thpos+1:thpos+n*l);
        thpos = thpos + n*l;
    else
        K = [];
    end
    if fx,
        x0 = theta(thpos+1:thpos+n);
    else
        x0 = [];
    end
end
end

function [A,C] = cthn2ac(thn,params)
partype = params.partype;
l = params.l;
n = params.n;
if strcmp(partype,'on');
    C = zeros(l,n);
    offset = 0;
    for j = 1:l
        C(j:l,j) = thn(offset+1:offset+l-j+1);
        offset = offset+l-j+1;
    end
    Ass = zeros(n);
    for j = 1:min(l,n-1)
        Ass = Ass+diag(thn(offset+1:offset+n-j),j)-diag(thn(offset+1:offset+n-j),-j);
        offset = offset+n-j;
    end
    A = -0.5*(C'*C)+Ass;
elseif strcmp(partype,'tr');
    % tri-diagonal
    A = zeros(n);
    C = zeros(l,n);
    for j = 1:n-1
        A(j,j+1) = thn(j);
        A(j,j) = thn(n-1+j);
        A(j+1,j) = thn(2*n-1+j);
    end
    A(n,n) = thn(2  *  n-1);
    C(:) = thn(3*n-1:3*n-2+n*l);
end
end

function [B,D,x0,K] = thl2bdxk(thl,params)
fB = params.fB;
fx = params.fx;
fD = params.fD;
fK = params.fK;
n = params.n;
m = params.m;
l = params.l;

thl_pos = 0;
if fB
    B = zeros(n,m);
    B(:) = thl(thl_pos+1:thl_pos+m*n);
    thl_pos = thl_pos + m*n;
else
    B = [];
end
if fD
    D = zeros(l,m);
    D(:) = thl(thl_pos+1:thl_pos+l*m);
    thl_pos = thl_pos + l*m;
else
    D = [];
end
if fx
    x0 = thl(thl_pos+1:thl_pos+n);
    thl_pos = thl_pos + n;
else
    x0 = [];
end
if fK
    K = zeros(n,l);
    K(:) = thl(thl_pos+1:thl_pos+n*l);
else
    K = [];
end
end