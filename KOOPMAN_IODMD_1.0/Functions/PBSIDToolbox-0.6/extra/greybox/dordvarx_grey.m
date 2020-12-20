function [S,X,Xd,Xs] = dordvarx(u,y,f,p,reg,opt,weight,noD)
%DORDVARX  Closed-loop LTI system identification using the PBSIDopt method.
%  [S,X]=dordvarx(u,y,f,p) delivers information about the order of the LTI
%  state-space model and acts as a pre-processor for dmodx. The latter is
%  used to identify an open-loop or closed-loop system for the N-by-r and
%  N-by-l and data vectors u, y. The input matrix u and output matrix y
%  must have the same number of observations but can have different numbers
%  of variables. The past and future window size p and f must be higher 
%  then the expected order n. The outputs are the singular values S, which 
%  can be used to determine the order of the identifiable system. Further 
%  is returned the state matrix X, wich have to be forwarded to dmodx.
%
%  [S,X]=dordvarx(u,y,f,p,reg,opt) adds a regularization to the
%  identification problem. The additional inputs are the regularization
%  method and selection parameters: reg = {'none', 'tikh', 'tsvd'} and opt
%  = {'gcv', 'lcurve', or any regularisation value as scalar}. With
%  regularisation, the solver can better deal with singular covariance
%  matrices. (default reg='none' and opt='gcv')
%
%  [S,X]=dordvarx(u,y,f,p,reg,opt,weight) if weight=1, then a left weigting
%  matrix is added to the lowrank decomposition problem, such that
%  resulting algorithm behaves more like the "open loop" CVA algorithm.
%  (default weight=0)
%
%  [S,X]=dordvarx(u,y,f,p,reg,opt,weight,noD) if noD=1, then the direct
%  feedtrough term D is not considered during estimatoion. Note the direct
%  feedtrough term D can improve the results when estimating low-order
%  models from high-order models. (default noD=0)
%
%  See also: dmodx.m, dx2abcdk.m, and dx2abck.m.
%
%  References:
%    [1] A. Chiuso, G. Picci, ``Consistency Analysis of Certain
%    Closed-Loop Subspace Identification Methods'', Automatica, Special
%    Issue on System Identification,  41(3), pp. 377--391, 2005.
%
%    [2] A. Chiuso, ``The role of vector auto regressive modeling in
%    predictor based subspace identification'', Automatica, 43(6), 
%    pp.1034–-1048, 2007.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check number if input arguments
if nargin < 4
    error('DORDVARX requires four or five input arguments.')
end

% assign default values to unspecified parameters
if (nargin < 8) || isempty(noD)
    noD = 0;
end
if (nargin < 7) || isempty(weight)
    weight = 0;
end
if (nargin < 6) || isempty(opt)
    opt = 'gcv';
end
if (nargin < 5) || isempty(reg)
    reg = 'none';
end

% check the size of the windows
if f > p
    error('Future window size f must equal or smaller then past window p. (f <= p)')
end


% check dimensions of inputs
if size(y,2) < size(y,1)
    y = y';
end
N = size(y,2);
l = size(y,1);
if isempty(u);
    r = 0;
    u = zeros(0,N);
else
    if size(u,2) < size(u,1)
        u = u';
    end
    r = size(u,1);
    if ~isequal(N,length(u))
        error('The number of rows of vectors/matrices u and y must be the same.')
    end
end
if l == 0
    error('DORDVARX requires an output vector y.')
end

% pre-allocate matrices
m = r+l;
UZ = zeros(p*r,N-p);
YZ = zeros(p*l,N-p);

% store the past and future vectors
for i = 1:N-p
    t = u(:,i:i+p-1);
    UZ(:,i) = t(:);
    t = y(:,i:i+p-1);
    YZ(:,i) = t(:);
end

% solve VARX problem
Y = y(:,p+1:N);
U = u(:,p+1:N);
Z = [UZ; YZ];
if ~noD
    Z = [Z; U];
end
if nargout > 4
    [VARX,reg_min,Zps] = regress(Y,Z,reg,opt);
else
    VARX = regress(Y,Z,reg,opt);
end

% construct LambdaKappa
LKU = zeros(f*l,p*r);
LKY = zeros(f*l,p*l);
for i = 1:f
    LKU((i-1)*l+1:i*l,p*r-(p-i+1)*r+1:p*r) = VARX(:,1:(p-i+1)*r);
    LKY((i-1)*l+1:i*l,p*l-(p-i+1)*l+1:p*l) = VARX(:,p*r+1:p*r+(p-i+1)*l);
end

% singular value decomposition
[U,S,V] = svd([LKU LKY]*Z(1:p*m,:),'econ');
X = (U*diag(sqrt(diag(S))))\([LKU LKY]*Z(1:p*m,:));
Xd = (U*diag(sqrt(diag(S))))\(LKU*UZ);
Xs = (U*diag(sqrt(diag(S))))\(LKY*YZ);
S = diag(S);