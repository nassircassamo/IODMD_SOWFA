function [S,X,FIR] = dordfir(u,y,f,p,reg,opt,noD)
%DORDFIR  Open-loop LTI system identification using the PBSIDopt method.
%  [S,X]=dordfir(u,y,f,p) delivers information about the order of the LTI
%  state-space model and acts as a pre-processor for dmodx. The latter is
%  used to identify an open-loop or closed-loop system for the N-by-r and
%  N-by-l and data vectors u, y. The input matrix u and output matrix y
%  must have the same number of observations but can have different numbers
%  of variables. The past and future window size p and f must be higher 
%  then the expected order n. The outputs are the singular values S, which 
%  can be used to determine the order of the identifiable system. Further 
%  is returned the state matrix X, wich have to be forwarded to dmodx.
%
%  [S,X]=dordfir(u,y,f,p,reg,opt) adds a regularization to the
%  identification problem. The additional inputs are the regularization
%  method and selection parameters: reg = {'none', 'tikh', 'tsvd'} and opt
%  = {'gcv', 'lcurve', or any regularisation value as scalar}.}. With
%  regularisation, the solver can better deal with singular covariance
%  matrices. (default reg='none' and opt='gcv')
%
%  See also: dmodx.m, dx2abcd.m, and dx2abc.m.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check number if input arguments
if nargin < 4
    error('DORDFIR requires four or five input arguments.')
end

% assign default values to unspecified parameters
if (nargin < 7) || isempty(noD)
    noD = 0;
end
if (nargin < 6) || isempty(opt)
    opt = 'gcv';
end
if (nargin < 5) || isempty(reg)
    reg = 'none';
end

% check for batches
if iscell(y)
    batch = length(y);
    ZZ = cell(1,batch);
    XX = cell(1,batch);
    yb = y;
    ub = u;
else
    batch = 1;
end

% do for all batches
for k = 1:batch
    if batch > 1
        y = yb{k};
        u = ub{k};
    end
    
    % check dimensions of inputs
    if size(y,2) < size(y,1)
        y = y';
    end
    if size(u,2) < size(u,1)
        u = u';
    end
    N = size(y,2);
    l = size(y,1);
    r = size(u,1);
    if r == 0
        error('DORDFIR requires an input vector u.')
    end
    if l == 0
        error('DORDFIR requires an output vector y.')
    end
    if ~isequal(N,length(u))
        error('The number of rows of vectors/matrices u and y must be the same.')
    end
    
    % store the past and future vectors
    m = r;
    Z = zeros(p*m,N-p);
    for i = 1:p
        Z((i-1)*m+1:i*m,:) = u(:,i:N+i-p-1);
    end
    
    % solve VARX problem
    Y = y(:,p+1:N);
    U = u(:,p+1:N);
    if ~noD
        Z = [Z; U];
    end
    if k == 1
        FIR = regress(Y,Z,reg,opt);
    else
        FIR = regress(Y,Z,reg,opt,FIR0);
    end
    
    if batch > 1
        FIR0 = FIR;
        ZZ{k} = Z;
    end
end

% construct LambdaKappa
LK = zeros(f*l,p*m);
for i = 1:f
    LK((i-1)*l+1:i*l,p*m-(p-i+1)*m+1:p*m) = FIR(:,1:(p-i+1)*m);
end

% singular value decomposition
if batch > 1
    Z = [ZZ{:}];
end
[U,S,V] = svd(LK*Z(1:p*m,:),'econ');
X = diag(sqrt(diag(S)))*V';
S = diag(S)';    
% do for all batches
for k = 1:batch
    if batch > 1
        XX{k} = X(:,1:size(ZZ{k},2));
        X = X(:,size(ZZ{k},2)+1:end); 
    end
end

% return batches
if batch > 1
    X = XX;
end





