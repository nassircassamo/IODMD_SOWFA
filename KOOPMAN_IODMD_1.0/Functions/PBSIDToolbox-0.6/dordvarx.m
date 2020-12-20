function [S,X,VARX,U,Zps] = dordvarx(u,y,f,p,reg,opt,weight,noD)
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
%  If reg = {'bpdn'} and opt = {'sv', 'aic', or a scalar between 0 and 1} 
%  then sparse estimation through Basis Pursuit Denoising is used. The 
%  solver can then better deal with a past window that is chosen too 
%  large with respect to N.
%  See help private/regress for more details.
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

%  Modifications by:
%
%  Pieter Gebraad
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2011

%  Christiaan Stoeckl Torres
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2011


% check number of input arguments
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

% check for batches
if iscell(y)
    batch = length(y);
    ZZ = cell(1,batch);
    X = cell(1,batch);
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
    
    % store the past and future vectors
    m = r+l;
    z = [u; y];
    Z = zeros(p*m,N-p);
    for i = 1:p
        Z((i-1)*m+1:i*m,:) = z(:,i:N+i-p-1);
    end
    
    % solve VARX problem
    Y = y(:,p+1:N);
    U = u(:,p+1:N);
    if ~noD
        Z = [Z; U];
    end
    if k == 1
        if nargout > 4
            [VARX,reg_min,Zps] = regress(Y,Z,reg,opt);
        else
            VARX = regress(Y,Z,reg,opt);
        end
    else
        if nargout > 4
            [VARX,reg_min,Zps] = regress(Y,Z,reg,opt,VARX0);
        else
            VARX = regress(Y,Z,reg,opt,VARX0);
        end
    end
    
    if batch > 1
        VARX0 = VARX;
        ZZ{k} = Z;
    end 
end

% construct LambdaKappa
LK = zeros(f*l,p*m);
if weight == 0
    for i = 1:f
        LK((i-1)*l+1:i*l,p*m-(p-i+1)*m+1:p*m) = VARX(:,1:(p-i+1)*m);
    end
elseif weight == 1
    for i = 0:f-1
        LK(i*l+1:(i+1)*l,i*m+1:p*m) = VARX(:,1:(p-i)*m);
        if i ~= 0
            for j = 0:i-1
                LK(i*l+1:(i+1)*l,:) = LK(i*l+1:(i+1)*l,:) + VARX(:,(p-i+j)*m+r+(1:l))*LK(j*l+1:(j+1)*l,:);
            end
        end
    end
end

% singular value decomposition
if batch > 1
    QQ = zeros(f*l,f*l);
    for k = 1:batch
       QQ = QQ + LK*(ZZ{k}(1:p*m,:))*(ZZ{k}(1:p*m,:))'*LK';
    end
    [U,S] = svd(QQ,'econ');
    S = sqrt(diag(S))';
    UiLK = (U*diag(sqrt(S)))\LK;
    for k = 1:batch
        X{k} = UiLK*ZZ{k}(1:p*m,:);
    end
else
    [U,S,V] = svd(LK*Z(1:p*m,:),'econ');
    S = diag(S)';
    X = diag(sqrt(S))*V';
end

if nargout > 3
    if weight == 0
        U = diag(1./sqrt(S))*U';
    elseif weight == 1
        LK0 = zeros(f*l,p*m);
        for i = 1:f
            LK0((i-1)*l+1:i*l,p*m-(p-i+1)*m+1:p*m) = VARX(:,1:(p-i+1)*m);
        end
        U = diag(1./sqrt(S))*U'*(LK*pinv(LK0));
    end
end

end


