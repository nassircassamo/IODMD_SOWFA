function [S,X,VARMAX] = dordvarmax(u,y,f,p,tol,reg,opt,weight)
%DORDVARMAX  Closed-loop LTI system identification using the PBSIDopt method.
%  [S,X]=dordvarmax(u,y,f,p) delivers information about the order of the
%  LTI state-space model and acts as a pre-processor for dmodx. The
%  latter is used to identify an open-loop or closed-loop system for the
%  N-by-r and N-by-l and data vectors u, y. The input matrix u and output
%  matrix y must have the same number of observations but can have
%  different numbers of variables. The past and future window size p and f
%  must be higher then the expected order n. The outputs are the singular
%  values S, which can be used to determine the order of the identifiable
%  system. Further is returned the state matrix X, wich have to be
%  forwarded to dmodx.
%
%  [S,X]=dordvarmax(u,y,f,p,tol) specifies the tolerance of the whitening
%  iterations.
%
%  [S,X]=dordvarmax(u,y,f,p,tol,reg,opt) adds a regularization to the
%  identification problem. The additional inputs are the regularization
%  method and selection parameters: reg = {'none', 'tikh', 'tsvd'} and opt
%  = {'gcv', 'lcurve'}. With regularisation, the solver can better deal
%  with singular covariance matrices. (default reg='none' and opt='gcv')
%
%  [S,X]=dordvarmax(u,y,f,p,tol,reg,opt,weight) if weight=1, then a left
%  weigting matrix is added to the lowrank decomposition problem, such that
%  resulting algorithm behaves more like the "open loop" CVA algorithm.
%  (default weight=0)
%
%  See also: dmodx.m, dx2abcdk.m, and dx2abck.m.
%
%  References:
%    [1] I. Houtzager, J. W. van Wingerden, M. Verhaegen, "VARMAX-based
%    closed-loop subspace model identification", in: 48th IEEE Conference
%    on Decision and Control, Shanghai, China, 2009.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check number if input arguments
if nargin < 4
    error('DORDVARMAX requires four or five input arguments.')
end

% assign default values to unspecified parameters
if (nargin < 8) || isempty(weight)
    weight = 0;
end
if (nargin < 7) || isempty(opt)
    reg = 'gcv';
end
if (nargin < 6) || isempty(reg)
    reg = 'none';
end
if (nargin < 5) || isempty(tol)
    tol = 1e-6;
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
        error('DORDVARMAX requires an output vector y.')
    end
    
    % store the past and future vectors
    m = r+2*l;
    z = [u; y; zeros(l,N)];
    Z = zeros(p*m,N-p);
    for i = 1:p
        Z((p-i)*m+1:(p-i+1)*m,:) = z(:,i:N+i-p-1);
    end
    
    % solve VARMAX problem
    Y = y(:,1:N-p);
    if k == 1
        [VARMAX,Z] = exls_back(Y,Z,p,r,'els',tol,reg,opt);
    else
        [VARMAX,Z] = exls_back(Y,Z,p,r,'els',tol,reg,opt,VARMAX0);
    end
    
    % calculate final VARMAX solution
%     if k == 1
%         VARMAX = regress(Y,Z,reg,opt);
%     else
%         VARMAX = regress(Y,Z,reg,opt,VARMAX0);
%     end
    
    if batch > 1
        VARMAX0 = VARMAX;
        ZZ{k} = Z;
    end
end

% construct LambdaKappa
LK = zeros(f*l,p*m);
if weight == 0
    for i = 1:f
        LK((i-1)*l+1:i*l,p*m-(p-i+1)*m+1:p*m) = VARMAX(:,1:(p-i+1)*m);
    end
elseif weight == 1
    for i = 0:f-1
        LK(i*l+1:(i+1)*l,i*m+1:p*m) = VARMAX(:,1:(p-i)*m);
        if i ~= 0
            for j = 0:i-1
                LK(i*l+1:(i+1)*l,:) = LK(i*l+1:(i+1)*l,:) + VARMAX(:,(p-i+j)*m+r+(1:l))*LK(j*l+1:(j+1)*l,:);
            end
        end
    end
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






