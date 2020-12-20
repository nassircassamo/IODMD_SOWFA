function [S,X] = lordvarx(u,y,mu,f,p,reg,opt,c,noD,ObsMatPoint)
%LORDVARX  Closed-loop LPV system identification using the PBSIDopt method.
%  [S,X]=lordvarx(u,y,mu,f,p) delivers information about the order of the
%  Linear Parameter Varing system and acts as a pre-processor for lmodx.
%  The latter is used to identify an open-loop or closed-loop system for
%  the N-by-r, N-by-l and N-by-m data matrices u, y and mu, where r, l and
%  m are the number of inputs, outputs and scheduling parameters. The input
%  matrix u, output matrix y and scheduling matrix mu must have the same
%  number of observations but can have different numbers of variables. The
%  past and future window size p and f must be higher then the expected
%  system order n. The outputs are the singular values S, which can be used
%  to determine the order of the identifiable system. Further the state 
%  matrix X is returned, which has to be forwarded to lmodx.
%
%  The data u,y,mu can be supplied in batches to prevent out-of-memory
%  errors when using a large number of samples: supply u,y,mu in a cell
%  array with different parts of the data in each cell.
%  See Example 5 ([PBSID toolbox path]/examples/ex05_lti_wts_batch.m)
%
%  [S,X]=lordvarx(u,y,mu,f,p,reg,opt) adds a regularization to the
%  identification problem. The additional inputs are the regularization
%  method and selection parameters: reg = {'none', 'tikh', 'tsvd'} and opt
%  = {'gcv', 'lcurve', or any regularisation value as scalar}. With
%  regularisation, the solver can better deal with singular covariance
%  matrices. (default reg='none' and opt='gcv')
%  If reg = {'bpdn'} and opt = {'sv', or a scalar between 0 and 1} 
%  then sparse estimation through Basis Pursuit Denoising is used. The 
%  solver can then better deal with a past window that is chosen too 
%  large.
%  See help private/regress for more details.
%
%  [S,X]=lordvarx(u,y,mu,f,p,reg,opt,c) specifies which of the system
%  matrices are constant and not parameter-varing. For each of the matrices
%  A, B, and K an 1 or 0 can be given in the vector c. 
%
%  [S,X]=lordvarx(u,y,mu,f,p,reg,opt,c,noD) if noD=1, then the direct
%  feedtrough term D is not considered during estimatoion. Note the direct
%  feedtrough term D can improve the results when estimating low-order
%  models from high-order models. (default noD=0)
%
%  [S,X]=lordvarx(u,y,mu,f,p,reg,opt,c,noD,ObsMatPoint) with ObsMatPoint=1, 
%  then in estimation of the state sequence, we consider the observability 
%  matrix in the operating point p = ones(1,m), which may yield better
%  results in some cases than with the default ObsMatPoint=0, which takes 
%  the observability matrix in the operating point p = [1,zeros(1,m-1)].
%
%  See also: lmodx.m, lx2abcdk.m, and lx2abck.m.
%
%  References:
%    [1] J.W. van Wingerden, and M. Verhaegen, ``Subspace identification
%    of Bilinear and LPV systems for open- and closed-loop data'',
%    Automatica 45, pp 372--381, 2009.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

%  Pieter Gebraad
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2011

% Jan-Willem van Wingerden
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2015


% check number if input arguments
if nargin < 4
    error('LORDVARX requires four or five input arguments.')
end

% assign default values to unspecified parameters
if (nargin < 10) || isempty(ObsMatPoint)
    ObsMatPoint = 0;
end
if (nargin < 9) || isempty(noD)
    noD = 0;
end
if (nargin < 8) || isempty(c)
    c = [0 0 0];
end
if (nargin < 7) || isempty(opt)
    opt = 'gcv';
end
if (nargin < 6) || isempty(reg)
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
if size(mu,2) < size(mu,1)
    mu = mu';
end
N = size(y,2);
l = size(y,1);
s = size(mu,1);
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
    error('LORDVARX requires an output vector y.')
end
if s == 0
    error('LORDVARX requires a scheduling sequence mu, use DORDVARX for LTI systems.')
end

% determine sizes
m = r+l;
k = r*s.^(1-c(2)+(1-c(1))*(p-1:-1:0))+ l*s.^(1-c(3)+(1-c(1))*(p-1:-1:0));
q = sum(k);
if q > (N-p)
    if ~strcmpi(reg,'bpdn')
        if ObsMatPoint == 1
            warning('lordvarx:ObsMatPoint1ThenNoKernel','Taking the observability matrix for p = ones(1,m) is not implemented for the kernel method. LORDVARX continues with ObsMatPoint=1, without kernel method.')
            kernel = 0;
        elseif ObsMatPoint == 0
            kernel = 1; 
        else 
            error('ObsMatPoint should be 0 or 1')
        end 
    else
        warning('lordvarx:BpdnThenNoKernel','The BPDN regularization is not implemented for the kernel method. LORDVARX continues with BPDN method, without kernel method.')
        kernel = 0;
    end
else
    kernel = 0;
end

% store the past and future vectors
if kernel
    Z = zeros(N-p,N-p);
    for j = 0:p-1
        Z = optkernel(Z,u,y,mu,p,c,0,j);
    end
    
else
    Z = zeros(q,N-p);
    if (c(2) == 0) && (c(3) == 0)
        z = [khatrirao(mu,u); khatrirao(mu,y)];
    elseif (c(2) == 1) && (c(3) == 0)
        z = [u; khatrirao(mu,y)];
    elseif (c(2) == 0) && (c(3) == 1)
        z = [khatrirao(mu,u); y];
    elseif (c(2) == 1) && (c(3) == 1)
        z = [u; y];
    end
    for i = 1:p
        Z(sum(k(1:i-1))+1:sum(k(1:i-1))+k(p),:) = z(:,i:N+i-p-1);
        if c(1) == 0
            for j = (i+1):p
                Z(sum(k(1:i-1))+1:sum(k(1:i-1))+k(p-j+i),:) = khatrirao(mu(:,j:N+j-p-1),Z(sum(k(1:i-1))+1:sum(k(1:i-1))+k(p-j+i+1),:));
            end
        end
    end
end

% solve VARX/KERNEL problem
if kernel
    Y = y(:,p+1:N);
    U = u(:,p+1:N);
    if ~noD
        Z = Z + U'*U;
    end 
    A = kernregress(Y,Z,reg,opt);
else
    Y = y(:,p+1:N);
    U = u(:,p+1:N);
    if ~noD
        Z = [Z; U];
    end
    VARX = regress(Y,Z,reg,opt);
end

% construct LambdaKappaZ
if kernel
    LKZ = zeros(f*l,N-p);
    for i = 0:f-1
        Z = zeros(N-p,N-p);
        for j = i:p-1
            Z = optkernel(Z,u,y,mu,p,c,i,j);
        end
        LKZ(i*l+1:(i+1)*l,:) = A*Z;
    end
    
    % singular value decomposition
    [~,S,V] = svd(LKZ,'econ');   
else
    if c(1) == 0
        if ObsMatPoint % consider the observability matrix in the operating point p = ones(1,m)
            LKZ = zeros(f*l,N-p); 
            for i = 1:f
                for j = i:p
                    for h = 1:s^(i-1)
                        LKZ((i-1)*l+1:i*l,:) = LKZ((i-1)*l+1:i*l,:) + VARX(:,sum(k(1:j-i))+((h-1)*k(j)+1:h*k(j)))*Z(sum(k(1:j-1))+1:sum(k(1:j)),:);
                    end
                end
            end
            % singular value decomposition
            [~,S,V] = svd(LKZ,'econ');
        else % consider the observability matrix in the operating point p = [1,zeros(1,m-1)]
            LK = zeros(f*l,q);        
            for i = 1:f
                for j = i:p
                    LK((i-1)*l+1:i*l,sum(k(1:j-1))+1:sum(k(1:j))) = VARX(:,sum(k(1:j-i))+1:sum(k(1:j-i))+k(j));
                end
            end
            % singular value decomposition
            [~,S,V] = svd(LK*Z(1:q,:),'econ');
        end
    else
        LK = zeros(f*l,q);
        for i = 1:f
            LK((i-1)*l+1:i*l,q-(p-i+1)*(q/p)+1:q) = VARX(:,1:(p-i+1)*(q/p));
        end
        % singular value decomposition
        [~,S,V] = svd(LK*Z(1:q,:),'econ');
    end
end
X = diag(sqrt(diag(S)))*V';
S = diag(S)';

end


function Z = optkernel(Z,u,y,mu,p,c,i,j)
N = size(y,2);
P = 1:1:N-p;
T = ones(N-p,N-p);
if all(c == 0)
    for v = 0:p-j-1
        T = T.*(mu(:,P+v+j-i)'*mu(:,P+v+j));
    end
    Z = Z + T.*([u(:,P+j-i); y(:,P+j-i)]'*[u(:,P+j); y(:,P+j)]);
else
    for v = 1:(1-c(1))*(p-j-1)
        T = T.*(mu(:,P+v+j-i)'*mu(:,P+v+j));
    end
    if c(2)
        Z = Z + T.*(u(:,P+j-i)'*u(:,P+j));
    else
        Z = Z + T.*(mu(:,P+j-i)'*mu(:,P+j)).*(u(:,P+j-i)'*u(:,P+j));
    end
    if c(3)
        Z = Z + T.*(y(:,P+j-i)'*y(:,P+j));
    else
        Z = Z + T.*(mu(:,P+j-i)'*mu(:,P+j)).*(y(:,P+j-i)'*y(:,P+j));
    end
end
end




