function [S,X,UM,K] = pordvarx(u,y,mu,f,p,pind,reg,opt,weight,c,ckeep)
%PORDVARX  Periodic LPV system identification using the PBSIDopt method.
%  [S,X,UM,K] = pordvarx(u,y,mu,f,p,pind) delivers information about the
%  order of the Periodic or Repeating Linear Parameter Varing system and
%  acts as a pre-processor for pmodx. The latter is used to identify an
%  open-loop or closed-loop system for the N-by-r, N-by-l and N-by-m data
%  matrices u, y and mu, where r, l and m are the number of inputs, outputs
%  and scheduling parameters. The input matrix u, output matrix y and
%  scheduling matrix mu must have the same number of observations but can
%  have different numbers of variables. The past and future window size p
%  and f must be higher then the expected system order n. The input pind
%  contains the period size and repeating index. The period size j must be
%  much larger then the number of scheduling parameters m (j>m). If window
%  size p is higher then period j, the period is multiplied until window
%  sizes p and f are lower then j.  The outputs are the singular values S,
%  which can be used to determine the order of the identifiable system.
%  Further is returned the state matrix X, wich have to be forwarded to
%  pmodx.
%
%  [S,X,UM,K] = pordvarx(u,y,mu,f,p,pind,reg,opt) adds a regularization to
%  the identification problem. The additional inputs are the regularization
%  method and selection parameters: reg = {'none', 'tikh', 'tsvd'} and opt
%  = {'gcv', 'lcurve', or any regularisation value as scalar}. With
%  regularisation, the solver can better deal with singular covariance
%  matrices. (default reg='none' and opt='gcv')
%
%  [S,X,UM,K] = pordvarx(u,y,mu,f,p,pind,reg,opt,weight) if weight=1, then
%  a left weigting matrix is added to the lowrank decomposition problem,
%  such that resulting algorithm behaves more like the "open loop" CVA
%  algorithm. (default weight=0)
%
%  [S,X,UM,K] = pordvarx(u,y,mu,f,p,pind,reg,opt,weight,c) specifies which
%  of the system matrices are constant and not parameter-varing. For each
%  of the matrices A, B, C and D an 1 or 0 can be given in the vector c. If
%  the A matrix is assumed constant the C matrix cannot be assumed constant
%  and also visa versa. If c(4) = -1 the direct feedtrough term D will not
%  be estimated. (default C=[0 0 0 0])
%
%  [S,X,UM,K] = pordvarx(u,y,mu,f,p,pind,reg,opt,weight,c,ckeep) uses a
%  L-by-M matrix CKEEP, which consists of 0(delete) or 1(keep), to reduce
%  the extended observability matrix for certain outputs. This is usefull
%  if the periodic scheduling signals P are not sufficiently exciting to
%  determine the LPV system. Currently this can only be used if A is
%  constant and C is varying. See also the wind turbine example.
%
%  References:
%    [1] van Wingerden, J.W., Houtzager, I., Verhaegen, M., Closed-loop
%    identification of the time-varying dynamics of variable-speed wind
%    turbines, Int. J. Robust Nonlinear Control 2008
%
%  See also: pmodx.m, px2abcdk.m, and px2abck.m.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check number if input arguments
if nargin < 5
    error('PORDVARX requires four or five input arguments.')
end

% assign default values to unspecified parameters
cred = true;
if (nargin < 11) || isempty(ckeep)
    cred = false;
    ckeep = 0;
end
if (nargin < 10) || isempty(c)
    c = [0 0 0 0];
end
if (nargin < 8) || isempty(weight)
    weight = 0;
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
    error('PORDVARX requires an output vector y.')
end
if ~isequal(N,length(y),length(mu))
    error('The number of rows of vectors/matrices U, Y and MU must be the same.')
end
if nargin == 8
    if (~isequal(size(ckeep,1),l) || ~isequal(size(ckeep,2),s))
        error('The matrix CKEEP must be of size L by M.');
    end
%    if c(1) ~= 1
%        error('The matrix A must be constant when using CKEEP');
%    end
    if c(3) == 1
        error('The matrix C cannot be constant when using CKEEP');
    end
end

% store the past and future vectors
m = r+l;
z = [u; y];
Z = zeros(p*m,N-p);
Y = y(:,p+1:N);
U = u(:,p+1:N);
Mu = mu(:,p+1:N);
for i = 1:p
    Z((i-1)*m+1:i*m,:) = z(:,i:N+i-p-1);
end
if c(4) ~= -1
    Z = [Z; U];
end

% check the number of repetitions (must be larger then past window)
j = 0;
for k = 1:size(pind,1)
    if length(pind{k,1}) >= p
        j = j + 1;
    end
end

% allocate matrices
X  = cell(j,1); 
UM = cell(j,1);
S  = zeros(f*l,j);
if sum(l.*s.^(1:f)) >= 1e4 || cred == true
    usparse = true;
else
    usparse = false;
end
if (c(1) == 1 && c(3) == 1) || c(1) == 2
    K = zeros(f*j*l,s*f*l);
else
    if usparse == true
        K = sparse(f*j*l,sum(l.*s.^(1:f)));
    else
        K = zeros(f*j*l,sum(l.*s.^(1:f)));
    end
end

% start for periodic identification
q = 0;
for k = 1:size(pind,1)
    nj = length(pind{k,1});
    if nj >= p
        q = q + 1;

        % solve VARX problem
        if c(4) == -1
            VARX = zeros(f*l,p*m);
        else
            VARX = zeros(f*l,p*m+r);
        end
        for i = 0:f-1
            VARX(i*l+1:(i+1)*l,:) = regress(Y(:,pind{k,1}+i),Z(:,pind{k,1}+i),reg,opt);
        end
        
        % construct LambdaKappa
        LK = zeros(f*l,p*m);
        if weight == 0
            for i = 1:f
                LK((i-1)*l+1:i*l,p*m-(p-i+1)*m+1:p*m) = VARX((i-1)*l+1:i*l,1:(p-i+1)*m);
            end
        elseif weight == 1
            for i = 0:f-1
                LK(i*l+1:(i+1)*l,i*m+1:p*m) = VARX(i*l+1:(i+1)*l,1:(p-i)*m);
                if i ~= 0
                    for j = 0:i-1
                        LK(i*l+1:(i+1)*l,:) = LK(i*l+1:(i+1)*l,:) + VARX(i*l+1:(i+1)*l,(p-i+j)*m+r+(1:l))*LK(j*l+1:(j+1)*l,:);
                    end
                end
            end
        end
        
        % singular value decomposition
        [UZ,SZ,VZ] = svd(LK*Z(1:p*m,pind{k,1}),'econ');
        S(:,q) = diag(SZ)';
        
        % build return matrices
        X{q,1} = diag(sqrt(S(1:(f*l),q)))*VZ(:,1:(f*l))';
        UM{q,1} = UZ(:,1:(f*l))*diag(sqrt(S(1:(f*l),q)));
        
        % building the matrix K
        if (c(1) == 1 && c(3) == 1) || c(1) == 2
            % for constant A
            M = zeros(f,s*f);
            for i = 1:f
                M(i,(i-1)*s+1:i*s) = Mu(:,pind{k,1}(1)+(i-1))';
            end
            K(f*l*(q-1)+1:f*l*q,:) = kron(M,eye(l));
        else
            % for varying A
            if usparse == true
                M = sparse(f,sum(s.^(1:f)));
            else
                M = zeros(f,sum(s.^(1:f)));
            end
            for i = 1:f
                d = s.^(1:f);
                if i == 1
                    if c(3) == 1
                        M(i,1) = 1;
                    else
                        M(i,1:d(1)) = Mu(:,pind{k,1}(1))';
                    end
                else
                    M(i,sum(d(1:i-1))+1:sum(d(1:i))) = kron(Mu(:,pind{k,1}(1)+(i-1))',M(i-1,sum(d(1:i-2))+1:sum(d(1:i-1))));
                end
            end
            if usparse == true
                K(f*l*(q-1)+1:f*l*q,:) = kron(M,speye(l));
            else
                K(f*l*(q-1)+1:f*l*q,:) = kron(M,eye(l));
            end
        end
    end
end

% apply reduction of extended observability matrix if requested
if cred == true
    if (c(1) == 1 && c(3) == 1) || c(1) == 2
        ckeep = kron(ones(1,f),ckeep(:)');
        ikeep = find(ckeep);
        K = K(:,ikeep);
    else
        ckeep = kron(ones(1,sum(s.^(1:f)./s)),ckeep(:)');
        ikeep = find(ckeep);
        K = K(:,ikeep);
    end
end

% apply kernel
K = full(K*K');




