function [S,X] = lordvarxydist(u,d,y,mu,f,p,reg,opt,c,noD,ObsMatPoint)
%LORDVARXYDIST is LORDVARX for a special case with output disturbances
%   x(k+1) = A kron(mu(k),x(k)) + B kron(mu(k),u(k)) + K kron(mu(k),e(k))
%   y(k)   = C x(k) + [Du, Dd] [u(k); kron(mu(k),d(k))] + e(k)
%   
%  if c(4)=1 then Dd is not varying with the scheduling mu
%  
%  See also: lordvarx.m and lx2abcdkydist.m.
%
%  References:
%    [1] J.W. van Wingerden, and M. Verhaegen, ``Subspace identification
%    of Bilinear and LPV systems for open- and closed-loop data'',
%    Automatica 45, pp 372--381, 2009.

%  Pieter Gebraad
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2011


% check number if input arguments
if nargin < 6
    error('LORDVARX requires four or five input arguments.')
end

% assign default values to unspecified parameters
if (nargin < 11) || isempty(ObsMatPoint)
    ObsMatPoint = 0;
end
if (nargin < 10) || isempty(noD)
    noD = 0;
end
if (nargin < 9) || isempty(c)
    c = [0 0 0 0];
end
if (nargin < 8) || isempty(opt)
    opt = 'gcv';
end
if (nargin < 7) || isempty(reg)
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
if size(d,2) < size(d,1)
    d = d';
end
N = size(y,2);
l = size(y,1);
rd = size(d,1);

if ~isequal(N,length(d))
    error('The number of rows of vectors/matrices d and y must be the same.')
end

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

if c(4)==0
    d = khatrirao(mu,d);
end

% determine sizes
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
        if c(1) ~= 0
            for j = (i+1):p
                Z(sum(k(1:i-1))+1:sum(k(1:i-1))+k(p-j+i),:) = khatrirao(mu(:,j:N+j-p-1),Z(sum(k(1:i-1))+1:sum(k(1:i-1))+k(p-j+i+1),:));
            end
        end
    end
end

Y = y(:,p+1:N);
U = u(:,p+1:N);
d = d(:,p+1:N);

% solve VARX/KERNEL problem
if kernel
    if ~noD
        Z = Z + U'*U + d'*d;
    end 
    A = kernregress(Y,Z,reg,opt);
else

    if ~noD
        Z = [Z; U; d];
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




