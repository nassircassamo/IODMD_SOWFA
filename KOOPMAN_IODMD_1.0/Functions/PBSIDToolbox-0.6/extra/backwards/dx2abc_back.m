function [A,B,C] = dx2abc(x,u,y,f,p,c)
%DX2ABC Estimates the matrices A, B, and C of the state space model
%  [A,B,C]=dx2abc(x,u,y,f,p) estimates the matrices A, B, and C of the
%  state space model:
%
%       x(k+1) = A x(k) + B u(k)
%       y(k)   = C x(k) + e(k)
%
%  using the knowledge of the state vector x, the input vector u and the
%  output vector u. The past window size p is recomended to be higher then
%  the expected system order n. Future window size f must equal or smaller
%  then past window size p.
%
%  [A,B,C]=dx2abcd(x,u,y,f,p,'stable') estimates a stable matrix A by using
%  the method in [1].
%
%  See also: dmodx.m, dordfir.m.
%
%  References:
%    [1] J.M. Maciejowski, "Guaranteed Stability with Subspace Methods",
%    Submitted to Systems and Control Letters, 1994.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check number if input arguments
if nargin == 5 || isempty(c)
    c = 'none';
end
if nargin < 5
    error('DX2ABC requires at least five input arguments.')
end

% check for batches
if iscell(y)
    batch = length(y);
    yb = y;
    ub = u;
    xb = x;
else
    batch = 1;
end

% do for all batches
for k = 1:batch
    if batch > 1
        y = yb{k};
        u = ub{k};
        x = xb{k};
    end
    
    % check dimensions of inputs
    if size(y,2) < size(y,1)
        y = y';
    end
    if size(u,2) < size(u,1)
        u = u';
    end
    if size(x,2) < size(x,1)
        x = x';
    end
    N = size(y,2);
    l = size(y,1);
    r = size(u,1);
    n = size(x,1);
    if r == 0
        error('DX2ABC requires an input vector u.')
    end
    if l == 0
        error('DX2ABC requires an output vector y.')
    end
    if n == 0
        error('DX2ABC requires an state vector x.')
    end
    if ~isequal(N,length(u))
        error('The number of rows of vectors/matrices u and y must be the same.')
    end
    
    % check if the input signal is sufficiently exciting
    if rank(u) < r
        warning('CLID:ranku','The input vector u is not sufficiently exciting. (rank(u) = r)')
    end
    
    % check if the state vector is full rank
    if rank(x) < n
        error('The state vector x is not full rank. (rank(x) = n)')
    end
    
    % check the size of the windows
    if f > p
        error('Future window size f must equal or smaller then past window p. (f <= p)')
    end
    
    % remove the window sizes from input and output vector
    u = u(:,1:size(x,2));
    y = y(:,1:size(x,2));
    
    % calculate the C matrix
    if k == 1
        C = y*pinv(x);
    else
        C = C0 + (y-C0*x)*pinv(x);
    end
    
    % calculate the A and B matrices
    z = vertcat(x(:,1:end-1),u(:,1:end-1));
    if k == 1
        AB = x(:,2:end)*pinv(z);
    else
        AB = [A0 B0] + (x(:,2:end)-[A0 B0]*z)*pinv(z);
    end
    A = AB(:,1:n);
    B = AB(:,n+1:n+r);
    
    % If selected, find a quaranteed stable A matrix
    if nargin > 5 && strcmpi(c,'stable')
        Gamma = zeros((p+1)*l,n);
        for i = 1:p
            if i == 1
                Gamma((i-1)*l+1:i*l,:) = C;
            else
                Gamma((i-1)*l+1:i*l,:) = Gamma((i-2)*l+1:(i-1)*l,:)*A;
            end
        end
        A = pinv(Gamma(1:p*l,:))*Gamma(l+1:(p+1)*l,:);
        
        % recalculate with stable A matrix
        z = vertcat(u(:,1:end-1));
        if k == 1
            B = (x(:,2:end) - A*x(:,1:end-1))*pinv(z);
        else
            B = B0 + (x(:,2:end) - A*x(:,1:end-1) - B0*z)*pinv(z);
        end
    end
    
    if batch > 1
        A0 = A;
        B0 = B;
        C0 = C;
    end
end