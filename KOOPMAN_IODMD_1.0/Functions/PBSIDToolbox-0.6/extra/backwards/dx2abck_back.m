function [A,B,C,K] = dx2abck(x,u,y,f,p,c)
%DX2ABCK Estimates the matrices A, B, and C of the state space model
%  [A,B,C] = dx2abck(x,u,y,f,p) estimates the matrices A, B, and C of the
%  state space model:
%
%       x(k+1) = A x(k) + B u(k) + K e(k)
%       y(k)   = C x(k) + e(k)
%
%  using the knowledge of the state vector x, the input vector u and the
%  output vector u. The past window size p is recomended to be higher then
%  the expected system order n. Future window size f must equal or smaller
%  then past window size p.
%
%  [A,B,C] = dx2abck(x,u,y,f,p,'stable') estimates a stable matrix A by
%  using the method in [1].
%
%  [A,B,C,K] = dx2abck(x,u,y,f,p) returns a Kalman matrix K calculated by
%  the discrete Riccati equation, which gives a guaranteed stable A-KC
%  (predictor from).
%
%  See also: dmodx.m, dordvarx.m, dordvarmax.m.
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
    if size(x,2) < size(x,1)
        x = x';
    end
    N = size(y,2);
    l = size(y,1);
    n = size(x,1);
    if l == 0
        error('DX2ABC requires an output vector y.')
    end
    if n == 0
        error('DX2ABC requires an state vector x.')
    end
    if ~isequal(N,length(u))
        error('The number of rows of vectors/matrices u and y must be the same.')
    end
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
    e = y - C*x;
    
    % calculate the A and B matrices
    z = vertcat(x(:,1:end-1),u(:,1:end-1),e(:,1:end-1));
    if k == 1
        ABK = x(:,2:end)*pinv(z);
    else
        ABK = [A0 B0 K0] + (x(:,2:end)-[A0 B0 K0]*z)*pinv(z);
    end
    A = ABK(:,1:n);
    B = ABK(:,n+1:n+r);
    %K = ABK(:,n+r+1:n+r+l);
    
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
        z = vertcat(u(:,1:end-1),e(:,1:end-1));
        if k == 1
            BK = (x(:,2:end) - A*x(:,1:end-1))*pinv(z);
        else
            BK = [B0 K0] + (x(:,2:end) - A*x(:,1:end-1) - [B0 K0]*z)*pinv(z);
        end
        B = BK(:,1:r);
    end
    
    % If selected, find a quaranteed stable A-KC matrix
    if nargout > 3
        % estimate covariance matrices
        VW = [x(:,2:end); y(:,1:end-1)] - [A B; C zeros(l,r)]*[x(:,1:end-1); u(:,1:end-1)];
        if k == 1
            QSR = (VW*VW')./length(VW);
        else
            VW = [VW VW0];
            QSR = (VW*VW')./length(VW);
        end
        Q = QSR(1:n,1:n);
        S = QSR(1:n,n+1:n+l);
        R = QSR(n+1:n+l,n+1:n+l);
        
        % calculate Kalman gain with discrete Riccati equation
        [~,~,K] = dare(A',C',Q,R,S);
        K = K';
    end
    
    if batch > 1
        A0 = A;
        B0 = B;
        C0 = C;
        K0 = K;
        VW0 = VW;
    end
end