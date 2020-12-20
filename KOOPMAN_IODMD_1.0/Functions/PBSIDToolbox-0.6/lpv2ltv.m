function [Altv,Bltv,Cltv,Dltv,Kltv] = lpv2ltv(varargin)
%LPV2LTV  Transforms a LPV system matrices to a LTV system matrices
%  [Altv,Bltv,Cltv] = lpv2ltv(Alpv,Blpv,Clpv,mu) estimates the linear
%  time-varying matrices A, B, and C from the periodic LPV state space
%  model: The inputs are the linear parameter-varying matrices A, B, and C,
%  where matrices have the form of A=[A(1) A(2) ... A(m)]. The outputs are
%  the linear parameter-varying matrices A, B, and C have the form of A{1},
%  A{2}, ..., A{N}.
%
%  [Altv,Bltv,Cltv,Dltv] = lpv2ltv(Alpv,Blpv,Clpv,Dlpv,mu) estimates the
%  linear time-varying matrices A, B, C, and D from the periodic LPV state
%  space model: The inputs are the linear parameter-varying matrices A, B,
%  C and D, where matrices have the form of A=[A(1) A(2) ... A(m)]. The
%  outputs are the linear parameter-varying matrices A, B, C, and D have
%  the form of A{1}, A{2}, ..., A{N}.
%
%  [Altv,Bltv,Cltv,Dltv,Kltv] = lpv2ltv(Alpv,Blpv,Clpv,Dlpv,Klpv,mu) returns
%  also the estimated Kalman matrix K. Note that this Kalman gain does not
%  give a guaranteed stable A-KC (predictor from).
%
%  [Altv,Bltv,Cltv,Dltv,Kltv] = lpv2ltv(Alpv,Blpv,Clpv,Dlpv,x,u,y,mu)
%  returns a Kalman matrix K calculated by the discrete-time non-steady
%  Riccati equation, which gives a guaranteed stable A-KC (predictor from).

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check input arguments
if nargin < 5
    error('LPV2LTV requires at least five input arguments!');
end

% assign values to unspecified parameters from command
Alpv = varargin{1};
Blpv = varargin{2};
Clpv = varargin{3};
mu = varargin{nargin};
N = length(mu);
Altv = cell(N,1);
Bltv = cell(N,1);
Cltv = cell(N,1);
Dltv = cell(N,1);
Kltv = cell(N,1);
switch nargin
    case 5
        Dlpv = varargin{4}; 
    case 6
        Dlpv = varargin{4};
        Klpv = varargin{5};
    case 8
        Dlpv = varargin{4};
        x     = varargin{5};
        u     = varargin{6};
        y     = varargin{7};
    otherwise
        error('Command not recognized.')
end
  
% determine lpv system sizes
if size(mu,1) < size(mu,2)
    mu = mu';
end
n = size(Alpv,1);          % The order of the system
m = size(Alpv,2)/n - 1;    % The number of scheduling parameters
r = size(Blpv,2)/(m + 1);  % The number of inputs
l = size(Clpv,1);          % The number of outputs

% build ltv cells
for i = 1:N
    Altv{i} = Alpv(1:n,1:n);
    Bltv{i} = Blpv(1:n,1:r);
    Cltv{i} = Clpv(1:l,1:n);
    if nargout > 3
        Dltv{i} = Dlpv(1:l,1:r);
    end
    if nargout > 4
        Kltv{i} = Klpv(1:n,1:l);
    end
    for j = 1:m
        Altv{i} = Altv{i} + Alpv(1:n,(j-1)*n+1+n:j*n+n).*mu(i,j);
        Bltv{i} = Bltv{i} + Blpv(1:n,(j-1)*r+1+r:j*r+r).*mu(i,j);
        Cltv{i} = Cltv{i} + Clpv(1:l,(j-1)*n+1+n:j*n+n).*mu(i,j);
        if nargout > 3
            Dltv{i} = Dltv{i} + Dlpv(1:l,(j-1)*r+1+r:j*r+r).*mu(i,j);
        end
        if nargout > 4
            Kltv{i} = Kltv{i} + Klpv(1:n,(j-1)*l+1+l:j*l+l).*mu(i,j);
        end
    end
end

if nargin == 9
    % check signal vector sizes
    if size(x,1) < size(x,2)
        x = x';
    end
    if size(u,1) < size(u,2)
        u = u';
    end
    if size(y,1) < size(y,2)
        y = y';
    end
    nn = floor(size(y,1)/p);
    
    % build periodic Kalman gain
    Ahat = zeros(n,n,p);
    Chat = zeros(n,l,p);
    Q = zeros(n,n,p);
    S = zeros(n,l,p);
    R = zeros(l,l,p);
    for i = 1:N
        % caculate residuals
        Ahat(:,:,i) = Altv{i}';
        Chat(:,:,i) = Cltv{i}';
        VW = [x(i+1:p:p*(nn-1),:)'; y(i:p:p*(nn-1),:)'] -...
            [Ahat(:,:,i)' Bltv{i}; Chat(:,:,i)' Dltv{i}]*...
            [x(i:p:p*(nn-1),:)'; u(i:p:p*(nn-1),:)'];
        
        % estimate covariance matrices
        QSR = (VW*VW')./(nn-1);
        Q(:,:,i) = QSR(1:n,1:n);
        S(:,:,i) = QSR(1:n,n+1:n+l);
        R(:,:,i) = QSR(n+1:n+l,n+1:n+l);
    end
    
    % calculate periodic Kalman gain with discrete Riccati equation
    F = zeros(n,n);
    Gp = dpkalm(Ahat,Chat,Q,R,S,F,N);
    for i = 1:N
        Kltv{i} = Gp(:,:,i)';
    end
end

end % end of function LPV2LTV
