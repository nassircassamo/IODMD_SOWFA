function [Altv,Bltv,Cltv,Dltv] = plpv2ltv(varargin)
%PLPV2LTV  Transforms a PLPV system matrices to a LTV system matrices
%  [Altv,Bltv,Cltv]=plpv2ltv(Aplpv,Bplpv,Cplpv,mu,j) estimates the linear
%  time-varying matrices A, B, and C from the periodic LPV state space
%  model: The inputs are the linear parameter-varying matrices A, B, and C,
%  where matrices have the form of A=[A(1) A(2) ... A(m)]. The outputs are
%  the linear parameter-varying matrices A, B, and C have the form of A{1},
%  A{2}, ..., A{j}. The scheduling sequence mu is periodic with period j.
%
%  [Altv,Bltv,Cltv,Dltv]=plpv2ltv(Aplpv,Bplpv,Cplpv,Dplpv,mu,j) estimates
%  the linear time-varying matrices A, B, C, and D from the periodic LPV
%  state space model: The inputs are the linear parameter-varying matrices
%  A, B, C and D, where matrices have the form of A=[A(1) A(2) ... A(m)].
%  The outputs are the linear parameter-varying matrices A, B, C, and D
%  have the form of A{1}, A{2}, ..., A{j}.
%
%  [Altv,Bltv,Cltv,Dltv,Kltv]=plpv2ltv(Aplpv,Bplpv,Cplpv,Dplpv,Klpv,mu,j)
%  returns also the estimated Kalman matrix K. Note that this Kalman gain
%  does not give a guaranteed stable A-KC (predictor from).
%
%  [Altv,Bltv,Cltv,Dltv,Kltv]=plpv2ltv(Aplpv,Bplpv,Cplpv,Dplpv,x,u,y,mu,j)
%  returns a Kalman matrix K calculated by the discrete-time periodic
%  Riccati equation, which gives a guaranteed stable A-KC (predictor from).

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check input arguments
if nargin < 5
    error('PLPV2LTV requires at least five input arguments!');
end

% assign values to unspecified parameters from command
Aplpv = varargin{1};
Bplpv = varargin{2};
Cplpv = varargin{3};
mu = varargin{nargin-1};
p = varargin{nargin};
Altv = cell(p,1);
Bltv = cell(p,1);
Cltv = cell(p,1);
switch nargin
    case 6
        Dplpv = varargin{4}; 
        Dltv  = cell(p,1);
    case 7
        Dplpv = varargin{4};
        Dltv  = cell(p,1);
        Kplpv = varargin{5};
        Kltv  = cell(p,1);
    case 9
        Dplpv = varargin{4};
        Dltv  = cell(p,1);
        Kltv  = cell(p,1);
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
n = size(Aplpv,1);          % The order of the system
m = size(Aplpv,2)/n - 1;    % The number of scheduling parameters
r = size(Bplpv,2)/(m + 1);  % The number of inputs
l = size(Cplpv,1);          % The number of outputs

% build ltv cells
for i = 1:p
    Altv{i} = Aplpv(1:n,1:n);
    Bltv{i} = Bplpv(1:n,1:r);
    Cltv{i} = Cplpv(1:l,1:n);
    if exist('Dlpv','var')
        Dltv{i} = Dplpv(1:l,1:r);
    end
    if exist('Klpv','var')
        Kltv{i} = Kplpv(1:l,1:r);
    end
    for j = 1:m
        Altv{i} = Altv{i} + Aplpv(1:n,(j-1)*n+1+n:j*n+n).*mu(i,j);
        Bltv{i} = Bltv{i} + Bplpv(1:n,(j-1)*r+1+r:j*r+r).*mu(i,j);
        Cltv{i} = Cltv{i} + Cplpv(1:l,(j-1)*n+1+n:j*n+n).*mu(i,j);
        if exist('Dlpv','var')
            Dltv{i} = Dltv{i} + Dplpv(1:l,(j-1)*r+1+r:j*r+r).*mu(i,j);
        end
        if exist('Klpv','var')
            Kltv{i} = Kltv{i} + Kplpv(1:l,(j-1)*r+1+r:j*r+r).*mu(i,j);
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
    for i = 1:p
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
    [Pp,Gp] = dpre(Ahat,Chat,Q,R,S);
    for i = 1:p
        Kltv{i} = Gp(:,:,i)';
    end
end

end % end of function PLPV2LTV