function X = nuclear(Y,Z,reg_min,tol,maxit,X0)
%NUCLEAR  Matrix regression with nuclear norm optimisation
%  X=NUCLEAR(Y,Z,reg_min) solves the least squares regression problem with
%  nuclear norm regularisation. The inputs are the matrices Y and Z, which
%  have the same number of columns. Output is the least squares estimate X.
%
%  X=NUCLEAR(Y,Z,reg_min,tol) specifies the tolerance of the nonsmooth
%  gradient method. If tol is [] then NUCLEAR uses the default, 1e-4.
%
%  X=NUCLEAR(Y,Z,reg_min,tol,maxit) specifies the maximum number of 
%  iterations. If MAXIT is [] then NUCLEAR uses the default, 200. Warning
%  is given when the maximum iterations is reached.
%
%  X=NUCLEAR(Y,Z,reg_min,tol,maxit,X0) specifies the initial least squares
%  estimate X.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% assign default values to unspecified parameters
if (nargin < 6) || isempty(X0)
    X0 = Y*pinv(Z);
end
if (nargin < 5) || isempty(maxit)
    maxit = 200;
end
if (nargin < 4) || isempty(tol)
    tol = 1e-4;
end

% do the nonsmooth gradient iterations
lambda = 1;
N = length(Y);
X = X0;
ok = true;
k = 1;
cost1 = 1e10;
while ok == true && k <= maxit
    % evaluate function
    E = Y - X*Z;
    cost = norm(E','fro')^2 + reg_min^2*sum(svd(X));
    
    % check residue
    if abs(cost1 - cost) <= tol^2*N
        ok = false;
    end
    
    % recalculate gradient step
    k = k + 1;
    lambda = fminbnd(@(x) (norm((Y-(X + x.*(E*Z'))*Z)','fro')^2 + reg_min^2*sum(svd(X + x.*(E*Z'))))/N,0,lambda);
    
    % step update with tresholding
    X = X + lambda.*(E*Z');
    [U,S,V] = svd(X);
    if size(Y,1) == 1
        X = U*diag(max(S(1,1)-reg_min^2*lambda,0))*V(:,1)';
    else
        X = U*diag(max(diag(S)-reg_min^2*lambda,0))*V(:,1:length(diag(S)))';
    end
    
    % swap
    cost1 = cost;
end