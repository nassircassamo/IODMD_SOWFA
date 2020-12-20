function [r,V,W] = rcca(X,Y,vx,vy)
%RCCA Regularized Canonical Correlation Analysis
%  [R,V,W]=RCCA(X,Y) solves the canonical correlation analysis problem. The
%  inputs are the matrices X and Y, which have the same number of rows.
%  Outputs are the canonical factor matrices V and W and the vector R with
%  the canonical correlation coefficents.
%
%  [R,V,W]=RCCA(X,Y,V1,V2) solves the regularized canonical correlation
%  analysis problem. The additional inputs are the regularized parameters
%  V1 and V2. If V > 0, the solver can better deal with singular covariance
%  matrices. (default V1=0 and V2=0)
%
%  See also CANONCORR.

%  References:
%    [1] Krzanowski, W.J., Principles of Multivariate Analysis,
%        Oxford University Press, Oxford, 1988.
%    [2] Seber, G.A.F., Multivariate Observations, Wiley, New York, 1984.

%  revised version of CANONCORR with the implementation of regularized 
%  parameters (original CANONCORR is in the statistical toolbox)

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% assign default values to unspecified parameters
if (nargin < 4) || isempty(vx)
    vx = 0;
end
if (nargin < 3) || isempty(vy)
    vy = 0;
end

% check input arguments
if nargin < 2
    error('RCCA requires at least two input arguments')
end

% sizes
[rx,cx] = size(X);
[ry,cy] = size(Y);

% check size matrices
if ~isequal(rx,ry);
    error('RCCA requires that the rows of X are equal to the rows of Y');
end

% center the row variables
X = X - repmat(mean(X,1), rx, 1);
Y = Y - repmat(mean(Y,1), ry, 1);

% decompose and find an orthonormal basis
[Qx,Rx,permx] = qr(X,0);
rankX = sum(abs(diag(Rx)) > eps(abs(Rx(1)))*max(rx,cx));
if rankX < cx
    %warning('RCCA:notFullRank','X is not full rank.');
    Qx = Qx(:,1:rankX);
    Rx = Rx(1:rankX,1:rankX);
end
[Qy,Ry,permy] = qr(Y,0);
rankY = sum(abs(diag(Ry)) > eps(abs(Ry(1)))*max(ry,cy));
if rankY < cy
    %warning('RCCA:notFullRank','Y is not full rank.');
    Qy = Qy(:,1:rankY);
    Ry = Ry(1:rankY,1:rankY);
end
d = min(rankX,rankY);
    
if vx == 0 && vy == 0
    % compute canonical coefficients and canonical correlations
    [L,D,M] = svd(Qx'*Qy,0);
    A = Rx\L(:,1:d)*sqrt(rx-1);
    B = Ry\M(:,1:d)*sqrt(ry-1);
else
    % add the regularized parameters
    Ax = real(sqrtm(Rx'*Rx + vx*eye(rankX)));
    Ay = real(sqrtm(Ry'*Ry + vy*eye(rankY)));

    % compute canonical coefficients and canonical correlations
    [L,D,M] = svd(Ax\(Qx*Rx)'*(Qy*Ry)/Ay,0);
    A = Ax\L(:,1:d)*sqrt(rx-1);
    B = Ay\M(:,1:d)*sqrt(ry-1);
end

% remove roundoff errors
r = min(max(diag(D(:,1:d))',0),1); 

% put coefficients back to their full size and their correct order
V(permx,:) = [A; zeros(cx-rankX,d)];
W(permy,:) = [B; zeros(cy-rankY,d)];


