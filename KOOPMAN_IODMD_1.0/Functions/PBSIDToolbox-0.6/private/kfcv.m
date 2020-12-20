function [vmx,vmy,CVm,CV] = kfcv(X,Y,k,vx,vy)
%KFCV K-Fold Cross Validation
%  [VX,VY]=KFCV(X,Y,K,VX,VY) determines the best regularized parameters VX 
%  and VY for the matrices X and Y. The scalar K is the divider which 
%  divides the X and Y in X1,X2...,Xk and Y1,Y2,...,Yk blocks. The best 
%  parameters are determined by cross validation of the X and Y matrices 
%  minus the blocks Xk and Yk. The input vector VX and VY contains all the 
%  evaluated values for these parameters. The best parameters with highest 
%  score are returned to the output VX and VY.   
%
%  [VX,VY,CVM,CV]=KFCV(X,Y,K,VX,VY) also returns the best score CVM. 
%
%  [VX,VY,CVM,CV]=KFCV(X,Y,K,VX,VY) returns the matrix CV which contains 
%  all the scores for the values in the input vectors VX and VY.
%
%  See also RCCA.

%  References:
%    [1] Krzanowski, W.J., Principles of Multivariate Analysis,
%        Oxford University Press, Oxford, 1988.
%    [2] Seber, G.A.F., Multivariate Observations, Wiley, New York, 1984.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% assign default values to unspecified parameters
if (nargin < 5) || isempty(vx)
    vy = logspace(-8,-2);
end
if (nargin < 4) || isempty(vy)
    vx = logspace(-8,-2);
end

% check input arguments
if nargin < 3
    error('KFCV requires at least three input arguments')
end

% sizes
[rx,cx] = size(X);
[ry,cy] = size(Y);
dl = rx/k; % block size

% check size matrices
if ~isequal(rx,ry);
    error('KFCV requires that the rows of X are equal to the rows of Y');
end

% center the variables
X = X - repmat(mean(X,1), rx, 1);
Y = Y - repmat(mean(Y,1), ry, 1);

% decompose and find an orthonormal basis
[Qx,Rx,permx] = qr(X,0);
rankX = sum(abs(diag(Rx)) > eps(abs(Rx(1)))*max(rx,cx));
if rankX < cx
    %warning('KFCV:notFullRank','X is not full rank.');
    Qx = Qx(:,1:rankX);
    Rx = Rx(1:rankX,1:rankX);
end
[Qy,Ry,permy] = qr(Y,0);
rankY = sum(abs(diag(Ry)) > eps(abs(Ry(1)))*max(ry,cy));
if rankY < cy
    %warning('KFCV:notFullRank','Y is not full rank.');
    Qy = Qy(:,1:rankY);
    Ry = Ry(1:rankY,1:rankY);
end

% full rank X and Y
X = Qx*Rx;
Y = Qy*Ry;

% for all parameters values, calculate the score
CV = zeros(length(vx),length(vy));
Cvx = [];
Cvy = [];
for h = 1:length(vx)
    for i = 1:length(vy)
        for j = 1:k
            % remove part K from X and Y 
            if j == 1
                Xk = X(k+1:dl*k,:);
                Yk = Y(k+1:dl*k,:);
            elseif j == k
                Xk = X(1:(dl-1)*k,:);
                Yk = Y(1:(dl-1)*k,:);
            else
                Xk = X([1:(i-1)*k i*k+1:dl*k],:);
                Yk = Y([1:(i-1)*k i*k+1:dl*k],:);
            end

            % add the regularized parameters    
            Ax = real(sqrtm(Xk'*Xk + vx(h)*eye(rankX)));
            Ay = real(sqrtm(Yk'*Yk + vy(i)*eye(rankY)));

            % compute canonical coefficients and canonical correlations
            [L,D,M] = svd(Ax\Xk'*Yk/Ay,0);
            a = Ax\L(:,1)*sqrt(rx-1);
            b = Ay\M(:,1)*sqrt(ry-1);
            Cvx(j,:) = (Xk*a)';
            Cvy(j,:) = (Yk*b)';
        end

        % calculate cross correlation coeffients
        Cov = corrcoef(Cvx,Cvy);
        CV(h,i) = Cov(2,1);
    end
end

% find maximum score
npoints = 200;
xi = logspace(log10(min(vx)),log10(max(vx)),npoints);
yi = logspace(log10(min(vy)),log10(max(vy)),npoints);
CV = interp2(vy,vx,CV,yi,xi,'spline');
[CVm,vh] = max(CV);
[CVm,vi] = max(CVm);

% return the parameter values for max score
vmx = xi(vh(1));
vmy = yi(vi(1));