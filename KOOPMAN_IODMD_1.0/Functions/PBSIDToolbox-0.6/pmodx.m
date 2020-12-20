function [X,CC,vx,vy] = pmodx(X,TT,K,n,vx,vy)
%PMODX  Periodic LPV system identification using the PBSIDopt method.
%  [X,CC] = pmodx(X,UM,K,n) estimates the state sequence X of the
%  identifiable system with order n. The order n can be determined from the
%  correlation coefficients C given by dordopt. The matrices UM and K are
%  also calculated by pordvarx.
%
%  [X,CC] = pmodx(X,UM,K,n,vx,vy) specifies the regularized parameters vx
%  and vy for the regularized canonical correlation analysis solver. For vx
%  > 0 and vy > 0, the solver can better deal with ill-conditioned
%  covariance matrices. To use the k-fold cross validation to calculate the
%  regularized parameter with best score, choose the vector vx as for
%  example vx = logspace(-8,-2). For vectors with lengths larger the two,
%  the code automatically uses cross validation. (default vx=0 and vy=0)
%
%  [X,CC,vxm,vym] = pmodx(X,UM,K,n,vx,vy) also returns the best regularized
%  parameter found by the k-fold cross validation.
%
%  References:
%    [1] van Wingerden, J.W., Houtzager, I., Verhaegen, M., Closed-loop
%    identification of the time-varying dynamics of variable-speed wind
%    turbines, Int. J. Robust Nonlinear Control 2008
%
%  See also: pordxarx, px2abcdk, and px2abck.m.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check number of arguments
if nargin < 3
    error('PMODX requires at least three input arguments.');
end

% cheack the dimensions of the inputs
if (n < 1) || isempty(n)
    error('System order of zero or lower does not make sense!')
end

% allocate matrices
j = size(X,1);
fl = size(TT{1},1);
UM = zeros(j*fl,j*n);
for q = 1:j
    UM((q-1)*fl+1:q*fl,(q-1)*n+1:q*n) = TT{q,1}(:,1:n);
end
clear TT;

% solve intersection problem
if (length(vx) > 1) || (length(vy) > 1)
    [vx,vy] = kfcv(UM,K,q,vx,vy);
end
[CC,TU] = rcca(UM,K,vx,vy);

% obtain the state sequence using the transformation matrix
TU = TU(1:j*n,1:n);
for q = 1:j
    X{q,1} = TU((q-1)*n+1:n*q,:)\X{q,1}(1:n,:);
end





