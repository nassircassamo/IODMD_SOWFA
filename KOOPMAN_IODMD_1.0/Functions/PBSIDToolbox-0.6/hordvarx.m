function [S,X,DU] = hordvarx(u,y,f,p,reg,opt,noD)
%HORDVARX  Closed-loop Hammerstein system identification using the PBSIDopt method.
%  [S,X,fu]=hordvarx(u,y,f,p) delivers information about the order of the
%  Hammerstein state-space model and acts as a pre-processor for hmodx. The
%  latter is used to identify an open-loop or closed-loop system for the
%  N-by-r and N-by-l and data vectors u, y. The input matrix u and output
%  matrix y must have the same number of observations but can have
%  different numbers of variables. The past and future window size p and f
%  must be higher then the expected order n. The outputs are the singular
%  values S, which can be used to determine the order of the identifiable
%  system. Further is returned the state matrix X and signal vector fu,
%  wich have to be forwarded to hmodx, hx2abcdk, and hx2abck.
%
%  [S,X,fu]=hordvarx(u,y,f,p,reg,opt) adds a regularization to the
%  identification problem. The additional inputs are the regularization
%  method and selection parameters: reg = {'none', 'tikh', 'tsvd'} and opt
%  = {'gcv', 'lcurve', or any regularisation value as scalar}. With
%  regularisation, the solver can better deal with singular covariance
%  matrices. (default reg='none' and opt='gcv')
%
%  [S,X,fu]=hordvarx(u,y,f,p,reg,opt,noD) if noD=1, then the direct
%  feedtrough term D is not considered during estimatoion. Note the direct
%  feedtrough term D can improve the results when estimating low-order
%  models from high-order models. (default noD=0)
%
%  See also: hmodx.m, hx2abcdk.m, and hx2abck.m.
%
%  Note: This function uses the RBF kernel as default.
%
%  References:
%    [1] J.W. van Wingerden, and M. Verhaegen, ``Closed loop
%    identification of MIMO Hammerstein models using LS-SVM'', 15th IFAC
%    Symposium on System Identification, Saint-Malo, France, July 6-8, 2009

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check number if input arguments
if nargin < 4
    error('HORDVARX requires four or five input arguments.')
end

% assign default values to unspecified parameters
if (nargin < 7) || isempty(noD)
    noD = 0;
end
if (nargin < 6) || isempty(opt)
    opt = 'gcv';
end
if (nargin < 5) || isempty(reg)
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
N = size(y,2);
l = size(y,1);
if size(u,2) < size(u,1)
    u = u';
end
r = size(u,1);
if ~isequal(N,length(u))
    error('The number of rows of vectors/matrices u and y must be the same.')
end
if r == 0
    error('HORDVARX requires an output vector u.')
end
if l == 0
    error('HORDVARX requires an output vector y.')
end

% store the past and future vectors
Up = zeros(p*r,N-p);
Yp = zeros(p*l,N-p);
for i = 1:p
    Up((i-1)*r+1:i*r,:) = u(:,i:N+i-p-1);
    Yp((i-1)*l+1:i*l,:) = y(:,i:N+i-p-1);
end

% solve VARX/KERNEL problem
Y = y(:,p+1:N);
Z = zeros(N-p,N-p);
for i=1:p
    for h = 1:r
        Z = Z + kernmatrix(Up((i-1)*r+h,:)','RBF',1,Up((i-1)*r+h,:)') - kernmatrix(Up((i-1)*r+h,:)','RBF',1,zeros(N-p,1));
    end
    Z = Z + Yp((i-1)*l+1:i*l,:)'*Yp((i-1)*l+1:i*l,:);
end
if ~noD
    for h = 1:r
        Z = Z + kernmatrix(u(h,p+1:N)','RBF',1) - kernmatrix(u(h,p+1:N)','RBF',1,zeros(N-p,1));
    end
end
A = regress(Y,Z,reg,opt);

% construct input sequence
LKU = zeros(f*l,N-p);
for i = 0:f-1
    U = zeros(N-p,N-p);
    for j = 1:p-i
        for h = 1:r
            U = U + kernmatrix(Up((i+j-1)*r+h,:)','RBF',1,u(h,p+1:N)') - kernmatrix(Up((i+j-1)*r+h,:)','RBF',1,zeros(N-p,1));
        end
    end
    LKU(i*l+1:(i+1)*l,:) = A*U;
end
[~,S,V] = svd(LKU,'econ');
DU = diag(sqrt(diag(S(1:r,1:r))))*V(:,1:r)';

% construct LambdaKappa
LKZ = zeros(f*l,N-p);
for i = 0:f-1
    Z = zeros(N-p,N-p);
    for j = 1:p-i
        for h = 1:r
            Z = Z + kernmatrix(Up((j-1)*r+h,:)','RBF',1,Up((i+j-1)*r+h,:)') - kernmatrix(Up((j-1)*r+h,:)','RBF',1,zeros(N-p,1));
        end
        Z = Z + Yp((j-1)*l+1:j*l,:)'*Yp((i+j-1)*l+1:(i+j)*l,:);
    end
    LKZ(i*l+1:(i+1)*l,:) = A*Z;
end

% singular value decomposition
[~,S,V] = svd(LKZ,'econ');
X = diag(sqrt(diag(S)))*V';
S = diag(S)';
end




