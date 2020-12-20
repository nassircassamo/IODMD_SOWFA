function [S,X,DY,cc] = wordvarx(u,y,f,p,vx,vy,noD)
%HORDVARX  Closed-loop Wiener system identification using the PBSIDopt method.
%  [S,X,g1y]=wordvarx(u,y,f,p) delivers information about the order of the
%  Wiener state-space model and acts as a pre-processor for wmodx. The
%  latter is used to identify an open-loop or closed-loop system for the
%  N-by-r and N-by-l and data vectors u, y. The input matrix u and output
%  matrix y must have the same number of observations but can have
%  different numbers of variables. The past and future window size p and f
%  must be higher then the expected order n. The outputs are the singular
%  values S, which can be used to determine the order of the identifiable
%  system. Further is returned the state matrix X and signal vector g1y,
%  wich have to be forwarded to wmodx, wx2abcdk, and wx2abck.
%
%  [S,X,g1y]=wordvarx(u,y,f,p,vx,vy)  specifies the regularized
%  parameters vx and vy for the regularized canonical correlation analysis
%  solver. For vx > 0 and vy > 0, the solver can better deal with
%  ill-conditioned covariance matrices. To use the k-fold cross validation
%  to calculate the regularized parameter with best score, choose the
%  vector vx as for example vx = logspace(-8,-2). For vectors with lengths
%  larger the two, the code automatically uses cross validation. (default
%  vx=0 and vy=0)
%
%  [S,X,g1y]=wordvarx(u,y,f,p,vx,vy,noD) if noD=1, then the direct
%  feedtrough term D is not considered during estimatoion. Note the direct
%  feedtrough term D can improve the results when estimating low-order
%  models from high-order models. (default noD=0)
%
%  See also: wmodx.m, wx2abcdk.m, and wx2abck.m.
%
%  Note: This function uses the RBF kernel as default.
%
%  References:
%    [1] J.W. van Wingerden, and M. Verhaegen, ``Closed-loop subspace
%    identification of Hammerstein-Wiener models'', Joint 48th IEEE
%    Conference on Decision and Control and 28th Chinese Control Conference
%    Shanghai, P.R. China, December 16-18, 2009.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010


% check number if input arguments
if nargin < 4
    error('WORDVARX requires four or five input arguments.')
end

% assign default values to unspecified parameters
if (nargin < 7) || isempty(noD)
    noD = 0;
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
    error('WORDVARX requires an output vector u.')
end
if l == 0
    error('WORDVARX requires an output vector y.')
end

% store the past and future vectors
Up = zeros(p*r,N-p);
Yp = zeros(p*l,N-p);
for i = 1:p
    Up((i-1)*r+1:i*r,:) = u(:,i:N+i-p-1);
    Yp((i-1)*l+1:i*l,:) = y(:,i:N+i-p-1);
end

% solve VARX/KERNEL problem
Y = zeros(N-p,N-p);
U = zeros(N-p,N-p);
for i=1:p
    U = U + Up((i-1)*r+1:i*r,:)'*Up((i-1)*r+1:i*r,:);
    for h = 1:l
        Y = Y + kernmatrix(Yp((i-1)*l+h,:)','RBF',1) - kernmatrix(Yp((i-1)*l+h,:)','RBF',1,zeros(N-p,1));
    end
end
if ~noD
    U = U + u(:,p+1:N)'*u(:,p+1:N);
end
Yf = zeros(N-p,N-p);
for h = 1:r
    Yf = Yf + kernmatrix(y(h,p+1:N)','RBF',1) - kernmatrix(y(h,p+1:N)','RBF',1,zeros(N-p,1));
end
Y = Yf - Y; 

% solve intersection problem
if (length(vx) > 1) || (length(vy) > 1)
    [vx,vy] = kfcv(Y',U',p,vx,vy);
    fprintf('K-fold selects the regularisation values vy = %6.4e and vu = %6.4e \n',vx,vy);
end
[cc,A,B] = rcca(Y',U',vx,vy);
A = A(:,1:l)';
B = B(:,1:l)';
clear Y

% construct output sequence
DY = A*Yf;
clear Yf;

% construct LambdaKappa
LKZ = zeros(f*l,N-p);
for i = 0:f-1
    U = zeros(N-p,N-p);
    Y = zeros(N-p,N-p);
    for j = 1:p-i
        U = U + Up((j-1)*r+1:j*r,:)'*Up((i+j-1)*r+1:(i+j)*r,:);
        for h = 1:r
            Y = Y + kernmatrix(Yp((j-1)*l+h,:)','RBF',1,Yp((i+j-1)*l+h,:)') - kernmatrix(Yp((j-1)*r+h,:)','RBF',1,zeros(N-p,1));
        end
    end
    LKZ(i*l+1:(i+1)*l,:) = B*U + A*Y;
end

% singular value decomposition
[~,S,V] = svd(LKZ,'econ');
X = diag(sqrt(diag(S)))*V';
S = diag(S)';
end




