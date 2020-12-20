function [P,sigma] = dvar4varx(u,y,p,VARX,Zps)
%DVAR4VARX Asymptotic variance of the VARX estimation
%  P=dvar4varx(u,y,p,VARX,Zps) returns the covariance of the VARX
%  estimation and acts as a pre-processor for dvar2frd. The latter is used
%  to calculate the probalistic error bounds around the identified bode
%  diagrams. The input matrix u, output matrix y, past window p, regression
%  model VARX, and data matrix Zps (obtained from dordvarx) must have the
%  same number of observations.
%
%  [P,sigma]=dvar4varx(u,y,p,VARX,Zps) also returns the covariance matrix
%  of the innovation noise. 

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check number if input arguments
if nargin < 5
    error('DVAR4VARX requires five input arguments.')
end

% check dimensions of inputs
if size(y,2) < size(y,1)
    y = y';
end
N = size(y,2);
l = size(y,1);
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
if l == 0
    error('DVAR4VARX requires an output vector y.')
end

% store the past and future vectors
m = r+l;
z = [u; y];
Z = zeros(p*m,N-p);
for i = 1:p
    Z((i-1)*m+1:i*m,:) = z(:,i:N+i-p-1);
end

% solve VARX problem
Y = y(:,p+1:N);
U = u(:,p+1:N);
if size(VARX,2)/m > p
    Z = [Z; U];
end

% calculate the innovation sequence
E = Y - VARX*Z;
sigma = (E*E')/length(E);

P = zeros(l*size(Zps,2),l*size(Zps,2));
% asymptotic variance
% P = kron(Zps',eye(l))*sparse(kron(speye(N-p),sigma))*kron(Zps',eye(l))';
% computed sample-by-sample for memory efficiency 
P = zeros(l*size(Zps,2),l*size(Zps,2));
for j = 1:size(Zps,1)
    P = P + kron(Zps(j,:)',eye(l))*sigma*kron(Zps(j,:),eye(l));        
end



