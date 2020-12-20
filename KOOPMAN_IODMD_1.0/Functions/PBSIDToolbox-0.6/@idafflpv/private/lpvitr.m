function x=lpvitr(A,B,p,u,w,x0)
%LPVITR     Iterate state equation of an LPV system.
%           lpvitr computes the state x(k) for k=1,2,...,N
%           satisfying the linear parameter-varying (LPV) state equation:
%
%           x(k+1) = A(:,1:n)*x(k) + A(:,n+1:n*(q+1))*kron(p(k),x(k))
%                  + B(:,1:m)*u(k) + B(:,m+1:m*(q+1))*kron(p(k),u(k))
%                  + w(k)
%
% Syntax: 
%           x = lpvitr(A,B,p,u,w,x0)
%
% Input:
%           A,B    System matrices of the LPV state equation
%                  with A: n x n(q+1) and B: n x m(q+1).
%           p      N x q matrix containing N data points of the q
%                  time varying parameters,
%           u      N x m matrix containing N data points of the m
%                  inputs.
%           w      N x n matrix containing N data points of the 
%                  process noise that has dimension n.
%           x0     column vector containing the initial state. 
%
% Output:
%           x     N x n matrix containing the computed state.
%
% See also: lpvsim, bilitr.

% Written by Vincent Verdult, February 2001.

% check input arguments
if nargin<6
    error('Not enough input arguments.')
end

[N,m]=size(u);
n=size(A,1);
s=size(p,2);

if size(A,2)~=n*(s+1)
    error('A matrix has wrong number of columns.');
end
if size(B,1)~=n
    error('B matrix has wrong number of rows.');
end
if size(B,2)~=m*(s+1)
    error('B matrix has wrong number of columns.');
end
if size(p,1)~=N
    error('p and u must have the same number of rows.');
end
if (~isempty(w) && size(w,2)~=n)
    error('w has wrong number of columns.');
end
if (~isempty(w) && size(w,1)~=N)
    error('p and w must have the same number of rows.');
end
if size(x0,1)~=n
    error('x0 and A must have the same number of rows.');
end
if size(x0,2)~=1
    error('x0 can only have one column.');
end

% allocate memory
x=zeros(n,N);
px=zeros(n*s,1);
pu=zeros(m*s,1);

if isempty(w)
    % without process noise
    for k=1:N
        x(:,k)=x0;
        for i=1:s
            px((i-1)*n+1:i*n,1)=p(k,i)*x0;
            pu((i-1)*m+1:i*m,1)=p(k,i)*u(k,:)';
        end
        x0=A(:,1:n)*x0+A(:,n+1:n*(s+1))*px+B(:,1:m)*u(k,:)'+B(:,m+1:m*(s+1))*pu;
    end
else
    % with process noise
    for k=1:N
        x(:,k)=x0;
        for i=1:s
            px((i-1)*n+1:i*n,1)=p(k,i)*x0;
            pu((i-1)*m+1:i*m,1)=p(k,i)*u(k,:)';
        end
        x0=A(:,1:n)*x0+A(:,n+1:n*(s+1))*px+B(:,1:m)*u(k,:)'+B(:,m+1:m*(s+1))*pu+w(k,:)';
    end
end
x=x';
