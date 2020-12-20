function [y,x]=lpvsim(A,B,C,D,p,u,w,x0)
%LPVSIM     Simulate LPV system.
%           lpvsim computes the output y(k) for k=1,2,...,N
%           of the linear parameter varying (LPV) system:
%
%           x(k+1) = A(:,1:n)*x(k) + A(:,n+1:n*(q+1))*kron(p(k),x(k))
%                  + B(:,1:m)*u(k) + B(:,m+1:m*(q+1))*kron(p(k),u(k)) + w(k)
%           y(k)   = C(:,1:n)*x(k) + C(:,n+1:n*(q+1))*kron(p(k),x(k))
%                  + D(:,1:m)*u(k) + D(:,m+1:m*(q+1))*kron(p(k),u(k))
%
% Syntax: 
%           x = lpvsim(A,B,C,D,p,u,w,x0)
%
% Input:
%           A,B,C,D  System matrices of the LPV state equation
%                    with A: n x n(q+1), B: n x m(q+1),
%                    C: l x n(q+1), and D: l x m(q+1).
%           p        N x q matrix containing N data points of the q
%                    time varying parameters,
%           u        N x m matrix containing N data points of the m
%                    inputs.
%           w        N x n matrix containing N data points of the 
%                    process noise that has dimension n.
%           x0       column vector containing the initial state. 
%
% Output:
%           y        N x l matrix containing the computed outputs.
%           x        N x n matrix containing the computed state.
%
% See also: lpvitr, bilitr.

% Written by Vincent Verdult, February 2001.

% check input arguments
if nargin<7
  error('Not enough input arguments.')
end

n=size(A,1);
l=size(C,1);
ns=size(A,2);
ms=size(B,2);
N=size(p,1);

if nargin<8
  x0=zeros(n,1);
end
if size(D,1)~=l
  error('C and D must have the same number of rows.')
end
if size(C,2)~=ns
  error('A and C must have the same number of columns.')
end
if size(D,2)~=ms
  error('B and D must have the same number of columns.')
end
% remaining checks will be dealt with in lpvitr

x=lpvitr(A,B,p,u,w,x0);

s=ns/n-1;
m=ms/(s+1);
pu=zeros(N,s*m);
px=zeros(N,s*n);
for i = 1:s
  pu(:,(i-1)*m+1:i*m) = u.*(p(:,i)*ones(1,m));
end 
for i = 1:s;
  px(:,(i-1)*n+1:i*n) = x.*(p(:,i)*ones(1,n));
end

y=[x, px]*C'+[u, pu]*D';


