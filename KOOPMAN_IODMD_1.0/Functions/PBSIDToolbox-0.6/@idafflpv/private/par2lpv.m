function [A,B,C,D]=par2lpv(th,d)
%PAR2LPV    Convert parameter vector to LPV system matrices.
%           par2lpv converts a parameter vector to the
%           system matrices A, B, C, D of the following LPV system:  
%
%           x(k+1) = A(:,1:n)*x(k) + A(:,n+1:n*(q+1))*kron(p(k),x(k))
%                  + B(:,1:m)*u(k) + B(:,m+1:m*(q+1))*kron(p(k),u(k))
%           y(k)   = C(:,1:n)*x(k) + C(:,n+1:n*(q+1))*kron(p(k),x(k))
%                  + D(:,1:m)*u(k) + D(:,m+1:m*(q+1))*kron(p(k),u(k))
%
% Syntax: 
%           [A,B,C,D] = par2lpv(th,dim)
% 
% Input:
%           th       Parameter vector created by lpv2par.
%           dim      Dimension vector created by lpv2par. 
%
% Output:
%           A,B,C,D  System matrices of the LPV state equation
%                    with A: n x n(q+1), B: n x m(q+1),
%                    C: l x n(q+1), and D: l x m(q+1).
% 
% See also: lpv2par, lpvopt.
%
% Written by Vincent Verdult, February 2001.

if size(th,2)>size(th,1)
  th=th';
end
if size(d,2)>size(d,1)
  d=d';
end
if size(th,2)~=1
  error('th must be a vector.')
end
if size(d,2)~=1
  error('dim must be a vector.')
end

n=d(1);
m=d(2);
l=d(3);
s=d(4);

if size(th,1)~=(n+l)*(n+m)*(s+1)
  error('th and dim are not compatible.')
end

temp=reshape(th,n+l,(n+m)*(s+1));
A=temp(1:n,1:n*(s+1));
B=temp(1:n,n*(s+1)+1:(n+m)*(s+1));
C=temp(n+1:n+l,1:n*(s+1));
D=temp(n+1:n+l,n*(s+1)+1:(n+m)*(s+1));

