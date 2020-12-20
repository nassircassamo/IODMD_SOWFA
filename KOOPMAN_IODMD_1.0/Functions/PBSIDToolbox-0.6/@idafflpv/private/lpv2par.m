function [th,d]=lpv2par(A,B,C,D)
%LPV2PAR    Convert  LPV system matrices to parameter vector.
%           lpv2par converts the system matrices A, B, C, D of 
%           the following LPV system to a parameter vector to
%           be used with lpvopt:
%
%           x(k+1) = A(:,1:n)*x(k) + A(:,n+1:n*(q+1))*kron(p(k),x(k))
%                  + B(:,1:m)*u(k) + B(:,m+1:m*(q+1))*kron(p(k),u(k))
%           y(k)   = C(:,1:n)*x(k) + C(:,n+1:n*(q+1))*kron(p(k),x(k))
%                  + D(:,1:m)*u(k) + D(:,m+1:m*(q+1))*kron(p(k),u(k))
%
% Syntax: 
%           [th,dim] = lpv2par(A,B,C,D)
% 
% Input:
%           A,B,C,D  System matrices of the LPV state equation
%                    with A: n x n(q+1), B: n x m(q+1),
%                    C: l x n(q+1), and D: l x m(q+1).
%
% Output:
%           th       Parameter vector.
%           dim      Dimension vector containing the dimensions of
%                    x(k), u(k), y(k) and p(k).
% 
% See also: lpv2par, lpvopt.

% Written by Vincent Verdult, February 2001.

n=size(A,1); 
sn=size(A,2);
sm=size(B,2);
l=size(C,1);

if size(B,1)~=n
  error('A and B must have the same number of rows.')
end
if size(C,2)~=sn
  error('A and C must have the same number of columns.')
end
if size(D,2)~=sm
  error('B and D must have the same number of columns.')
end
if size(D,1)~=l
  error('C and D must have the same number of rows.')
end
s=sn/n-1;
if s~=fix(s)
  error('A matrix is of wrong size.')
end
m=sm/(s+1);
if m~=fix(m)
  error('B matrix is of wrong size.')
end

d=[n;m;l;s];
th=reshape([A B; C D],(n+l)*(sm+sn),1);



