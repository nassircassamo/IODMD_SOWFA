function X = kpinv(A,B,varargin)
%KPINV Pseudoinverse of Kronecker product.
%  X=KPINV(A,B) produces a matrix X of the same dimensions as KRON(A',B')
%  so that KRON(A,B)*X*KRON(A,B)=KRON(A,B), X*KRON(A,B)*X=X and KRON(A,B)*X
%  and X*KRON(A,B) are Hermitian. The computation is based on A=U1*S1*V1'
%  and B=U2*S2*V2', because KRON(A,B)=KRON(U2,U2)*KRON(S1,S2)*KRON(V1,V2)'.
%  Any singular values less than a tolerance are treated as zero. The
%  default tolerances are MAX(SIZE(KRON(A,B)))*EPS(class(A)).
%
%  KPINV(A,B,TOL) uses the tolerance TOL instead of the default.
%
%  See also PINV.

%  Revised version of PINV for kronecker products, is faster and needs less
%  memory especially when large matrices are being used.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% sizes
[m1,n1] = size(A);
[m2,n2] = size(B);
m = m1*m2;
n = n1*n2;

% quick return if matrices are empty
if isempty(A)
  X = zeros(n,m,class(A));
  return
end

if n > m
   X = kpinv(A',B',varargin{:})';
else
   [U1,S1,V1] = svd(A,0);
   [U2,S2,V2] = svd(B,0);

   if m > 1,
       s1 = diag(S1);
       s2 = diag(S2);
       s = kron(s1,s2);
   elseif m == 1,
       s1 = S1(1);
       s2 = S2(1);
       s = s1*s2;
   else
       s1 = 0;
       s2 = 0;
       s = 0;
   end

   if nargin == 3
      tol1 = varargin{1};
      tol2 = tol1;
   else
      tol1 = max(m1,n1)*eps(max(s1));
      tol2 = max(m2,n2)*eps(max(s2));
   end

   r1 = sum(s1 > tol1);
   r2 = sum(s2 > tol2);
   r = r1*r2;
   if (r == 0)
      X = zeros(n,m,class(A));
   else
      s = diag(ones(r,1)./s(1:r));
      X = kron(V1(:,1:r1),V2(:,1:r2))*s*kron(U1(:,1:r1)',U2(:,1:r2)');
   end
end
