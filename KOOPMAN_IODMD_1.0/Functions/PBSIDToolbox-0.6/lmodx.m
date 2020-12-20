function x = lmodx(X,n)
%LMODX  Closed-loop LPV system identification using the PBSIDopt method.
%  x=lmodx(X,n) estimates the state sequence x of the identifiable system
%  with order n. The order n can be determined from the singular values
%  given by lordvarx. The matrix X is also calculated by hordvarx.
%
%  See also: lordvarx, lx2abcdk.m, and lx2abck.m.
%
%  References:
%    [1] J.W. van Wingerden, and M. Verhaegen, ``Subspace identification
%    of Bilinear and LPV systems for open- and closed-loop data'',
%    Automatica 45, pp 372--381, 2009.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check number of arguments
if nargin < 2
    error('LMODX requires at least two input arguments.');
end

% cheack the dimensions of the inputs
mx = size(X,1);
if (n < 1) || isempty(n)
    error('System order of zero or lower does not make sense!')
end
if mx < n
    error('The number of rows of matrix X must be equal or higher then the order n.')
end

x = X(1:n,:);

% scale the state sequence
%x = diag(1./sqrt(var(x,0,2)))*x;





