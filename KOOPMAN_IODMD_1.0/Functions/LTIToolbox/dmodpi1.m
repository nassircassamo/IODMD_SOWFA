function [A,B,C,D] = dmodpi4(R,n,stable)
%DMODPI4    Estimates the system matrices of a LTI state space model in
%           innovation form using the output of the dord4 routine. 
%           A and C are found via MOESP, B and D via the impulse matrix.
%           
%
% Model structure:
% 	 	x(k+1) = Ax(k) + Bu(k)
% 		y(k)   = Cx(k) + Du(k) + v(k)
%           where v(k) is zero-mean noise of arbitrary color,
%           independent of the noise-free input u(k).
% 					
% Syntax:
% 	    [A,B,C,D] = dmodpi4(R,n)
% 	    [A,B,C,D] = dmodpi4(R,n,'stable')
% 					
% Input:
%   R       Triangular factor from order detection step.
%   n       Order of system to be estimated. 
%   stable  Estimates a stable A matrix.
% 			
% Output:
%   A,B,C,D Estimated system matrices.
% 					
% See also: dmodpo, dordpo, dmodpi, dordpi

% The modpi routine corresponds to the PI scheme
% derived and analyzed in VERHAEGEN: "Subspace Model
% Identification. Part 3" Int. J. Control, Vol. 57.
%
% Michel Verhaegen December 1994
% Revised by Ivo Houtzager, 2010
% Copyright (c) 1994-2010, Delft Center of Systems and Control 

% Check number of arguments
if nargin < 2
    error('DMODPI1 requires at least two input arguments.');
end
if (nargin < 3) || (isempty(stable))
    stable = 0;
elseif strcmpi(stable,'stable')
    stable = 1;
else
    stable = 0;
end

if (size(R,1)<2) || (size(R,2)<3)
    error('Matrix R has wrong size')
end
m = R(1,2);
l = R(1,3);
s = R(2,3);

if m < 1
    error('Illegal value for number of inputs')
end
if l < 1
    error('Illegal value for number of outputs')
end
if s < 2
    error('Illegal value  for ''s'' in R matrix')
end
if ~((size(R,1)==s*(2*m+l)) && (size(R,2)==2*s*(m+l)))
    error('R-matrix has unexpected size.')
end

if n < 1
    error('System order of zero or lower does not make sense!')
end
if  n >= s
    error('n chosen too large, it should be smaller  than ''s''.')
end

Un = R(2*m*s+1:(2*m+l)*s,(2*m+l)*s+1:2*(m+l)*s);
if stable
    un1 = Un(1:s*l,1:n);
    un2 = [Un(l+1:s*l,1:n); zeros(l,n)];
else
    un1 = Un(1:(s-1)*l,1:n);
    un2 = Un(l+1:s*l,1:n);
end
A = un1\un2;
C = un1(1:l,:);

% least square problems to obtain the B and D matrices
R11 = tril(R(1:s*m,1:s*m));
R31 = R(2*m*s+1:(2*m+l)*s,1:s*m);
Xn = Un(:,n+1:end)'*R31*pinv(R11);
LL = zeros(s*(l*s-n),n+l);
XI = zeros(s*(l*s-n),m);
for i = 1:s
    LL((i-1)*(l*s-n)+1:i*(l*s-n),:) = Un(:,n+1:end)'*[zeros(l*(i-1),n+l); zeros(l,n) eye(l); Un(1:l*(s-i),1:n) zeros(l*(s-i),l)];
    XI((i-1)*(l*s-n)+1:i*(l*s-n),:) = Xn(:,(i-1)*m+1:i*m);
end
S = LL\XI;
B = S(1:n,1:m);
D = S(n+1:n+l,1:m);




