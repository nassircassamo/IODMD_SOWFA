function [A,B,C,D] = dmodpi4(R,u,y,n,stable)
%DMODPI4    Estimates the system matrices of a LTI state space model in
%           innovation form using the output of the dord4 routine. 
%           A and C are found via MOESP, B and D via N4SID.
%           
%
% Model structure:
% 	 	x(k+1) = Ax(k) + Bu(k)
% 		y(k)   = Cx(k) + Du(k) + v(k)
%           where v(k) is zero-mean noise of arbitrary color,
%           independent of the noise-free input u(k).
% 					
% Syntax:
% 	    [A,B,C,D] = dmodpi4(R,u,y,n)
% 	    [A,B,C,D] = dmodpi4(R,u,y,n,'stable')
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

% Written by Ivo Houtzager, 2010
% Copyright (c) 1994-2010, Delft Center of Systems and Control 

% Check number of arguments
if nargin < 4
    error('DMODPI4 requires at least four input arguments.');
end
if (nargin < 5) || (isempty(stable))
    stable = 0;
elseif strcmpi(stable,'stable')
    stable = 1;
else
    stable = 0;
end
if size(u,2) > size(y,1)
    y = y';
end
if size(u,2) > size(u,1)
    u = u';
end
N = size(u,1);

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
if n >= s
    error('n chosen too large, it should be smaller  than ''s''.')
end

% get A and C
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

% get B and D using the prediction of the state
NN = N-2*s+1;
Up = zeros(NN,m*s);
for i = (1:s),
    Up(:,(i-1)*m+1:i*m) = u(i:NN+i-1,:);
end 
R32 = R(2*m*s+1:(2*m+l)*s,m*s+1:2*m*s);
R22 = tril(R(m*s+1:2*m*s,m*s+1:2*m*s));
Xn = pinv(Un(1:s*l,1:n))*R32*pinv(R22)*Up';
B = (Xn(:,2:end) - A*Xn(:,1:end-1))*pinv(u(s+1:NN+s-1,:)');
D = (y(s+1:NN+s,:)' - C*Xn)*pinv(u(s+1:NN+s,:)');





