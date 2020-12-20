function [A,B,C,D,K] = dmodpo1(R,n,stable)
%DMODPO4    Estimates the system matrices of a LTI state space model in
%           innovation form using the output of the dordpo routine. 
%           A and C are found via MOESP, B and D via the impulse matrix.
%           
%
% Model structure:
% 	 	x(k+1) = Ax(k) + Bu(k) + w(k)
% 		y(k)   = Cx(k) + Du(k) + v(k)
%           where w(k) and v(k) are zero-mean white noise sequences,
%           independent of the noise-free input u(k). When a fifth 
%           output variable is specified, the Kalman gain K is
%           calculated. The Kalman gain can be used to construct the
%           one-step ahead predictor:
%             x'(k+1)  =  Ax'(k) + (Bu(k) + K(y(k)-yest(k))
%             yest(k)  =  Cx'(k) + Du(k)
% 					
% Syntax:
% 	    [A,B,C,D] = dmodpo1(R,n)
% 	    [A,B,C,D] = dmodpo1(R,n,'stable')
% 	    [A,B,C,D,K] = dmodpo1(R,n)
% 					
% Input:
%   R       Triangular factor from order detection step.
%   n       Order of system to be estimated. 
%   stable  Estimates a stable A matrix.
% 			
% Output:
%   A,B,C,D Estimated system matrices.
%   K       Estimated Kalman gain.
% 					
% See also: dmodpo, dordpo, dmodpi, dordpi

% Michel Verhaegen December 1994
% Revised by Ivo Houtzager, 2010
% Copyright (c) 1994-2010, Delft Center of Systems and Control 

% Check number of arguments
if nargin < 2
    error('DMODPO1 requires at least two input arguments.');
end
if (nargin < 3) || (isempty(stable))
    stable = 0;
elseif strcmpi(stable,'stable')
    stable = 1;
else
    stable = 0;
end
m = R(1,2);
l = R(1,3);
s = R(1,4);

if (size(R,1)<2*m || size(R,2)<3) && (m ~= 0) 
    error('Matrix R has wrong size')
end
if m <= 0
    error('Illegal value for number of inputs in R matrix')
end
if l < 1
    error('Illegal value for number of outputs in R matrix ')
end
if s < 0
    error('Illegal value  for ''s'' in R matrix')
end
if ~((size(R,1)==(2*m+2*l)*s) && (size(R,2)==max(4,(2*m+3*l)*s)))
    error('R-matrix has unexpected size.')
end
if n < 1
    error('System order of zero or lower does not make sense!')
end
if n >= s
    error('n chosen too large, it should be smaller than s.')
end

Un = R((2*m+l)*s+1:(2*m+2*l)*s,(2*m+2*l)*s+1:(2*m+3*l)*s);
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
R31 = R((2*m+l)*s+1:2*(m+l)*s,1:s*m);
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

% calculate Kalman gain
if nargout== 3
    R1 = R(1:m*s,1:m*s);
    R2 = R((2*m+l)*s+1:2*(m+l)*s,1:(2*m+l)*s+l);
    Zip1 = R2(l+1:s*l,1:(2*m+l)*s+l);
    Zi   = [R2(1:s*l,1:(2*m+l)*s) zeros(s*l,l)];
    Gam  = Un(:,1:n);

    LHS    = [Gam(1:(s-1)*l,1:n)\Zip1;R2(1:l,1:(2*m+l)*s+l)];
    RHS    = [Gam\Zi  ;[R1 zeros(m*s,(m+l)*s+l)]];
    K      = LHS/RHS;
    Resid  = LHS - K*RHS;
    
    % Q,R,S matrices
    SOL = Resid*Resid';

    Q = SOL(1:n,1:n);
    S = SOL(1:n,n+1:n+l);
    R = SOL(n+1:n+l,n+1:n+l);
    lambda = min(eig([Q S; S' R]));
    if lambda < 0
        if lambda < -1e-6
            warning('LTI:zeroK','Outputing zero K')
            K = zeros(n,l);
        else
            warning('LTI:adjustR','Adjusting R, to make it positive definite')
            R = R-eye(length(R))*lambda*2;
            Q = Q-eye(length(Q))*lambda*2;
            [X,L,K] = dare(A',C',Q,R,S);
            K = K';
        end
    else
        [X,L,K] = dare(A',C',Q,R,S);
        K = K';
    end
end


