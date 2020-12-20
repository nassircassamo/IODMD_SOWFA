function x = ltiitr(A,B,u,w,x0) 
%LTIITR    Iterates the state equation of an LPV system 
%          Computes the state x(k) for k=1,2,...,N 
%          satisfying the linear time-inavriant (LPV) 
%          state equation: 
% 
%          x(k+1) = A(:,1:n)*x(k) + B(:,1:m)*u(k) + w(k) 
% 
%          This function is used internally by dfunlti and is not 
%          meant for stand-alone use. 
% 
% Syntax: 
%          x=ltiitr(A,B,u,w,x0) 
% 
% Inputs: 
%  A       An LTI state-transition matrix of size n-by-n
%  B       An LTI input matrix of size n-by-m. 
%  u       A N-by-m matrix containing N samples of the m inputs. 
%  w       (optional) A N-by-n matrix containing the process noise. 
%  x0      The initial state, a n-by-1 vector 
% 
% Outputs: 
%  x       The computed state, an N-by-n matrix. 
% 
% See also: dfunlti 

% Based on LPVITR
% Written by Vincent Verdult, February 2001 
% Revised by Niek Bergboer, October 2001 
% Revised by Ivo Houtzager, 2007
% Copyright (c) 2001-2007, Delft Center of Systems and Control 

if nargin < 5
    error('LTIITR requires five input arguments.')
end
[N,m] = size(u);
n = size(A,1);
if size(A,2) ~= n
    error('A matrix has wrong number of columns.');
end
if size(B,1) ~= n
    error('B matrix has wrong number of rows.');
end

if (~isempty(w) && size(w,2)~=n)
    error('W has wrong number of columns.');
end
if ~isempty(x0)
    if size(x0,1)~=n
        error('X0 and A must have the same number of rows.');
    end
    if size(x0,2)~=1
        error('X0 can only have one column.');
    end
else
    % Set the initial state to zero is x0 is empty
    x0 = zeros(n,1);
end

if isempty(w),
    w = zeros(N,n);
end

x = zeros(n,N);
for k = 1:N
    x(:,k) = x0;
    x0 = A(:,1:n)*x0 + B(:,1:m)*u(k,:)' + w(k,:)';
end
x = x';
 
 
 
 

