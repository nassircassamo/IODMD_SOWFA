function [E,covE] = dvar2eig(P,A)
%DVAR2EIG Eigenvalues and its covariance estimation
%  [E,covE]=dvar2frd(P,A) returns the estimated eigenvalues and its
%  covariance for the state space matrix A and its covariance P.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

n = size(A,1);
P = P(1:n^2,1:n^2);
J = jacobianest(@(x) eigen(x,n),A(:));
covE = J*P*J';
E = eig(A);   
end

function dE = eigen(A,n)
A = reshape(A,n,n);
E = eig(A);
dE = zeros(2*n,1);
dE(1:2:end) = real(E);
dE(2:2:end) = imag(E);
end    
  