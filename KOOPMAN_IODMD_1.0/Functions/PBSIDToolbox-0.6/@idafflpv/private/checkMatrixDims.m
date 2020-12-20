function [Ns,Ny,Nu,Np] = checkMatrixDims(a,b,c,d,k,x0)
% Determines and checks the system matrix sizes.

% Check sizes
if size(a,2) ~= size(c,2)
    error('Number of columns of A and C must be equal.');
end
if size(a,1) ~= size(b,1)
    error('Number of rows of A and B must be equal.');
end
if size(a,1) ~= size(k,1)
    error('Number of rows of A and K must be equal.');
end
if size(a,1) ~= size(x0,1)
    error('Number of rows of A and x0 must be equal.');
end
if size(b,2) ~= size(d,2)
    error('Number of columns of B and D must be equal.');
end
if size(c,1) ~= size(d,1)
    error('Number of rows of C and D must be equal.');
end

% Number of states
Ns = size(a,1);

% Determine number scheduling parameters
Np = size(a,2)/Ns;

% Number of inputs
Nu = size(b,2)/Np;

% Number of outputs
Ny = size(c,1);

Np = Np - 1;

