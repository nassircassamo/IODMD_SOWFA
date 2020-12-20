function C  =  cholicm(Af,Ab,Sf,Sb,N) 
%CHOLICM    This function calculates a Cholesky factor of 
%           the Inverse Covariance Matrix of a multivariable 
%           autoregressive noise process. 
% 
% Syntax: 
%           C = cholicm(Af,Ab,Sf,Sb,N); 
% 
% Description: 
%           The Inverse Covariance Matrix S of a multivariable 
%           autoregressive noise process according to [1] is 
%           calculated. The Cholesky factor C is returned such 
%           that C'*C = S 
% 
%           The noise model contains a causal and an auticausal 
%           part, both of which describe the actual noise v(k). 
%           If e(k) is a Gaussian white innovation, the model is 
%           given by: 
% 
%           v(k) = e_f(k) - Af_1 v(k-1) - ... - Af_d v(k-d) 
%           v(k) = e_b(k) - Ab_1 v(k+1) - ... - Ab_d v(k+d) 
% 
%           Subscripts f denote the causal (Forward) components while 
%           subscripts d denote the anti-causal (Backward) ones. 
% 
% Inputs: 
%  Af       A l x ld matrix containing the causal part of the noise 
%           process. Af = [Af1 Af2 ... Afd]. 
%  Ab       A l x ld matrix containing the anti-causal part of the noise 
%           process. Ab = [Ab1 Ab2 ... Abd]. 
%  Sf       A l x l matrix describing the covariance E[e_f e_f'] 
%  Sb       A l x l matrix describing the covariance E[e_b e_b'] 
%  N        The number of samples. 
% 
% Outputs: 
%  C        The Cholesky factor of the ICM. This matrix is stored in 
%           LAPACK/BLAS band-storage; its size is (d+1)*l x N, and the 
%           bottom row contains the diagonal of C. The row above contains 
%           a zero, and then the first superdiagonal of C. The row above 
%           contains two zeros, and then the second superdiagonal. Etc. 
%           The top row contains (d+1)*l-1 zeros, and then the 
%           ((d+1)*l-1)th superdiagonal. 
% 
% Limitations: 
%           A covariance matrix of a stationary process is always positive 
%           definite. However, it is very well possible to specify filter 
%           coefficients Af, Ab and covariances Sf and Sb such that the 
%           theoretical ICM calculated per [1] is not positive definite. 
%           In such cases, no Cholesky factor can be calculated, and an 
%           identity matrix will be returned along with a warning message. 
%           The filter should be checked and adjusted in these cases. 
% 
% References: 
%  [1]      Benoit David, 'Parameter Estimation in Nonlinear Dynamical 
%           Systems with Correlated Noise', Ph.D. thesis, Universite 
%           Catholique de Louvain, November 2001. 

% Revised by Ivo Houtzager, 2007
% Copyright (c) 1996-2007, Delft Center of Systems and Control 
 
% Check number of arguments
if nargin ~= 5
    error('CHOLICM requires five input arguments.');
end

% Extract dimensions
l = size(Af,1);
d = size(Af,2)/l;

% Check input parameters
if l < 1
    error('The number of variables cannot be zero.');
end
if round(d) ~= d
    error('Af should contain an integer number of AR matrices.');
end
if ~all(size(Af) == size(Ab))
    error('Af and Ab must be the same dimensions.');
end
if ~all(size(Sf) == [l l])
    error('The dimension of Sf should be the same as the number of parameters.');
end
if ~all(size(Sb) == [l l])
    error('The dimension of Sb should be the same as the number of parameters.');
end
if N < 2*d
    error('The number of samples should be at least 2d.');
end

% Preprocess the filter coefficient matrices in
% order to incorporate the block-diagonal Sf and Sb matrices
% to reduce the calculation of S to L^T*L - M^T*M
[Tf,pdflag] = chol(inv(Sf));
if pdflag~=0,
    warning('LTI:useIdentity','Causal innovation covariance matrix should be positive definite: using identity');
    Tf = eye(l);
end;
[Tb,pdflag] = chol(inv(Sb));
if pdflag~=0,
    warning('LTI:useIdentity','Anti-causal innovation covariance matrix should be positive definite: using identity');
    Tb = eye(l);
end;
Af = Tf*Af;
Ab = Tb*Ab;

% Transpose and flip Af and Ab
Af = [Tf; Af'];
Ab = Ab';
for i = 1:d
    % Tf (A0) has already been transposed, transpose the rest of Af
    Af(i*l+1:(i+1)*l,:) = Af(i*l+1:(i+1)*l,:)';

    % Transpose Ab
    Ab((i-1)*l+1:i*l,:) = flipud(Ab((i-1)*l+1:i*l,:)');
end
Ab = flipud(Ab);

% Allocate memory
S = sparse([],[],[],N*l,N*l,l*l*(N*(d+1)-d*(d+1)/2)-N*l*(l-1)/2);

% Calculate one block-row of the upper triangular part of L'*L
BRow = zeros(l,(d+1)*l);
for i = 0:d
    % Calculate the (i+1)th block in the block-row
    BRow(:,i*l+1:(i+1)*l) = Af(i*l+1:(d+1)*l,:)'*Af(1:(d+1-i)*l,:);
end;

% Make the (l-1)*(l-1) lower triangular part of BRow zero
BRow(:,1:l) = triu(BRow(:,1:l));

% Copy the data to the the block-rows the ICM
for i = 1:N-d
    % Copy BRow to block row i
    S((i-1)*l+1:i*l,(i-1)*l+1:(i+d)*l) = BRow;
end

% Calculate the lower-right corner the ICM.
BBlock = zeros(l*d,l*d);
for i = 1:d,
    BBlock((i-1)*l+1:d*l,(i-1)*l+1:i*l) = Af(1:(d+1-i)*l,:);
end;
S((N-d)*l+1:N*l,(N-d)*l+1:N*l) = triu(BBlock'*BBlock);

% Subtract V^T V from the ICM.
BBlock = zeros(l * d,l * d);
for i = 1:d,
    BBlock((i-1)*l+1:d*l,(i-1)*l+1:i*l) = Ab(1:(d+1-i)*l,:);
end;
S(1:l*d,1:l*d) = S(1:l*d,1:l*d)-triu(BBlock'*BBlock);

% Calculate the Cholesky factorization
[CC,pdflag] = chol(S);
if pdflag~=0,
    warning('LTI:useIdentity','Supplied noise model causes ICM to be non positive definite, using identity matrix');
    C = ones(1,N*l);
else
    % Return a LAPACK/BLAS style triangular band matrix
    C = spdiags(CC,(d+1)*l-1:-1:0)';
end

