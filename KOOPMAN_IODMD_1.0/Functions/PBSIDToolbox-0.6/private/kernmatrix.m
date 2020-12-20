function omega = kernmatrix(Xtrain,kernel_type,kernel_pars,Xt)
%KERNMATRIX Construct the positive (semi-) definite and symmetric kernel matrix
%
%  Omega = kernel_matrix(X, kernel_fct, sig2)
%
%  This matrix should be positive definite if the kernel function
%  satisfies the Mercer condition. Construct the kernel values for
%  all test data points in the rows of Xt, relative to the points of X.
%
%  Omega_Xt = kernel_matrix(X, kernel_fct, sig2, Xt)
%
%
%  Full syntax:
%
%  Omega = kernel_matrix(X, kernel_fct, sig2)
%  Omega = kernel_matrix(X, kernel_fct, sig2, Xt)
%
%  Outputs:
%   Omega  : N x N (N x Nt) kernel matrix
%  Inputs:
%   X      : N x d matrix with the inputs of the training data
%   kernel : Kernel type (by default 'RBF_kernel')
%   sig2   : Kernel parameter (bandwidth in the case of the 'RBF_kernel')
%   Xt(*)  : Nt x d matrix with the inputs of the test data

% Copyright (c) 2002,  KULeuven-ESAT-SCD,
% License & help @ http://www.esat.kuleuven.ac.be/sista/lssvmlab

nb_data = size(Xtrain,1);

if nb_data> 3000,
    error('Too memory intensive, the kernel matrix is restricted to size 3000 x 3000 ');
end

if strcmpi(kernel_type,'rbf'),
    if nargin<4,
        XXh = sum(Xtrain.^2,2)*ones(1,nb_data);
        omega = (XXh+XXh') - 2*(Xtrain*Xtrain');
        omega = exp(-omega./kernel_pars(1));
    else
        XXh1 = sum(Xtrain.^2,2)*ones(1,size(Xt,1));
        XXh2 = sum(Xt.^2,2)*ones(1,nb_data);
        omega = XXh1+XXh2' - 2*Xtrain*Xt';
        omega = exp(-omega./kernel_pars(1));
    end
else
    if nargin<4,
        omega = zeros(nb_data,nb_data);
        for i=1:nb_data,
            omega(i:end,i) = feval(lower(kernel_type),Xtrain(i,:),Xtrain(i:end,:),kernel_pars);
            omega(i,i:end) = omega(i:end,i)';
        end
    else
        if size(Xt,2)~=size(Xtrain,2),
            error('dimension test data not equal to dimension traindata;');
        end
        omega = zeros(nb_data, size(Xt,1));
        for i=1:size(Xt,1),
            omega(:,i) = feval(lower(kernel_type),Xt(i,:),Xtrain,kernel_pars);
        end
    end
end
end

function x = lin(a,b,c)
% kernel function for implicit higher dimension mapping, based on
% the standard inner-product
%
%   x = lin_kernel(a,b)
%
% 'a' can only contain one datapoint in a row, 'b' can contain N
% datapoints of the same dimension as 'a'.
%
% see also:
%    poly_kernel, RBF_kernel, MLP_kernel, trainlssvm, simlssvm

% Copyright (c) 2002,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.ac.be/sista/lssvmlab


x = zeros(size(b,1),1);
for i=1:size(b,1),
    x(i,1) = a*b(i,:)';
end
end

function x = poly(a,b,d)
% polynomial kernel function for implicit higher dimension mapping
%
%  X = poly_kernel(a,b,[t,degree])
%
% 'a' can only contain one datapoint in a row, 'b' can contain N
% datapoints of the same dimension as 'a'.
%
% x = (a*b'+t^2).^degree;
%
% see also:
%    RBF_kernel, lin_kernel, MLP_kernel, trainlssvm, simlssvm
%

% Copyright (c) 2002,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.ac.be/sista/lssvmlab

if length(d)>1, d=d(2); t=d(1); else d = d(1);t=1; end
d = (abs(d)>=1)*abs(d)+(abs(d)<1); % >=1 !!

x = zeros(size(b,1),1);
for i=1:size(b,1),
    x(i,1) = (a*b(i,:)'+t^2).^d;
end
end

function x = mlp(a,b, par)
% Multi Layer Perceptron kernel function for implicit higher dimension mapping
%
%   x = MLP_kernel(a,b,[s,t])
%
% 'a' can only contain one datapoint in a row, 'b' can contain N
% datapoints of the same dimension as 'a'.
%
%   x = tanh(s*a'b+t^2)
%
% see also:
%    poly_kernel, lin_kernel, RBF_kernel, trainlssvm, simlssvm

% Copyright (c) 2002,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.ac.be/sista/lssvmlab



if length(par)==1, par(2) = 1; end
x = zeros(size(b,1),1);
for i=1:size(b,1),
    dp = a*b(i,:)';
    x(i,1) = tanh(par(1)*dp + par(2)^2);
end
end