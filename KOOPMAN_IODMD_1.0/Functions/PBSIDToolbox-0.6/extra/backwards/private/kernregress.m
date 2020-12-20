function A = kernregress(Y,Z,reg,opt)
%REGRESS  Matrix regression with regularisation
%  X = REGRESS(Y,P) solves the least squares regression problem. The inputs
%  are the matrices Y and P, which have the same number of columns. Output
%  is the least squares estimate X.
%
%  X = REGRESS(Y,P,REG,OPT) solves the regularized least squares problem.
%  The additional inputs are the regularization method and selection
%  parameters: REG = {'NONE', 'TIKH', TSVD'} and OPT = {'GCV', 'LCURVE'}.
%  With regularisation, the solver can better deal with singular covariance
%  matrices. (default REG='NONE' and OPT='GCV')
%
%  See also LSCCA.

%  Ivo Houtzager
%
%  Delft Center of Systems and Control
%  The Netherlands, 2009

if nargin < 4 || isempty(opt)
    opt = 'gcv';
end
if nargin < 3 || isempty(reg)
    reg = 'none';
end

% center the variables
%X = X - repmat(mean(X,1), rx, 1);
%Y = Y - repmat(mean(Y,1), ry, 1);

if strcmpi(reg,'none')
    A = Y*pinv(Z);
elseif strcmpi(reg,'tikh')
    [V,s] = eig(Z);
    s = abs(diag(s));
    [s,I] = sort(s,'descend');
    V = real(V(:,I));

    if strcmpi(opt,'lcurve')
        reg_min = reglcurve(Y',V,s);
    elseif strcmpi(opt,'gcv')
        reg_min = reggcv(Y',V,s);
    end
    A = (V*(diag(s./(s.^2 + reg_min^2)))*V'*Y')';
elseif strcmpi(reg,'tsvd')
    [V,s] = eig(Z);
    s = abs(diag(s));
    [s,I] = sort(s,'descend');
    V = real(V(:,I));
    
    if strcmpi(opt,'lcurve')
        k_min = reglcurve(Y',V,s,'tsvd');
    elseif strcmpi(opt,'gcv')
        k_min = reggcv(Y',V,s,'tsvd');
    end
    A = (V(:,1:k_min)*diag(1./s(1:k_min))*V(:,1:k_min)'*Y')';
end


