function [A,reg_min] = kernregress(Y,Z,reg,opt)
%KERNREGRESS  Matrix regression with regularisation
%  X = REGRESS(Y,P) solves the least squares regression problem. The inputs
%  are the matrices Y and P, which have the same number of columns. Output
%  is the least squares estimate X.
%
%  X = REGRESS(Y,P,REG,OPT) solves the regularized least squares problem.
%  The additional inputs are the regularization method and selection
%  parameters: REG = {'NONE', 'TIKH', TSVD'} and OPT = {'GCV', 'LCURVE', or
%  any regularisation value as scalar}. With regularisation, the solver can
%  better deal with singular covariance matrices. (default REG='NONE' and
%  OPT='GCV')
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

if strcmpi(reg,'none')
    s = warning('query', 'MATLAB:lscov:RankDefDesignMat');
    warning off MATLAB:lscov:RankDefDesignMat
    warning off regress:lastbeforelscov
    warning('regress:lastbeforelscov','Last Before')
    A = lscov(Z',Y')';
    [msgstr, msgid] = lastwarn;
    if strcmp(msgid,'MATLAB:lscov:RankDefDesignMat')
        warning('regress:RankDefDataMat', 'The data matrix Z in the regression problem solved by KERNREGRESS is rank deficient, the solution may be inaccurate.\n Consider using a regularization method. See HELP for KERNREGRESS or for the identification function you used, i.e. LORDVARX/DORDVARX/etc.')
    end
    warning(s.state,'MATLAB:lscov:RankDefDesignMat');
    warning on regress:lastbeforelscov
    reg_min = 0;
elseif strcmpi(reg,'nuclear') || strcmpi(reg,'nuc')
    A = nuclear(Y,Z,opt);  
    reg_min = opt;
elseif strcmpi(reg,'tikh')
    [V,s] = eig(Z);
    s = abs(diag(s));
    [s,I] = sort(s,'descend');
    V = real(V(:,I));

    if isscalar(opt)
        reg_min = opt;
    elseif strcmpi(opt,'lcurve')
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
    
    if isscalar(opt)
        k_min = opt;
    elseif strcmpi(opt,'lcurve')
        k_min = reglcurve(Y',V,s,'tsvd');
    elseif strcmpi(opt,'gcv')
        k_min = reggcv(Y',V,s,'tsvd');
    end
    A = (V(:,1:k_min)*diag(1./s(1:k_min))*V(:,1:k_min)'*Y')';
end


