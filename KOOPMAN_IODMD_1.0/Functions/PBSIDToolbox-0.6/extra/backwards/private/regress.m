function [X,reg_min,Z] = regress(Y,P,reg,opt,X0)
%REGRESS  Matrix regression with regularisation
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

N = size(Y,2);
if nargin <= 4
    if strcmpi(reg,'none')
        if nargout == 3
            Z = pinv(P);
            X = Y*Z;
        else
            X = Y*pinv(P);
        end
        reg_min = 0;
    elseif strcmpi(reg,'tikh')
        [U,S,V] = svd(P*P'./N);
        s = diag(S);
        YP = (Y*P'./N)';
        if isscalar(opt)
            reg_min = opt;
        elseif strcmpi(opt,'lcurve')
            reg_min = reglcurve(YP,U,s);
        elseif strcmpi(opt,'gcv')
            reg_min = reggcv(YP,U,s);
        end
        if nargout == 3
            Z = (V*(diag(s./(s.^2 + reg_min^2)))*U'*(P./N))';
            X = Y*Z;
        else
            X = (V*(diag(s./(s.^2 + reg_min^2)))*U'*YP)';
        end
    elseif strcmpi(reg,'tsvd')
        [U,S,V] = svd(P*P'./N);
        s = diag(S);
        YP = (Y*P'./N)';
        if isscalar(opt)
            k_min = opt;
        elseif strcmpi(opt,'lcurve')
            k_min = reglcurve(YP,U,s,'tsvd');
        elseif strcmpi(opt,'gcv')
            k_min = reggcv(YP,U,s,'tsvd');
        end
        if nargout == 3
            Z = (V(:,1:k_min)*diag(1./s(1:k_min))*U(:,1:k_min)'*(P./N))';
            X = Y*Z;
        else
            X = (V(:,1:k_min)*diag(1./s(1:k_min))*U(:,1:k_min)'*YP)';
        end
    end
else
    if strcmpi(reg,'none')
        X = X0 + (Y-X0*P)*pinv(P);
        reg_min = 0;
    elseif strcmpi(reg,'tikh')
        [U,S,V] = svd(P*P'./N);
        s = diag(S);
        YP = ((Y-X0*P)*P'./N)';
        if isscalar(opt)
            reg_min = opt;
        elseif strcmpi(opt,'lcurve')
            reg_min = reglcurve(YP,U,s);
        elseif strcmpi(opt,'gcv')
            reg_min = reggcv(YP,U,s);
        end
        if nargout == 3
            Z = (V*(diag(s./(s.^2 + reg_min^2)))*U'*(P./N))';
            X = Y*Z;
        else
            X = X0 + (V*(diag(s./(s.^2 + reg_min^2)))*U'*YP)';
        end
    elseif strcmpi(reg,'tsvd')
        [U,S,V] = svd(P*P'./N);
        s = diag(S);
        YP = ((Y-VARX0*P)*P'./N)';
        if isscalar(opt)
            k_min = opt;
        elseif strcmpi(opt,'lcurve')
            k_min = reglcurve(YP,U,s,'tsvd');
        elseif strcmpi(opt,'gcv')
            k_min = reggcv(YP,U,s,'tsvd');
        end
        if nargout == 3
            Z = (V(:,1:k_min)*diag(1./s(1:k_min))*U(:,1:k_min)'*(P./N))';
            X = X0 + Y*Z;
        else
            X = X0 + (V(:,1:k_min)*diag(1./s(1:k_min))*U(:,1:k_min)'*YP)';
        end
    end
end
reg_min = reg_min*N;
