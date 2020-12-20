function [X,reg_min,Z] = regress(Y,P,reg,opt,X0)
%REGRESS  Matrix regression with regularisation
%  X = REGRESS(Y,P) solves the least squares regression problem. The inputs
%  are the matrices Y and P, which have the same number of columns. Output
%  is the least squares estimate X.
%
%  X = REGRESS(Y,P,REG,OPT) solves the regularized least squares problem.
%  The additional inputs are the regularization method (REG) and selection
%  parameters (OPT):
%
%  - If REG = {'TIKH', TSVD'}, i.e. Tikhonov regularization or  Truncated 
%    Singular Value Decomposition, we can choose 
%    OPT = {'GCV', 'LCURVE', or any scalar value}, i.e. Generalized Cross 
%    Validation or L-Curve regularization.
%    With regularisation, the solver can better deal with singular 
%    covariance matrices. (default REG='NONE' and OPT='GCV') [1]
%
%  - If REG = {'BPDN'} then sparse estimation through Basis Pursuit 
%    Denoising is used. The solver can then better deal with a past window 
%    that is chosen (too) large. Choose OPT = {'SV', or any scalar value} 
%    to select trade-off parameter reg_min in the BPDN problem:
%       min ||X||_1 s.t. ||Y-X*P||_F < E0*(1+reg_min)
%    where E0 is the residual norm ||Y-XP||_F of the regression problem
%    solved without regularization. (default OPT = 'SV')
%    For OPT = {'AIC'} the trade-off point is chosen as:
%       reg_min = 2*size(P,1)/size(P,2)
%    following Aikaike Information Criterion, see [2,3] for more 
%    explanation.
%    For OPT = {'SV'} the trade-off point is chosen on the basis of the 
%    prediction error on portion of the data not used for the regression 
%    (i.e. validation data). reg_min is stepwise decreased until an
%    increase in the error on validation data (i.e. overfitting) is 
%    detected. [2]
%    Also, we can choose a structure array for OPT with:
%    * OPT.opt = {'SV', or any scalar value}
%    * OPT.verbosity = {1,2} to display iterations of the BPDN solver SPGL1
%      (default: verbosity = 0)
%    * OPT.valportion is a scalar value between 0 and 1 setting the portion 
%      of the data used as validation data when OPT.opt = 'SV'.
%      (default: valportion = 1/4) 
%
%  See also KERNREGRESS, REGLCURVE, REGGCV.

%  References:
%    [1] ``Regularization Tools'' by P.C. Hansen, 
%        http://www2.imm.dtu.dk/~pch/Regutools/
%    [2] P.M.O. Gebraad, J.W. van Wingerden and M. Verhaegen,
%        ``Sparse Estimation for Predictor-Based Subspace Identification of
%        LPV Systems'' submitted for SYSID 2012, 16th IFAC Symposium on 
%        System Identification, Brussels, Belgium, July 11-13, 2012
%    [3] C. Rojas, and H. Hjalmarsson,
%        ``Sparse estimation based on a validation criterion'', proc. pf
%        the 50th IEEE Conf. on Decision and Control (CDC). Orlando,
%        Florida, USA.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  The Netherlands, 2009

%  Pieter Gebraad
%  Delft Center of Systems and Control
%  The Netherlands, 2011

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
            s = warning('query', 'MATLAB:lscov:RankDefDesignMat');
            warning off MATLAB:lscov:RankDefDesignMat
            warning off regress:lastbeforelscov
            warning('regress:lastbeforelscov','Last Before')
            X = lscov(P',Y')'; % edit by Pieter *, used to be X = Y*pinv(P);
            [msgstr, msgid] = lastwarn;
            if strcmp(msgid,'MATLAB:lscov:RankDefDesignMat')
            	warning('regress:RankDefDataMat', 'The data matrix P in the regression problem solved by REGRESS is rank deficient, the solution may be inaccurate.\n Consider using a regularization method.\n See HELP for REGRESS or for the identification function you used, i.e. LORDVARX/DORDVARX/etc.')
            end
            warning(s.state,'MATLAB:lscov:RankDefDesignMat');
            warning on regress:lastbeforelscov
        end
        reg_min = 0;
	elseif strcmpi(reg,'nuclear') || strcmpi(reg,'nuc')
        if isfield(opt,'opt')
           opt = opt.opt;
        end
        X = nuclear(Y,P,opt);
		reg_min = opt;
    elseif strcmpi(reg,'tikh')
        if isfield(opt,'opt')
           opt = opt.opt;
        end
        [V,S] = eig(P*P'./N);
        S = abs(diag(S));
        [S,I] = sort(S,'descend');
        V = real(V(:,I));     
        YP = (Y*P'./N)';
        if isscalar(opt)
            reg_min = opt;
        elseif strcmpi(opt,'lcurve')
            reg_min = reglcurve(YP,V,S);
        elseif strcmpi(opt,'gcv')
            reg_min = reggcv(YP,V,S);
        end
        if nargout == 3
            Z = (V*(diag(S./(S.^2 + reg_min^2)))*V'*(P./N))';
            X = Y*Z;
        else
            X = (V*(diag(S./(S.^2 + reg_min^2)))*V'*YP)';
        end
		reg_min = reg_min*N;
    elseif strcmpi(reg,'tsvd')
        if isfield(opt,'opt')
           opt = opt.opt;
        end
        [V,S] = eig(P*P'./N);
        S = abs(diag(S));
        [S,I] = sort(S,'descend');
        V = real(V(:,I));
        YP = (Y*P'./N)';
        if isscalar(opt)
            k_min = opt;
        elseif strcmpi(opt,'lcurve')
            k_min = reglcurve(YP,V,S,'tsvd');
        elseif strcmpi(opt,'gcv')
            k_min = reggcv(YP,V,S,'tsvd');
        end
        if nargout == 3
            Z = (V(:,1:k_min)*diag(1./S(1:k_min))*V(:,1:k_min)'*(P./N))';
            X = Y*Z;
        else
            X = (V(:,1:k_min)*diag(1./S(1:k_min))*V(:,1:k_min)'*YP)';
        end
		reg_min = k_min;
    elseif strcmpi(reg,'bpdn')
        if isfield(opt,'verbosity')
           options = spgSetParms('verbosity',opt.verbosity);
        else
           options = spgSetParms('verbosity',0);
        end
        if isfield(opt,'valportion')
            if isscalar(opt.valportion) && opt.valportion < 1 && opt.valportion > 0
                valportion = opt.valportion;
            else
                error('In regress.m, valportion should be chosen as a scalar between 0 and 1.')
            end
        else
            valportion = 0.25;
        end
        if isfield(opt,'opt')
           opt = opt.opt;
        end
        if isscalar(opt)
            X = regress(Y,P);
            E = norm(Y-X*P,'fro');
            eps = E*(1+opt);
            reg_min = opt;
            X = spg_mmv(P',Y',eps,options)';
        elseif strcmpi(opt,'aic')
            X = regress(Y,P);
            E = norm(Y-X*P,'fro');
            reg_min = 2*size(P,1)/size(Y,2);
            eps = E*(1+reg_min);
            X = spg_mmv(P',Y',eps,options)';
        elseif strcmpi(opt,'sv')
            N = size(Y,2);
            Nid = ceil((1-valportion)*N);
            X = spg_mmv_stopvali(P(:,1:Nid)',Y(:,1:Nid)',P(:,Nid+1:N)',Y(:,Nid+1:N)',0,options)';
            if nargout>1
                X0 = regress(Y,P);
                E0 = norm(Y-X0*P,'fro');
                E = norm(Y-X*P,'fro');
                reg_min = E/E0-1;
            end
        end
    else
        error(['Regularization method ''', reg,''' not available.']);
    end
else
    if strcmpi(reg,'none')
        X = X0 + lscov(P',(Y-X0*P)')'; % edit by Pieter *, used to be X = X0 + (Y-X0*P)*pinv(P);  
        reg_min = 0;
	elseif strcmpi(reg,'nuclear') || strcmpi(reg,'nuc')	
		X = nuclear(Y,P,opt,[],X0);
		reg_min = opt;
    elseif strcmpi(reg,'tikh')
        [V,S] = eig(P*P'./N);
        S = abs(diag(S));
        [S,I] = sort(S,'descend');
        V = real(V(:,I));
        YP = ((Y-X0*P)*P'./N)';
        if isscalar(opt)
            reg_min = opt;
        elseif strcmpi(opt,'lcurve')
            reg_min = reglcurve(YP,V,S);
        elseif strcmpi(opt,'gcv')
            reg_min = reggcv(YP,V,S);
        end
        if nargout == 3
            Z = (V*(diag(S./(S.^2 + reg_min^2)))*V'*(P./N))';
            X = Y*Z;
        else
            X = X0 + (V*(diag(S./(S.^2 + reg_min^2)))*V'*YP)';
        end
		reg_min = reg_min*N;
    elseif strcmpi(reg,'tsvd')
        [V,S] = eig(P*P'./N);
        S = abs(diag(S));
        [S,I] = sort(S,'descend');
        V = real(V(:,I));
        YP = ((Y-X0*P)*P'./N)';
        if isscalar(opt)
            k_min = opt;
        elseif strcmpi(opt,'lcurve')
            k_min = reglcurve(YP,V,S,'tsvd');
        elseif strcmpi(opt,'gcv')
            k_min = reggcv(YP,V,S,'tsvd');
        end
        if nargout == 3
            Z = (V(:,1:k_min)*diag(1./S(1:k_min))*V(:,1:k_min)'*(P./N))';
            X = X0 + Y*Z;
        else
            X = X0 + (V(:,1:k_min)*diag(1./S(1:k_min))*V(:,1:k_min)'*YP)';
        end
		reg_min = k_min;
    else
        error(['Regularization method ''', reg, ''' not available for batchwise regression.'])
    end
end
end
% * LSCOV (uses QR) is faster than PINV (uses SVD), and finds 
%   basic solution instead of minimum norm solution, which may better
%   represent the (decaying) structure of the VARX parameter