function [x,resnorm,residual,exitflag,output,lambda,jacobian] = lmmore(fun,xinit,LowerBound,UpperBound,options,varargin)
%LMMORE    This function is a More-Hebden implementation of the 
%          Levenberg-Marquardt non-linear least-squares 
%          optimization algorithm. 
% 
% Syntax: 
%          x=lmmore('func',xinit,lb,ub,options,arg2,...) 
% 
%          [x,resnorm,residual,exitflag,output,lambda,jacobian] 
%              =lmmore('func',xinit,lb,ub,options,arg2,...) 
% 
% Inputs: 
%  'func'  The cost-function that is to be used. 
%  xinit   The parameter-vector's starting point in the non-linear 
%          optimization. 
%  lb      Lower-bound on the parameters 
%  ub      Upper-bound on the parameters 
%  options A MATLAB 6 compatible optimset-structure that contains 
%          options for the optimization algorithm. 
%  arg2    This will be passed as second argument to the cost-function 
%          'func'. 
%          Arguments 3 to N may be appended after arg2. 
% 
% Outputs: 
%  x       Result of the optimization. The solution x is guaranteed to 
%          have an equal or smaller cost than xinit. 
% 
%          All other parameters are compatible with MATLAB 6's lsqnonlin 
%          function. 
% 
% Remarks: 
%          The interface to lmmore has been made compatible with the lsqnonlin 
%          optimization function in the MATLAB 6 Optimization Toolbox. 
%          Note that although a lower and upper bound are given (consistent 
%          with lsqnonlin's interface), they are NOT used internally. 
% 
%          This optimization implementation supports overparametrized cost- 
%          functions. If options.Manifold (not part of optimset's normal 
%          structures) is passed and set to 'on', lmmore expects the cost- 
%          function to be able to return three arguments: an error-vector, 
%          a projected Jacobian and a projection matrix. The column-space 
%          of this  matrix form a basis orthogonal to the subspace in which 
%          the cost-function does not change due to the overparametrization. 
% 
%          This optimization implementation supports cost-functions that 
%          return the R-factor of the (projected) Jacobian and the error-vector: 
% 
%          [ J E ] = Q R 
% 
%          Cost-functions may use this functionality, e.g. to build up the 
%          R-factor in such a way that less memory is required. 
% 
% See also: lsqnonlin, optimset 
 
% Niek Bergboer, 2001 
% Revised by Ivo Houtzager, 2007
% Copyright (c) 2001-2007, Delft Center of Systems and Control  

% Check number of arguments
if nargin < 2
    error('LMMORE requires at least two input arguments.');
end

xinit = xinit(:); 
xsize = size(xinit,1); 
ytest = feval(fun,xinit,varargin{:}); 
ysize = size(ytest,2); 

if nargin > 3
    if ~isempty(UpperBound) || ~isempty(LowerBound)
        warning('LTI:ignoreBounds','Upper and lower parameter bounds are ignored by lmmore.');
    end
end

defaultopt  =  struct('Display','final', ...
    'LargeScale','off',  ...
    'TolX',1e-6,'TolFun',1e-6, ...
    'Jacobian','on', ...
    'MaxFunEvals',100 * xsize, ...
    'TypicalX',ones(xsize,1), ...
    'MaxIter',400, ...
    'LevenbergMarquardt','on');

MatlabVersion = version;
if MatlabVersion(1)<'6',
    % Use MATLAB 5 syntax
    curroptions = defaultopt;
    if isfield(options,'Display') && ~isempty(options.Display)
        curroptions.Display = options.Display;
    end
    if isfield(options,'LargeScale') && ~isempty(options.LargeScale)
        curroptions.LargeScale = options.LargeScale;
    end
    if isfield(options,'TolX') && ~isempty(options.TolX)
        curroptions.TolX = options.TolX;
    end
    if isfield(options,'TolFun') && ~isempty(options.TolFun)
        curroptions.TolFun = options.TolFun;
    end
    if isfield(options,'Jacobian') && ~isempty(options.Jacobian)
        curroptions.Jacobian = options.Jacobian;
    end
    if isfield(options,'MaxFunEvals') && ~isempty(options.MaxFunEvals)
        curroptions.MaxFunEvals = options.MaxFunEvals;
    end
    if isfield(options,'TypicalX') && ~isempty(options.TypicalX)
        curroptions.TypicalX = options.TypicalX;
    end
    if isfield(options,'MaxIter') && ~isempty(options.MaxIter)
        curroptions.MaxIter = options.MaxIter;
    end
    if isfield(options,'LevenbergMarquardt') && ~isempty(options.LevenbergMarquardt)
        curroptions.LevenbergMarquardt = options.LevenbergMarquardt;
    end
else
    % Use MATLAB 6 syntax
    curroptions = optimset(defaultopt,options);
end

if ~isfield(options,'Manifold')
    curroptions.Manifold = 'off';
else
    curroptions.Manifold = options.Manifold;
end
if ~isfield(options,'RFactor')
    curroptions.RFactor = 'off';
else
    curroptions.RFactor = options.RFactor;
end

% Extract parameters
TolX = curroptions.TolX;
TolFun = curroptions.TolFun;
MaxIter = curroptions.MaxIter;
MaxFunEvals = curroptions.MaxFunEvals;

if ~strcmp(curroptions.Jacobian,'on')
    error('This Levenberg-Marquardt algorithm requires a Jacobian');
end
% if ~strcmp(curroptions.LevenbergMarquardt,'on')
%     error('Cannot do Levenberg-Marquardt if OPTIONS.LevenbergMarquardt=''off''');
% end
if strcmp(curroptions.Display,'off'),
    verbosity = 0;
elseif strcmp(curroptions.Display,'final'),
    verbosity = 1;
elseif strcmp(curroptions.Display,'iter') || strcmp(curroptions.Display,'notify'),
    verbosity = 2;
elseif strcmp(curroptions.Display,'testing'),
    verbosity = Inf;
else
    verbosity = 0;
end

ready = 0;
iter = 0;
fCount = 0;
gCount = 0;
firsttime = 1;
xcurr = xinit;
Deltak = 1e-1;
% Initial size of trust region
Deltakold = Deltak;
maxDeltak = 1e10;

if verbosity >= 2
    if strcmp(curroptions.Manifold,'on')
        fprintf('Iter  fCount      Residual     norm(SD)       Lambda Trust-region Linearity  NFP\n');
    else
        fprintf('Iter  fCount      Residual     norm(SD)       Lambda Trust-region Linearity\n');
    end
end

while ~ready
    iter = iter+1;
    %xold = xcurr;

    thisiterok = 0;
    %subiter = 0;

    if firsttime
        if strcmp(curroptions.Manifold,'on')
            [En,PsiN,U2] = feval(fun,xcurr,varargin{:});
        else
            [En,PsiN] = feval(fun,xcurr,varargin{:});
        end
        fCount = fCount+1;
        gCount = gCount+1;
        costcurr = (0.5/ysize)*(En'*En);

        % Set variables for all subsequent iterations
        lambda = 0;
        firsttime = 0;
    end
    if strcmp(curroptions.Manifold,'on')
        % Check whether the cost-function returns PsiN*U2, or
        % only the R-factor from the QR-factorization of PsiN*U2
        if strcmp(curroptions.RFactor,'on')
            R11 = PsiN(1:size(U2,2),1:size(U2,2));
            R12 = PsiN(1:size(U2,2),size(U2,2)+1);
            clear PsiN;
        else
            M = [PsiN En];
            clear PsiN;
            R = triu(qr(M,0));
            clear M;
            R11 = R(1:size(U2,2),1:size(U2,2));
            R12 = R(1:size(U2,2),size(U2,2)+1);
            clear R;
        end
    else
        % Check whether the cost-function returns PsiN*U2, or
        % only the R-factor from the QR-factorization of PsiN*U2
        if strcmp(curroptions.RFactor,'on')
            R11 = PsiN(1:xsize,1:xsize);
            R12 = PsiN(1:xsize,xsize+1);
            clear PsiN;
        else
            M = [PsiN En];
            clear PsiN;
            R = triu(qr(M,0));
            clear M;
            R11 = R(1:xsize,1:xsize);
            R12 = R(1:xsize,xsize+1);
            clear R;
        end
    end
    [Un,Sn,Vn] = svd(R11);
    sn = diag(Sn);
    UnTR12 = Un'*R12;

    while ~thisiterok
        if strcmp(curroptions.Manifold,'on')
            dnew = -U2*Vn*diag(sn.^(-1))*UnTR12;
        else
            dnew = -Vn*diag(sn.^(-1))*UnTR12;
        end

        if norm(dnew)>Deltak || (sn(1) >1e8  * sn(length(sn)))
            % No; larger than trust region or (projected) Jacobian
            % is too ill-conditioned.
            % Initialize Inner Iteration
            piter = 0;
            pready = 0;

            % Initialize vairables
            ub = norm(R11'*R12,2)/Deltak;

            % Need the singular values of R11 for now...
            if (sn(1)>1e8*sn(length(sn))),
                lb = 0;
            else
                % Now, dnew is the Gauss-Newton step, so indeed phi(0)
                lb = (norm(dnew,2)-Deltak);
                lb = lb*norm(dnew)/(UnTR12' * diag(sn.^(-4))*UnTR12);
            end;

            if lambda == 0,
                lambda = eps;
            else
                lambda = lambda-((phi+Deltak)/(Deltak))*(phi/phia);
            end

            while ~pready
                lambda = max(lambda,lb);
                lambda = min(lambda,ub);

                if strcmp(curroptions.Manifold,'on')
                    dnew = -U2*Vn*diag(sn./(sn.^2+lambda))*UnTR12;
                else
                    dnew = -Vn*diag(sn./(sn.^2+lambda))*UnTR12;
                end
                phi = norm(dnew,2)-Deltak;
                phia = -(UnTR12'*diag((sn.^2)./((sn.^2 + lambda).^3))*UnTR12)/norm(dnew,2);

                % Check termination
                if abs(phi)<0.1*Deltak
                    % Inner Iteration is ready
                    pready = 1;
                else
                    if phi<0,
                        ub = lambda;
                    end;
                    lb = max(lb,lambda-phi/phia);
                    lambda = lambda-((phi+Deltak)/(Deltak))*(phi/phia);
                end
                if verbosity == Inf
                    fprintf('piter: %d, %8g, %8g, %8g, %8g\n',piter,norm(dnew)/Deltak,lb,ub,lambda);
                end

                % Check Inner Convergence
                if ub==0 && lb==0
                    warning('LTI:converge','Inner-iteration for lambda did not converge.');
                    pready = 1;
                end
            end
        else
            % Gauss-Newton step (no regularization)
            lambda = 0;
        end
        xnew = xcurr+dnew;
        Ennew = feval(fun,xnew,varargin{:});
        fCount = fCount+1;
        costnew = (0.5/ysize)*(Ennew'*Ennew);

        % Predict the new costs
        dTPsiNTEn = UnTR12'*diag((sn.^2)./(sn.^2+lambda))*UnTR12;
        costpred = costcurr-(1/ysize)*dTPsiNTEn+(0.5/ysize)*UnTR12'*diag((sn.^4)./((sn.^2+lambda).^2))*UnTR12;

        if (costpred-costcurr) ~=0
            r = (costnew-costcurr)/(costpred-costcurr);
        else
            r = 1e20;
        end

        if costnew<=costcurr,
            xcurr = xcurr+dnew;
            %d = dnew;
            thisiterok = 1;
        end
        Deltakold = Deltak;
        if r < 0.25,
            % Poor linearity
            if strcmp(curroptions.Manifold,'on')
                lfactor = 2-(costnew-costcurr)/((0.5/ysize)*dnew'*U2*R11'*R12);
            else
                lfactor = 2-(costnew-costcurr)/((0.5/ysize)*dnew'*R11'*R12);
            end
            lfactor = min(lfactor,10);
            lfactor = max(lfactor,2);
            Deltak = Deltak/lfactor;
        elseif r > 0.75
            % Good linearity
            Deltak = Deltak*2;
            Deltak = min(Deltak,maxDeltak);
        elseif (r > 0.25 && r < 0.75) && lambda==0
            % Reasonable linearity and Gauss-Newton step
            Deltak = Deltak*2;
            Deltak = min(Deltak,maxDeltak);
        end
    end
    
    % Calculate New Error and Jacobian
    if strcmp(curroptions.Manifold,'on')
        [En,PsiN,U2] = feval(fun,xcurr,varargin{:});
    else
        [En,PsiN] = feval(fun,xcurr,varargin{:});
    end

    % Update function and gradient evaluation counts
    fCount = fCount+1;
    gCount = gCount+1;

    % Print information
    if verbosity>=2
        if strcmp(curroptions.Manifold,'on'),
            fprintf('%4d %3d,%3d  %12g %12g %12g %12g  %8g %4d\n',iter,fCount,gCount,Ennew' * Ennew,norm(dnew),lambda,Deltakold,r,size(U2,2));
        else
            fprintf('%4d %3d,%3d  %12g %12g %12g %12g  %8g\n',iter,fCount,gCount,Ennew' * Ennew,norm(dnew),lambda,Deltakold,r);
        end
    end
    
    if max(abs(dnew)) < TolX
        if verbosity>0,
            fprintf('Optimization terminated successfully:\n Search direction less than tolX\n');
        end;
        ready = 1;
        exitflag = 1;
    end

    if strcmp(curroptions.RFactor,'on')
        if strcmp(curroptions.Manifold,'on')
            EnTPsiN = PsiN(1:size(U2,2),size(U2,2)+1)'*PsiN(1:size(U2,2),1:size(U2,2))*U2';
        else
            EnTPsiN = PsiN(1:xsize,xsize+1)'*PsiN(1:xsize,1:xsize);
        end
    else
        if strcmp(curroptions.Manifold,'on')
            EnTPsiN = (En'*PsiN)*U2';
        else
            EnTPsiN = En'*PsiN;
        end
    end

    if abs(EnTPsiN*dnew) < TolFun && max(abs(EnTPsiN)) < 10*(TolX+TolFun)
        if verbosity > 0
            fprintf('Optimization terminated successfully:\n Gradient in the search direction less than tolFun\n');
        end
        ready = 1;
        exitflag = 1;
    end

    if iter >= MaxIter
        if verbosity > 0
            fprintf('Maximum number of iterations exceeded\n Increase OPTIONS.MaxIter\n');
        end
        ready = 1;
        exitflag = 0;
    end

    if fCount >= MaxFunEvals
        if verbosity > 0
            fprintf('Maximum number of function evaluations exceeded\n Increase OPTIONS.MaxFunEvals\n');
        end
        ready = 1;
        exitflag = 0;
    end

    if abs(costnew-costcurr) < TolFun*abs(costcurr)
        if verbosity > 0
            fprintf('Optimization terminated:\n Relative decrease of cost-function less than TolFun\n');
        end
        ready = 1;
        exitflag = 1;
    end
    costcurr = costnew;
end
x = xcurr;
resnorm = costcurr;
residual = En;
output = struct('iterations',iter,'funcCount',fCount, ...
    'algorithm','Levenberg-Marquardt More-Hebden trust region');
lambda = struct('lower',zeros(xsize,1),'upper',zeros(xsize,1));

if strcmp(curroptions.RFactor,'off')
    jacobian = PsiN;
else
    jacobian = [];
end




