function reg_c=reglcurve(Y,Vn,Sn,method,show)
%REGLCURVE   Compute regularization using L-curve criterion.
%            Determine the regularization parameter for ordkernel
%            using  L-curve criterion. It plots the L-curve and
%            find its corner. If the regularization method is
%            'tsvd' then the Spline Toolbox is needed to determine
%            the corner. If this toolbox is not available NaN is
%            returned.
%
% Syntax:
%            reg=reglcurve(Y,V,S)
%            reg=reglcurve(Y,V,S,method,show)
%
% Input:
%            Y,V,S    Data matrices from lpvkernel or bilkernel.
%            method   Regularization method to be used.
%                     'Tikh' - Tikhonov regularization (default).
%                     'tsvd' - Truncated singular value decomposition.
%            show     Display intermediate steps of the algorithm.
%
% Output:
%            reg      Regularization paramater for the kernel
%                     subspace identification method of ordkernel.

% Written by Vincent Verdult, May 2004.
% Based on Regularization Tools by P. C. Hansen

% default method
if nargin<4
    method='Tikh';
end
if nargin<5
    show=0;
end
if size(Y,1)~=size(Vn,1)
    error('The number of rows in Y must equal the number of rows in V.')
end
if size(Vn,1)~=size(Vn,2)
    error('V must be a square matrix.')
end
if size(Sn,2)~=1
    error('S must be a column vector.')
end
if size(Vn,1)~=size(Sn,1)
    error('The number of rows in S must equal the number of rows in V.')
end

% Initialization.
N = size(Sn,1);
s=sqrt(Sn);
beta = Vn'*Y;
xi = diag(1./s)*beta;

%%%%%%%%%%%%%%%%

% Tikhonov regularization
if (strncmp(method,'Tikh',4) || strncmp(method,'tikh',4))
    
    SkipCorner=0;
    txt = 'Tikh.';
    marker='-';
    npoints = 200;        % Number of points on the curve.
    smin_ratio = 16*eps;  % Smallest regularization parameter.
    eta = zeros(npoints,1);
    rho = zeros(npoints,1);
    reg_param = zeros(npoints,1);
    reg_param(npoints) = max([s(N),s(1)*smin_ratio]);
    ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));
    if show==1
        disp('Calculation points on L-curve')
    end
    for i=npoints-1:-1:1
        reg_param(i) = ratio*reg_param(i+1);
    end
    n=size(xi,2);
    for i=1:npoints
        f = Sn./(Sn + reg_param(i)^2);
        eta(i) = norm((f*ones(1,n)).*xi,'fro');
        rho(i) = norm(((1-f)*ones(1,n)).*beta,'fro');
    end
    
    % locate corner
    if show==1
        disp('Calculating curvature of L-curve')
    end
    % The L-curve is differentiable; computation of curvature in
    % log-log scale is easy.
    
    % Compute g = - curvature of L-curve.
    g = reglcfun(reg_param,Sn,beta,xi);
    
    % Locate the corner.  If the curvature is negative everywhere,
    % then define the leftmost point of the L-curve as the corner.
    if show==1
        disp('Searching for corner in L-curve')
        OPT=optimset('Display','iter');
    else
        OPT=optimset('Display','off');
    end
    [gmin,gi] = min(g);
    reg_c = fminbnd(@reglcfun,...
        reg_param(min(gi+1,length(g))),reg_param(max(gi-1,1)),...
        OPT,Sn,beta,xi); % Minimizer.
    kappa_max = - reglcfun(reg_c,Sn,beta,xi); % Maximum curvature.
    
    if (kappa_max < 0)
        lr = length(rho);
        reg_c = reg_param(lr);
        rho_c = rho(lr);
        eta_c = eta(lr);
    else
        f = Sn./(Sn + reg_c^2);
        eta_c = norm((f*ones(1,n)).*xi,'fro');
        rho_c = norm(((1-f)*ones(1,n)).*beta,'fro');
    end
    
    %%%%%%%%%%%%%%%%
    
    % Truncated SVD
elseif (strncmp(method,'tsvd',4) || strncmp(method,'TSVD',4))
    
    % spline toolbox needed for determination of the corner.
    SkipCorner = exist('splines','dir')~=7;
    txt = 'TSVD';
    marker='o';
    eta = zeros(N,1);
    rho = zeros(N,1);
    eta(1) = sum(xi(1,:).^2);
    for k=2:N
        eta(k) = eta(k-1) + sum(xi(k,:).^2);
    end
    eta = sqrt(eta);
    rho(N) = eps^2;
    for k=N-1:-1:1
        rho(k) = rho(k+1) + sum(beta(k+1,:).^2);
    end
    rho = sqrt(rho);
    reg_param = (1:N)';
    % Determine corner using Splines
    if (SkipCorner)
        reg_c = NaN;
    else
        % The L-curve is discrete and may include unwanted fine-grained
        % corners.  Use local smoothing, followed by fitting a 2-D spline
        % curve to the smoothed discrete L-curve.
        
        % Set default parameters for treatment of discrete L-curve.
        deg   = 2;  % Degree of local smooting polynomial.
        q     = 2;  % Half-width of local smoothing interval.
        order = 4;  % Order of fitting 2-D spline curve.
        
        % Neglect singular values less than s_thr.
        s_thr = eps;
        index = find(s > s_thr);
        rho_t = rho(index);
        eta_t = eta(index);
        reg_param_t = reg_param(index);
        
        % Convert to logarithms.
        lr = length(rho_t);
        lrho = log(rho_t);
        leta = log(eta_t);
        slrho = lrho;
        sleta = leta;
        
        % For all interior points k = q+1:length(rho)-q-1 on the discrete
        % L-curve, perform local smoothing with a polynomial of degree deg
        % to the points k-q:k+q.
        v = (-q:q)';
        A = zeros(2*q+1,deg+1);
        A(:,1) = ones(length(v),1);
        for j = 2:deg+1
            A(:,j) = A(:,j-1).*v;
        end
        for k = q+1:lr-q-1
            cr = A\lrho(k+v); slrho(k) = cr(1);
            ce = A\leta(k+v); sleta(k) = ce(1);
        end
        
        % Fit a 2-D spline curve to the smoothed discrete L-curve.
        sp = spmak(1:lr+order,[slrho';sleta']);
        pp = ppbrk(sp2pp(sp),[4,lr+1]);
        
        % Extract abscissa and ordinate splines and differentiate them.
        % Compute as many function values as default in spleval.
        P     = spleval(pp);  dpp   = fnder(pp);
        D     = spleval(dpp); ddpp  = fnder(pp,2);
        DD    = spleval(ddpp);
        ppx   = P(1,:);       ppy   = P(2,:);
        dppx  = D(1,:);       dppy  = D(2,:);
        ddppx = DD(1,:);      ddppy = DD(2,:);
        
        % Compute the corner of the discretized .spline curve via max. curvature.
        % No need to refine this corner, since the final regularization
        % parameter is discrete anyway.
        % Define curvature = 0 where both dppx and dppy are zero.
        k1    = dppx.*ddppy - ddppx.*dppy;
        k2    = (dppx.^2 + dppy.^2).^(1.5);
        I_nz  = find(k2 ~= 0);
        kappa = zeros(1,length(dppx));
        kappa(I_nz) = -k1(I_nz)./k2(I_nz);
        [kmax,ikmax] = max(kappa);
        x_corner = ppx(ikmax); y_corner = ppy(ikmax);
        
        % Locate the point on the discrete L-curve which is closest to the
        % corner of the spline curve.  Prefer a point below and to the
        % left of the corner.  If the curvature is negative everywhere,
        % then define the leftmost point of the L-curve as the corner.
        if (kmax < 0)
            reg_c = reg_param_t(lr);
            rho_c = rho_t(lr);
            eta_c = eta_t(lr);
        else
            index = find(lrho < x_corner & leta < y_corner);
            if ~isempty(index)
                [dummy,rpi] = min((lrho(index)-x_corner).^2 + (leta(index)-y_corner).^2);
                rpi = index(rpi);
            else
                [dummy,rpi] = min((lrho-x_corner).^2 + (leta-y_corner).^2);
            end
            reg_c = reg_param_t(rpi); rho_c = rho_t(rpi); eta_c = eta_t(rpi);
        end
    end
    
else
    error('Illegal method')
end

%%%%%%%%%%%%%%%%

% Plot
if show==1
    N=length(rho);
    loglog(rho(2:end-1),eta(2:end-1))
    ax = axis;
    ni = round(N/10);
    if (max(eta)/min(eta) > 10 || max(rho)/min(rho) > 10)
        loglog(rho,eta,marker,rho(ni:ni:N),eta(ni:ni:N),'x')
    else
        plot(rho,eta,marker,rho(ni:ni:N),eta(ni:ni:N),'x')
    end
    HoldState = ishold;
    hold on;
    for k = ni:ni:N
        text(rho(k),eta(k),num2str(reg_param(k)));
    end
    if ~(SkipCorner)
        loglog([min(rho)/100,rho_c],[eta_c,eta_c],':r',...
            [rho_c,rho_c],[min(eta)/100,eta_c],':r')
        title(['L-curve, ',txt,' corner at ',num2str(reg_c)]);
    else
        title('L-curve')
    end
    axis(ax)
    if (~HoldState)
        hold off
    end
    xlabel('residual norm || A x - b ||_2')
    ylabel('solution norm || x ||_2')
end
end

function g = reglcfun(lambda,Sn,beta,xi)
% reglcfun  Computes L-curve for reglcurve.
%           Auxiliary function for reglcurve.
%
% Written by Vincent Verdult, May 2004.
% Based on Regularization Tools by P. C. Hansen


% Initialization.
L=size(lambda,1);
n=size(xi,2);
phi  = zeros(L,1);
dphi = zeros(L,1);
psi  = zeros(L,1);
dpsi = zeros(L,1);
eta  = zeros(L,1);
rho  = zeros(L,1);

% Compute some intermediate quantities.
for i = 1:L
    f  = Sn./(Sn + lambda(i)^2);
    cf = 1 - f;
    eta(i) = norm((f*ones(1,n)).*xi,'fro');
    rho(i) = norm((cf*ones(1,n)).*beta,'fro');
    f1 = -2*f.*cf/lambda(i);
    f2 = -f1.*(3-4*f)/lambda(i);
    phi(i)  = sum(f.*f1.*sum(xi.^2,2));
    psi(i)  = sum(cf.*f1.*sum(beta.^2,2));
    dphi(i) = sum((f1.^2 + f.*f2).*sum(xi.^2,2));
    dpsi(i) = sum((-f1.^2 + cf.*f2).*sum(beta.^2,2));
end

% Now compute the first and second derivatives of eta and rho
% with respect to lambda;
deta  =  phi./eta;
drho  = -psi./rho;
ddeta =  dphi./eta - deta.*(deta./eta);
ddrho = -dpsi./rho - drho.*(drho./rho);

% Convert to derivatives of log(eta) and log(rho).
dlogeta  = deta./eta;
dlogrho  = drho./rho;
ddlogeta = ddeta./eta - (dlogeta).^2;
ddlogrho = ddrho./rho - (dlogrho).^2;

% Let g = curvature.
g = - (dlogrho.*ddlogeta - ddlogrho.*dlogeta)./...
    (dlogrho.^2 + dlogeta.^2).^(1.5);
end
