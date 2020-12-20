function sys = c2d(sys,Ts,method,approx,varargin)
%IDAFFLPV/C2D  Continuous-to-discrete conversion of state-space models.
%  MD = C2D(MC,Ts,METHOD,APPROX) converts the continuous-time IDAFFLPV
%  model MC to a discrete-time model MD with sample time Ts. The string
%  METHOD selects the discretization method among the following:
%    'zoh'       Zero-order hold on the inputs
%    'euler'     Euler approximation
%    'tustin'    Bilinear (Tustin) approximation
%    'prewarp'   Tustin approximation with frequency prewarping
%  The critical frequency Wc (in rad/sec) is specified as fifth input by
%  MD = C2D(MC,Ts,'prewarp',APPROX,Wc). The string APPROX selects the
%  approximation for the Taylor approximation method to retrieve an affine
%  linear parameter varying system among the following:
%    'linear'    linear approximation
%    'quadratic' quadratic approximation
%
%  References:
%    [1] Tóth, R., P. S. C. Heuberger and P. M. J. Van den Hof, "On the 
%    Discretization of LPV state-space representations, 
%    IET Control Theory & Applications, 2010.

% Check number of inputs
ni = nargin;
error(nargchk(2,4,ni))

% Error checking
if ~isct(sys)
   error('System is already discrete')
elseif ~(length(Ts)==1 || Ts > 0)
   error('Sample time TS must be a positive scalar.')
elseif ni==2 || isempty(method),
   method = 'euler';
elseif ~ischar(method)
   error('METHOD must be a string.')
else
   m1 = lower(method);
   if strcmp(m1,'prewarp')
      if ni<4
         error('The critical frequency must be specified when using prewarp method.')
      elseif ~(isPositiveScalar(varargin{1}) && varargin{1}<pi/Ts)
         error('Prewarping frequency must be a positive scalar smaller than the Nyquist frequency.')
      end
   end
   method = m1;
end
if nargin < 4 || isempty(approx)
    approx = 'linear';
end

% Get system matrices and sizes
[a,b,c,d,k] = getABCDK(sys);
[ny,nu,nx,np] = size(sys);

% Discretize using method(1)
switch lower(method)
    case 'zoh'
        % Create a symbolic system matrices
        as = a(:,1:nx);
        bs = b(:,1:nu);
        ks = k(:,1:ny);
        sched = cell(1,np);
        for p = 1:np
            sched{p} = sym(['p' int2str(p)],'real');
            as = as + sched{p}.*a(:,p*nx+1:(p+1)*nx);
            bs = bs + sched{p}.*b(:,p*nu+1:(p+1)*nu);
            ks = ks + sched{p}.*k(:,p*ny+1:(p+1)*ny);
        end
        
        % ZOH discretization method
        ae = expm(as.*Ts);
        be = as\(ae - eye(nx))*bs;
        ke = as\(ae - eye(nx))*ks;

        % Multivariable Taylor approximation
        var = '[';
        for p = 1:np
            if p == np
                var = [var char(sched{p}) ']'];
            else
                var = [var char(sched{p}) ','];
            end
        end
        z = eye(np);
        if strcmpi(approx,'linear')
            for i = 1:nx
                for j = 1:nx
                    ar = feval(symengine,'mtaylor',ae(i,j),var,2);
                    %ar = maple('mtaylor',ae(i,j),var,2);
                    a(i,j) = subs(ar,sched,zeros(1,np));
                    for p = 1:np
                        a(i,p*nx+j) = subs(ar,sched,z(p,:)) - a(i,j);
                    end
                end
                for j = 1:nu
                    br = feval(symengine,'mtaylor',be(i,j),var,2);
                    %br = maple('mtaylor',be(i,j),var,2);
                    b(i,j) = subs(br,sched,zeros(1,np));
                    for p = 1:np
                        b(i,p*nu+j) = subs(br,sched,z(p,:)) - b(i,j);
                    end
                end
                for j = 1:ny
                    kr = feval(symengine,'mtaylor',ke(i,j),var,2);
                    %br = maple('mtaylor',be(i,j),var,2);
                    k(i,j) = subs(kr,sched,zeros(1,np));
                    for p = 1:np
                        k(i,p*nu+j) = subs(kr,sched,z(p,:)) - k(i,j);
                    end
                end
            end
        elseif strcmpi(approx,'quadratic')
            a = zeros(nx,(1+np+sum(1:np))*nx);
            b = zeros(nx,(1+np+sum(1:np))*nu);
            k = zeros(nx,(1+np+sum(1:np))*ny);
            c = [c zeros(ny,sum(1:np)*nx)];
            d = [d zeros(ny,sum(1:np)*(nu+ny))];       
            for i = 1:nx
                for j = 1:nx
                    a(i,j) = double(subs(ae(i,j),sched,zeros(1,np)));
                    h = 1;
                    for p = 1:np
                        Pdiff = diff(ae(i,j),sched{p});
                        a(i,p*nx+j) = double(subs(Pdiff,sched,zeros(1,np)));
                        for g = p:np
                            if g == p
                                a(i,(1+np)*nx+(h-1)*nx+j) = double(subs(diff(Pdiff,sched{g}),sched,zeros(1,np)))./2;
                            else
                                a(i,(1+np)*nx+(h-1)*nx+j) = double(subs(diff(Pdiff,sched{g}),sched,zeros(1,np)));
                            end
                            h = h + 1;
                        end
                    end
                end
                for j = 1:nu
                    b(i,j) = double(subs(be(i,j),sched,zeros(1,np)));
                    h = 1;
                    for p = 1:np
                        Pdiff = diff(be(i,j),sched{p});
                        b(i,p*nu+j) = double(subs(Pdiff,sched,zeros(1,np)));
                        for g = p:np
                            if g == p
                                b(i,(1+np)*nu+(h-1)*nu+j) = double(subs(diff(Pdiff,sched{g}),sched,zeros(1,np)))./2;
                            else
                                b(i,(1+np)*nu+(h-1)*nu+j) = double(subs(diff(Pdiff,sched{g}),sched,zeros(1,np)));
                            end
                            h = h + 1;
                        end
                    end
                end
                for j = 1:ny
                    k(i,j) = double(subs(ke(i,j),sched,zeros(1,np)));
                    h = 1;
                    for p = 1:np
                        Pdiff = diff(ke(i,j),sched{p});
                        k(i,p*ny+j) = double(subs(Pdiff,sched,zeros(1,np)));
                        for g = p:np
                            if g == p
                                k(i,(1+np)*ny+(h-1)*ny+j) = double(subs(diff(Pdiff,sched{g}),sched,zeros(1,np)))./2;
                            else
                                k(i,(1+np)*ny+(h-1)*ny+j) = double(subs(diff(Pdiff,sched{g}),sched,zeros(1,np)));
                            end
                            h = h + 1;
                        end
                    end
                end
            end
        else
            error('APPROX not recognized!')
        end
              
    case {'euler'}
        bt = zeros(nx,(np+1)*(nu+ny));
        for i = 1:np+1
            bt(:,(i-1)*(nu+ny)+1:i*(nu+ny)) = [b(:,(i-1)*nu+1:i*nu) k(:,(i-1)*ny+1:i*ny)];
        end
        b = bt;
        
        if strcmpi(approx,'linear')
            % Rectangular approximation method (Euler)
            a = [eye(nx) zeros(nx,np*nx)] + a.*Ts;
            b = b.*Ts;
            
            bt = zeros(nx,(np+1)*nu);
            k = zeros(nx,(np+1)*ny);
            for i = 1:np+1
                bt(:,(i-1)*nu+1:i*nu) = b(:,(i-1)*(nu+ny)+1:(i-1)*(nu+ny)+nu);
                k(:,(i-1)*ny+1:i*ny) = b(:,(i-1)*(nu+ny)+nu+1:i*(nu+ny));
            end
            b = bt;

        elseif strcmpi(approx,'quadratic')
            % Quadratic approximation method (Hanselmann)
            a = [eye(nx) zeros(nx,(np+sum(1:np))*nx)] + [a zeros(nx,sum(1:np)*nx)].*Ts + mconv(a,a,nx,nx,np).*((Ts^2)/2);
            b = mconv(a.*(Ts/2),b,nx,(nu+ny),np).*Ts + [b zeros(nx,sum(1:np)*(nu+ny))].*Ts;
            c = [c zeros(ny,sum(1:np)*nx)];
            d = [d zeros(ny,sum(1:np)*nu)];
            
            bt = zeros(nx,(np+sum(1:np)+1)*nu);
            k = zeros(nx,(np+sum(1:np)+1)*ny);
            for i = 1:np+sum(1:np)+1
                bt(:,(i-1)*nu+1:i*nu) = b(:,(i-1)*(nu+ny)+1:(i-1)*(nu+ny)+nu);
                k(:,(i-1)*ny+1:i*ny) = b(:,(i-1)*(nu+ny)+nu+1:i*(nu+ny));
            end
            b = bt;
        else
            error('APPROX not recognized!')
        end
    case {'tustin','prewarp'}
        % Create a symbolic system matrices
        as = a(:,1:nx);
        bs = b(:,1:nu);
        ks = b(:,1:ny);
        cs = c(:,1:nx);
        ds = d(:,1:nu);
        sched = cell(1,np);
        for p = 1:np
            sched{p} = sym(['p' int2str(p)],'real');
            as = as + sched{p}.*a(:,p*nx+1:(p+1)*nx);
            bs = bs + sched{p}.*b(:,p*nu+1:(p+1)*nu);
            ks = ks + sched{p}.*k(:,p*ny+1:(p+1)*ny);
            cs = cs + sched{p}.*c(:,p*nx+1:(p+1)*nx);
            ds = ds + sched{p}.*d(:,p*nu+1:(p+1)*nu);
        end
        
        % Tustin discretization method
        if strncmp(method,'t',1)
            T = Ts;
        else
            % Handle prewarping
            w = varargin{1};
            T = 2*tan(w*Ts/2)/w;
        end
        at = (eye(nx) - (T/2).*as)\(eye(nx) + (T/2).*as);
        bt = ((eye(nx) - (T/2).*as)\bs);
        kt = ((eye(nx) - (T/2).*as)\ks);
        ct = T.*(cs/(eye(nx) - (T/2).*as));
        dt = (T/2).*(cs*bt) + ds;

        % Multivariable Taylor approximation
        var = '[';
        for p = 1:np
            if p == np
                var = [var char(sched{p}) ']'];
            else
                var = [var char(sched{p}) ','];
            end
        end
        z = eye(np);
        if strcmpi(approx,'linear')
            for i = 1:nx
                for j = 1:nx
                    ar = feval(symengine,'mtaylor',at(i,j),var,2);
                    %ar = maple('mtaylor',ae(i,j),var,2);
                    a(i,j) = subs(ar,sched,zeros(1,np));
                    for p = 1:np
                        a(i,p*nx+j) = subs(ar,sched,z(p,:)) - a(i,j);
                    end
                end
                for j = 1:nu
                    br = feval(symengine,'mtaylor',bt(i,j),var,2);
                    %br = maple('mtaylor',be(i,j),var,2);
                    b(i,j) = subs(br,sched,zeros(1,np));
                    for p = 1:np
                        b(i,p*nu+j) = subs(br,sched,z(p,:)) - b(i,j);
                    end
                end
                for j = 1:ny
                    kr = feval(symengine,'mtaylor',kt(i,j),var,2);
                    %br = maple('mtaylor',be(i,j),var,2);
                    k(i,j) = subs(kr,sched,zeros(1,np));
                    for p = 1:np
                        k(i,p*nu+j) = subs(kr,sched,z(p,:)) - k(i,j);
                    end
                end
            end
            for i = 1:ny
                for j = 1:nx
                    cr = feval(symengine,'mtaylor',ct(i,j),var,2);
                    %cr = maple('mtaylor',ct(i,j),var,2);
                    c(i,j) = subs(cr,sched,zeros(1,np));
                    for p = 1:np
                        c(i,p*nx+j) = subs(cr,sched,z(p,:)) - c(i,j);
                    end
                end
                for j = 1:nu
                    dr = feval(symengine,'mtaylor',dt(i,j),var,2);
                    %dr = maple('mtaylor',dt(i,j),var,2);
                    d(i,j) = subs(dr,sched,zeros(1,np));
                    for p = 1:np
                        d(i,p*nu+j) = subs(dr,sched,z(p,:)) - d(i,j);
                    end
                end
            end
        elseif strcmpi(approx,'quadratic')
            a = zeros(nx,(1+np+sum(1:np))*nx);
            b = zeros(nx,(1+np+sum(1:np))*nu);
            k = zeros(nx,(1+np+sum(1:np))*ny);
            c = [c zeros(ny,sum(1:np)*nx)];
            d = [d zeros(ny,sum(1:np)*nu)];       
            for i = 1:nx
                for j = 1:nx
                    a(i,j) = double(subs(at(i,j),sched,zeros(1,np)));
                    h = 1;
                    for p = 1:np
                        Pdiff = diff(at(i,j),sched{p});
                        a(i,p*nx+j) = double(subs(Pdiff,sched,zeros(1,np)));
                        for g = p:np
                            if g == p
                                a(i,(1+np)*nx+(h-1)*nx+j) = double(subs(diff(Pdiff,sched{g}),sched,zeros(1,np)))./2;
                            else
                                a(i,(1+np)*nx+(h-1)*nx+j) = double(subs(diff(Pdiff,sched{g}),sched,zeros(1,np)));
                            end
                            h = h + 1;
                        end
                    end
                end
                for j = 1:nu
                    b(i,j) = double(subs(bt(i,j),sched,zeros(1,np)));
                    h = 1;
                    for p = 1:np
                        Pdiff = diff(bt(i,j),sched{p});
                        b(i,p*nu+j) = double(subs(Pdiff,sched,zeros(1,np)));
                        for g = p:np
                            if g == p
                                b(i,(1+np)*nu+(h-1)*nu+j) = double(subs(diff(Pdiff,sched{g}),sched,zeros(1,np)))./2;
                            else
                                b(i,(1+np)*nu+(h-1)*nu+j) = double(subs(diff(Pdiff,sched{g}),sched,zeros(1,np)));
                            end
                            h = h + 1;
                        end
                    end
                end
                for j = 1:ny
                    k(i,j) = double(subs(kt(i,j),sched,zeros(1,np)));
                    h = 1;
                    for p = 1:np
                        Pdiff = diff(kt(i,j),sched{p});
                        k(i,p*ny+j) = double(subs(Pdiff,sched,zeros(1,np)));
                        for g = p:np
                            if g == p
                                k(i,(1+np)*ny+(h-1)*ny+j) = double(subs(diff(Pdiff,sched{g}),sched,zeros(1,np)))./2;
                            else
                                k(i,(1+np)*ny+(h-1)*ny+j) = double(subs(diff(Pdiff,sched{g}),sched,zeros(1,np)));
                            end
                            h = h + 1;
                        end
                    end
                end
            end
            for i = 1:ny
                for j = 1:nx
                    c(i,j) = double(subs(ct(i,j),sched,zeros(1,np)));
                    h = 1;
                    for p = 1:np
                        Pdiff = diff(ct(i,j),sched{p});
                        c(i,p*nx+j) = double(subs(Pdiff,sched,zeros(1,np)));
                        for g = p:np
                            if g == p
                                c(i,(1+np)*nx+(h-1)*nx+j) = double(subs(diff(Pdiff,sched{g}),sched,zeros(1,np)))./2;
                            else
                                c(i,(1+np)*nx+(h-1)*nx+j) = double(subs(diff(Pdiff,sched{g}),sched,zeros(1,np)));
                            end
                            h = h + 1;
                        end
                    end
                end
                for j = 1:nu
                    d(i,j) = double(subs(dt(i,j),sched,zeros(1,np)));
                    h = 1;
                    for p = 1:np
                        Pdiff = diff(dt(i,j),sched{p});
                        d(i,p*nu+j) = double(subs(Pdiff,sched,zeros(1,np)));
                        for g = p:np
                            if g == p
                                d(i,(1+np)*nu+(h-1)*nu+j) = double(subs(diff(Pdiff,sched{g}),sched,zeros(1,np)))./2;
                            else
                                d(i,(1+np)*nu+(h-1)*nu+j) = double(subs(diff(Pdiff,sched{g}),sched,zeros(1,np)));
                            end
                            h = h + 1;
                        end
                    end
                end
            end 
        else
            error('APPROX not recognized!');
        end
    otherwise
        error('METHOD not recognized!')
end

% change scheduling names
if strcmpi(approx,'quadratic')
    Names = sys.SchedulingName;
    NewNames = cell(np+sum(1:np),1);
    h = 1;
    for p = 1:np
        if isempty(Names{p,1})
            P = sprintf('p%d',p);
        else
            P = Names{p,1};
        end
        NewNames{p,1} = P;
        for j = p:np
            if isempty(Names{j,1})
                K = sprintf('p%d',j);
            else
                K = Names{j,1};
            end
            if p == j
                NewNames{np+h,1} = [P '^2'];
            else
                NewNames{np+h,1} = [P '*' K];
            end
            h = h + 1;
        end
    end
    sys.SchedulingName = NewNames;
end

% Store results
sys.a = a;
sys.b = b;
sys.c = c;
sys.d = d;
sys.k = k;
sys.Ts = Ts;

end % end of function C2D 

% matrix convolution
function c = mconv(a,b,nx,nu,np)
    c = zeros(nx,(1+np+sum(1:np))*nu);
    h = 1;
    for p = 1:np+1
        for k = p:np+1
            if p == k
                c(:,(h-1)*nu+1:h*nu) = a(:,(p-1)*nx+1:p*nx)*b(:,(k-1)*nu+1:k*nu);
            else
                c(:,(h-1)*nu+1:h*nu) = a(:,(p-1)*nx+1:p*nx)*b(:,(k-1)*nu+1:k*nu) + a(:,(k-1)*nx+1:k*nx)*b(:,(p-1)*nu+1:p*nu);     
            end
            h = h + 1;
        end
    end    
end








