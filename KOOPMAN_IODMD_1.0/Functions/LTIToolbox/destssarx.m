function [A,B,C,D,K,Re] = destssarx(varargin)
%DESTSSARX estimates the system matrices of a LTI state space model using 
%  the output from DORDSSARX.
%               
%          
% SYNOPSYS
%   [a,b,c,d,k,re] = destssarx(info,n);
%   [a,c,k,re] = destssarx(info,n);
%
% DESCRIPTION
%  Estimates the system matrices of a LTI state space model using
%  the output from DORDSSARX. The functions DORDSSARX and DESTSSARX
%  provide an modified implementation of the SSAR(X) subspace
%  identification algorithm.
%
%

% (C) Karel Hinnen, February 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse input arguments                                           % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error(nargchk(2,2,nargin));

nout = nargout;
if not(ismember(nout,[3,4,5,6]))
  error('Wrong number of output arguments');
end

n = varargin{2};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unpack data structure produced by ord_ssarx                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info = varargin{1};

weight = info.weight;
p = info.p; s = info.s;
m = info.m; l = info.l;  
U = info.U; S = info.S; V = info.V;
Up = info.Up;  Yp = info.Yp;
R = info.R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine Kappa                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch lower(weight);
 case 'lsq'     
  Kappa   = sqrt(S(1:n,1:n))*V(:,1:n)';
 case 'rq'    
  Kappa   = sqrt(S(1:n,1:n))*V(:,1:n)'/[Up;Yp];
 case 'cca' 
  Kappa   = sqrt(S(1:n,1:n))*V(:,1:n)'/[Up;Yp];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate state sequence and shifted version                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Up  = R(1:p*m,1:(s+p)*m+p*l+l);
Upp = R(m+1:p*m+m,1:(s+p)*m+p*l+l);
Uc  = R(p*m+1:p*m+m,1:(s+p)*m+p*l+l);

Yp  = R((s+p)*m+1:(s+p)*m+p*l,1:(s+p)*m+p*l+l);
Ypp = R((s+p)*m+l+1:(s+p)*m+p*l+l,1:(s+p)*m+p*l+l);
Yc  = R((s+p)*m+p*l+1:(s+p)*m+p*l+l,1:(s+p)*m+p*l+l);

Xp = Kappa*[Up;Yp];
Xf = Kappa*[Upp;Ypp];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine state space matrices by linear regression             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ABK = Xf/[Xp;Uc;Yc];

Ae = ABK(:,1:n);
Be = ABK(:,n+1:n+m);
K  = ABK(:,n+m+1:n+m+l);

CD = Yc/[Xp;Uc];
C  = CD(:,1:n);
D  = CD(:,n+1:n+m);

A = Ae + K*C;
B = Be + K*D;

if (nout == 4 && m == 0) || (nout == 6 && m > 0),
   E  = Yc(:,1:end) - CD*[Xp;Uc];
   Re = cov(E');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return only A, K and C matrices in output only case             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nout == 3 || nout == 4) && (m == 0),
   B = C;   
   C = K;

   if (nout == 4 && m == 0),
      D = Re;
   end
   
end
