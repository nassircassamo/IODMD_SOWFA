function [S,info] = dordssarx(varargin)
%DORDSSARX delivers information about the order of the LTI state 
%   space model and acts as a pre-processor for DESTSSARX, which 
%   estimates the system matrices.
% 
% SYNOPSYS
%   [S,info] = dordssarx(u,y,s);
%   [S,info] = dordssarx(u,y,s,p);
%   [S,info] = dordssarx(u,y,s,p,options);
%   [S,info] = dordssarx([],y,s,n);
%
% DESCRIPTION
%   This function delivers information about the order of the LTI
%   state space model and acts as a pre-processor for DESTSSARX. 
%   The functions DORDSSARX and DESTSSARX provide an modified
%   implementation of the SSAR(X) subspace identification
%   algorithm. The parameters s and p represent the number of
%   block-rows of future and past data that are taken into
%   account. If the parameter p is not specified it is assumed that
%   the number of past and future block-rows are the same
%   (s=p). The additional options can be used to specify the
%   weighting of the matrix M that is used to estimate the row
%   space of the extended controllability matrix.
%
% OPTIONS
%   'weight'  - Different type of weighting the extended 
%               controllability matrix:
%
%       'lsq'    least squares     M = Rzp inv(Rzz)  
%       'rq'     rq weighting      M = Rzp inv(Rpp^1/2)
%       'cca'    corr analysis     M = inv(Rpp^1/2) Rzp inv(Rzz^1/2)
%
%

% (C) Karel Hinnen, February 2006



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse input arguments                                           % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weight = 'cca';

nin = nargin;
while ischar(varargin{nin-1}),
   switch lower(varargin{nin-1}),
      case 'weight', weight = varargin{nin};
   end
   nin = nin - 2;
end
error(nargchk(3,4,nin));

nout = nargout;
if not(ismember(nout,2))
  error('Wrong number of output arguments');
end

u = varargin{1};
y = varargin{2};
s = varargin{3};

switch nin,
 case 3,
  p = s; 
 case 4;
  p = varargin{4}; 
end

l = size(y,2);   % dimension of output signal
m = size(u,2);   % dimension of input signal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of Hankel matrices and RQ factorization            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~mod(s+p,2) && (exist('fastr','file') == 3)
    R = fastr(u,y,(s+p)/2,'F')';
else
    if (exist('fastr','file') == 3)
        warning('SSARX:evenRQ','In order to perform fast RQ factorization s+p should be even');
    end

    Upf = bhankel(u',(s+p));
    Ypf = bhankel(y',(s+p));

    H = [Upf' Ypf'];
    R = triu(qr(H,0))';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate a high order AR model                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Uc = R((p+s-1)*m+1:(p+s)*m,1:(p+s)*m+(p+s-1)*l);
Up = R(1:(p+s-1)*m,1:(p+s)*m+(p+s-1)*l);

Yc = R((p+s)*m+(p+s-1)*l+1:(p+s)*(m+l),1:(p+s)*m+(p+s-1)*l);
Yp = R((p+s)*m+1:(p+s)*m+(p+s-1)*l,1:(p+s)*m+(p+s-1)*l);


ARX = Yc/[Yp; Up; Uc];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct PHI1=(I-tilde PSI) matrix                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PHI1 = eye(s*l,s*l);

for ii=(1:min(s,p)-1),
   PHI1(ii*l+1:(ii+1)*l,1:ii*l) = -ARX(:,(s+p-ii-1)*l+1:(s+p-1)*l);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct PHI2 = tilde PHI matrix                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PHI2 = zeros(s*l,s*m);

for jj = 1:min(s,p), 
   PHI2((jj-1)*l+1:jj*l,1:jj*m) = ...
       ARX(:,(s+p-1)*l+(s+p-jj)*m+1:(s+p-1)*(l+m)+m);   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of the column space of Gamma                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch lower(weight);
   case 'lsq'      % M = R_zp (R_pp)^(-1)
      fprintf('Specified weighting: lsq \n');

      Uf = R(p*m+1:(s+p)*m,1:(s+p)*m+p*l);
      Up = R(1:p*m,1:(s+p)*m+p*l);
      Yf = R((s+p)*m+p*l+1:(s+p)*(m+l),1:(s+p)*m+p*l);
      Yp = R((s+p)*m+1:(s+p)*m+p*l,1:(s+p)*m+p*l);
      
      Zf = PHI1*Yf - PHI2*Uf;

      M  = Zf/([Up;Yp]);      
      [U,S,V] = svd(M);
      
   case 'rq'    % M = R_zp (R_pp)^(-1/2)
      fprintf('Specified weighting: rq \n');

      Uf = R(p*m+1:(s+p)*m,1:(s+p)*m+p*l);
      Up = R(1:p*m,1:(s+p)*m+p*l);
      Yf = R((s+p)*m+p*l+1:(s+p)*(m+l),1:(s+p)*m+p*l);
      Yp = R((s+p)*m+1:(s+p)*m+p*l,1:(s+p)*m+p*l);
      
      Zf = PHI1*Yf - PHI2*Uf;
      
      M = Zf;            
      [U,S,V] = svd(M);

   case 'cca'   % M = (R_zz)^(-1/2) R_zp (R_pp)^(-1/2)
      fprintf('Specified weighting: cca \n');

      Uf = R(p*m+1:(s+p)*m,1:(s+p)*(m+l));
      Up = R(1:p*m,1:(s+p)*m+p*l);
      Yf = R((s+p)*m+p*l+1:(s+p)*(m+l),1:(s+p)*(m+l));
      Yp = R((s+p)*m+1:(s+p)*m+p*l,1:(s+p)*m+p*l);
      
      Zf  = PHI1*Yf - PHI2*Uf;
      Rzf = triu(qr(Zf',0))';
      
      M  = Rzf(:,1:s*l)\Zf(:,1:(s+p)*m+p*l);

      [U,S,V] = svd(M);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return stucture which is used as input for est_ssarx            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info.weight = weight;
info.p = p; info.s = s;
info.m = m; info.l = l;  
info.U = U; info.S = S; info.V = V;
info.Up = Up; info.Yp = Yp;
info.R = R;

S = diag(S);
