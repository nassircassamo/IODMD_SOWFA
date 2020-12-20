function [Sn,R,Znew] = dordeiv(u,y,r,s,Zold)
%DORDEIV    Delivers information about the order of the LTI
%           state space model and acts as a pre-processor for
%           dmodeiv and dac2bd_eiv/dac2bd_b. The latter functions estimate
%           the quadruple of system matrices.
% 
% Model structure for the errors-in-variables problem: 
%              x(k+1) = Ax(k) + Bu~(k) + f(k)
%              y~(k)  = Cx(k) + Du~(k) 
%           with measurements 
%              u(k) = u~(k) + w(k)
%              y(k) = y~(k) + v(k)
%           where f(k), w(k) and v(k) are zero-mean white noise sequences
%           independent of the input u~(j) for k >= j. 
%
%           Notes: 
%           1. The plant can be operated under either open-loop or 
%              closed-loop. For closed-loop operation, r(k) is an 
%              external reference input. 
%           2. For open-loop operation, do NOT use this routine if 
%              u~(k) is white noise, use dordpo and dmodpo instead.
% 
% Syntax:   For open-loop data: 
%           [Sn,R,Znew]=dordeiv(u,y,[],s);
%           [Sn,R,Znew]=dordeiv(u,y,[],s,Zold);
% 
%           For closed-loop data: 
%           If r is not available, use the same syntax as above, otherwise
%           [Sn,R,Znew]=dordeiv(u,y,r,s);
%           [Sn,R,Znew]=dordeiv(u,y,r,s,Zold);
%           
% Input:
%   u,y     Input and output data of the system to be identified
%   r       Closed-loop reference input
%   s       The dimension parameter that determines the number
%           of block rows in the Hankel matrices
%   Zold    Information on system given by previous calls of 
%           dordeiv. See smidemo2 or the manual for more details.
%
% Output:
%   Sn      Singular values bearing information on the order
%           of the system.
%   R       dmodeiv uses R to estimate A and C matrices 
%   Znew    Information on system which can be used by 
%           future calls of dordeiv. See smidemo2 or the 
%           manual for more details. 
%
% See also: dmodeiv, dac2bd_eiv, dac2b_eiv, dac2bd_cl, dac2b_cl
% 
% Reference: C.T. Chou and M. Verhaegen 
%            Subspace algorithms fpr the identification of
%            multivariable dynamic errors-in-variables models
%            Automatica vol 33, no 10, pp. 1857-1869, October, 1997

% C.T. Chou, Oct 1997
% Revised by Ivo Houtzager, 2009
% Copyright (c) 1997-2009, Delft Center of Systems and Control 

[N,L] =size(y);
m=size(u,2);
cr=size(r,2);
nZ = 2*s*(m+L+cr);

if ~(size(u,1)==N) 
  error('Input and output should have same length')
end

if (~isempty(r) && ~(size(r,1)==N))
  error('Input and reference should have same length')
end   

if s>=N/2
  error('s is chosen too large or number of datapoints is too small')
end

N=N-2*s+1;

if nargin < 5, 
  Zold = []; 
elseif ~((size(Zold,1)==nZ) && (size(Zold,2)==nZ))
  error('Z-matrix has unexpected size.')
end

%
% construction of Hankel matrices
%
Y=[];U=[];E=[];
if isempty(r)
  for i=(1:2*s) 
    U(:,(i-1)*m+1:i*m)=u(i:N+i-1,:);
    Y(:,(i-1)*L+1:i*L)=y(i:N+i-1,:);
  end    
else   
  for i=(1:2*s) 
    U(:,(i-1)*m+1:i*m)=u(i:N+i-1,:);
    Y(:,(i-1)*L+1:i*L)=y(i:N+i-1,:);
    E(:,(i-1)*m+1:i*m)=r(i:N+i-1,:);
  end    
end       

upr = 1:m*s;
ufr = m*s+1:2*m*s;
ypr = 1:L*s;
yfr = L*s+1:2*L*s;

Z = triu(qr([U(:,ufr) U(:,upr) Y(:,ypr) E Y(:,yfr); Zold]));
Znew = Z(1:nZ,1:nZ);
Zx = Znew';

F0 = Zx([1:s*m s*2*m+s*L+2*s*cr+1:nZ],1:s*2*m+s*L+2*s*cr);
B0 = Zx(s*m+1:s*2*m+s*L+2*s*cr,1:s*2*m+s*L+2*s*cr)';

[Q1,B1,E1] = qr(B0);
B1 = B1(1:size(B1,2),:);
rB1 = max(rank(B1),(m+L)*s);
B1 = B1(1:rB1,:);

F1 = F0*Q1;
F1 = F1(:,1:rB1);

F2 = triu(qr(F1'))';
rr = (s*m+1):(s*(m+L));
R = F2(m*s+1:(m+L)*s,m*s+1:(m+L)*s);
[Un,Sn] = svd(R);
Sn = diag(Sn);
Sn = Sn(1:s);
R(1,2) = m;
R(1,3) = L;
R(2,3) = s;

% END OF THE CALCULATIONS






