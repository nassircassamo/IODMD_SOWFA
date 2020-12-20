function [B,D,Rnew] = dac2bd_eiv(A,C,u,y,r,NN,Rold)
%DAC2BD_EIV A routine to solve the dynamic errors-in-variables problem. 
%           Estimates the B and D matrices of LTI state 
%           space model in innovations form 
%              x(k+1) = Ax(k) + Bu~(k) + f(k)
%              y~(k)  = Cx(k) + Du~(k) 
%           given noisy input/output measurements 
%              u(k) = u~(k) + w(k)
%              y(k) = y~(k) + v(k)
%           where f(k), w(k) and v(k) are zero-mean white noise sequences
%           independent of the input u~(j) for k >= j. The plant can
%           be operated under either open-loop or closed-loop. For
%           closed-loop operation, r(k) is an external reference input.
% 
% Syntax:   
% 	    [B,D,Rnew] = dac2bd_eiv(A,C,u,y)
% 	    [B,D,Rnew] = dac2bd_eiv(A,C,u,y,r,NN,Rold)
% 					
% Input:
% A,C       Estimated A and C matrices from dmodeiv
% u,y       Noisy input and output data 
% r         External reference signal (only for closed-loop data) 
%           Optional argument. Default = []. 
% NN        Number of extended instrumental variables. 
%           If u and y are specfied, it should be a 2-element vector.
%           If r is also present, then it should be a 3-element vector. 
%           Roughly speaking, NN should be chosen such that 
%           NN(1) x #outputs + NN(2) x #inputs > 
%           #states + #states x #inputs + #inputs x #outputs 
%           If not specified or entered as [], the program 
%           determines a suitable NN based on given data
% Rold      Compressed data matrix from the previous call of eiv_bd. 
% 
% Output:
% B,D       Estimated system matrices.
% Rnew      Compressed data matrix for the next call of eiv_bd. 
%           Used when analyzing multiple data sequences.
% 
% See also: dac2bd_eiv, dac2b_cl, dac2bd_cl, dordeiv, dmodeiv 
% 
% Reference: C.T. Chou and M. Verhaegen 
%            Subspace algorithms fpr the identification of
%            multivariable dynamic errors-in-variables models
%            Automatica vol 33, no 10, pp. 1857-1869, October, 1997
% 

% C.T. Chou, Oct 1997
% Revised by Ivo Houtzager, 2009
% Copyright (c) 1997-2009, Delft Center of Systems and Control 

% error checks
if nargin < 4, error('Too few inputs'), end

n=size(A,1);
N=size(y,1);
L=size(y,2);
m=size(u,2);

if ~(size(A,2)==n)
    error('A should be square')
end

if ~(size(C,2)==n)
    error('C matrix should have same number of columns as A')
end

if ~(size(u,1)==N)
    error('Input and output should have same length')
end

if nargin < 5
    r = [];
    nr = 0;
else
    nr = size(r,2);
end

if (~isempty(r) && ~(size(u,1)==N))
    error('Input and reference should have same length')
end

if nargin < 6
    NN = [];
end

Nparaa = m*(n+L);
Npara  = Nparaa + n;

if isempty(NN)
    N1 = ceil(2*Npara/(L^2+0.5*L*m));
    N2 = N1;
    if isempty(r)
        N3 = 0;
    else
        N3 = N2;
    end
else
    NNN = length(NN);
    if (NNN < 2) || (NNN > 3)
        error('NN should have at least 2 elements')
    elseif (NNN == 2)
        N1 = NN(1);
        N2 = NN(2);
        if isempty(r),
            N3 = 0;
        else
            disp('The number of blocks of R is not specified. Default is used.')
            N3 = N2;
        end
    elseif (NNN == 3)
        N1 = NN(1);
        N2 = NN(2);
        N3 = NN(3);
    end
end
Nf = floor((N-N1+1)/2);
tau = N - Nf;

if nargin < 7,
    Rold = [];
elseif ~((size(Rold,1)==Nparaa) && (size(Rold,2)==Nparaa+1))
    error('Rold matrix has unexpected size.')
end

% check if A unstable, if so change the algorithm
mp = max(abs(eig(A)));
if (mp >= 1)
    disp('Note from eiv_bd: A matrix is unstable so closed-loop data is assumed.')
    K = place(A',C',(0:n-1)/n);
    K = K';
    A = A - K * C;
    z = dlsim(A,K,C,zeros(L),y);
    y = y - z;
end

% form the data equation and the instrumental variable matrices
% LongY = [Gamma  S1  S2] x [x1; vec(B); vec(D)]
% IVy, IVu
LongY = zeros(Nf*L,1);
Gamma = zeros(Nf*L,n);
S1 = zeros(Nf*L,n*m);
S2 = zeros(Nf*L,L*m);
IVy = zeros(Nf*L,N1*L^2);
IVu = zeros(Nf*L,N2*L*m);
Ivr = zeros(Nf*L,N3*L*nr);

for ri = 1:Nf
    rr = ((ri-1)*L+1):(ri*L);   % row range
    tpi = ri;                   % time index --- past
    tfi = tau + tpi;            % time index --- future
    
    LongY(rr,1) = y(tfi,:)';
    
    if (ri > 1)
        Block1 = Block1 * A;
    else
        Block1 = C;
    end
    Gamma(rr,:) = Block1;
    
    if (ri > 1)
        for ii = 1:m
            iic1 = (ii-1)*n+1;
            iic2 = ii*n;
            Block2(:,iic1:iic2) = Block2(:,iic1:iic2) * A + C * u(tfi-1,ii);
        end
        S1(rr,:) = Block2;
    else
        Block2 = zeros(L,m*n);
    end
    
    ii = 1; iicr = 1:L:m*L; DBlock = zeros(L,m*L);
    DBlock(ii,iicr) = u(tfi,:);
    for ii = 2:L
        iicr = iicr + 1;
        DBlock(ii,iicr) = u(tfi,:);
    end
    S2(rr,:) = DBlock;
    
    for j = 1:N1,
        ii = 1; iicr = 1:L:L*L; DBlock = zeros(L,L*L);
        DBlock(ii,iicr) = y(tpi+j-1,:);
        for ii = 2:L
            iicr = iicr + 1;
            DBlock(ii,iicr) = y(tpi+j-1,:);
        end
        cr = ((j-1)*L^2+1):(j*L^2);
        IVy(rr,cr) = DBlock;
    end
    
    for j = 1:N2,
        ii = 1; iicr = 1:L:m*L; DBlock = zeros(L,m*L);
        DBlock(ii,iicr) = u(tpi+j-1,:);
        for ii = 2:L
            iicr = iicr + 1;
            DBlock(ii,iicr) = u(tpi+j-1,:);
        end
        cr = ((j-1)*L*m+1):(j*L*m);
        IVu(rr,cr) = DBlock;
    end
    
    for j = 1:N3,
        ii = 1; iicr = 1:L:nr*L; DBlock = zeros(L,nr*L);
        DBlock(ii,iicr) = r(tpi+j-1,:);
        for ii = 2:L
            iicr = iicr + 1;
            DBlock(ii,iicr) = r(tpi+j-1,:);
        end
        cr = ((j-1)*L*nr+1):(j*L*nr);
        IVr(rr,cr) = DBlock;
    end
    
end  % end `ri' loop

Rx = triu(qr([IVy IVu Ivr Gamma S1 S2 LongY]));
part1 = N1*L^2+N2*L*m+N3*L*nr;
part2 = n;
part3 = n*m+L*m+1;
Rx = Rx(1:part1+part2,part1+1:part1+part2+part3);

R = triu(qr(Rx));
R = R(n+1:Npara,n+1:Npara+1);
Rnew = triu(qr([R ; Rold]));
Rnew = R(1:Nparaa,:);
RHS = Rnew(:,1:Nparaa);
LHS = Rnew(:,Nparaa+1);
para = RHS\LHS;

B = para(1:n*m);
B = reshape(B,n,m);
para(1:n*m) = [];
D = reshape(para,L,m);

if (mp >= 1)
    B = B + K * D;
end



