function [P,sigma] = dvar4varxbatch(u,y,p,VARX,Zps)
%DVAR4VARX Asymptotic variance of the VARX estimation
%  P=dvar4varx(u,y,p,VARX,Zps) returns the covariance of the VARX
%  estimation and acts as a pre-processor for dvar2frd. The latter is used
%  to calculate the probalistic error bounds around the identified bode
%  diagrams. The input matrix u, output matrix y, and data matrix Zps 
%  (obtained from DORDVARX) must have the same number of observations N. 
%  p is the past window used in DORDVARX and VARX the regression model 
%  obtained from DORDVARX.
%
%  [P,sigma]=dvar4varx(u,y,p,VARX,Zps) also returns the covariance matrix
%  of the innovation noise. 

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check number if input arguments
if nargin < 5
    error('DVAR4VARX requires five input arguments.')
end

% check presence and proper structure of data batches
if iscell(y) && iscell(u)
    if ~iscell(u)
         error('Output y is batch data. Input u must also be batch data.')
    end
    if length(y) ~= length(u)
        error('The amount of batches for u and y must be the same.')
    end
    batch = length(y);
    yb = y;
    ub = u;  
elseif ~iscell(y) && ~iscell(u)
        batch = 1;
else
    error('If output y is batch data, input u must also be batch data, and visa versa')
end
for k = 1:batch
    if batch>1
        y = yb{k};
        u = ub{k};
    end
        
    % check dimensions of inputs   
    if size(y,2) < size(y,1)
        y = y';
    end
    
    N = size(y,2);
    l = size(y,1);
    
    if isempty(u);
        r = 0;
        u = zeros(0,N);
    else
        if size(u,2) < size(u,1)
            u = u';
        end
        r = size(u,1);
        if ~isequal(N,length(u))
            error('The number of samples in u and y must be the same.')
        end
    end
    
    if k == 1
        l_check = l;
        r_check = r;
        P = zeros(l*size(Zps,2),l*size(Zps,2)); % reserve memory for P = lq x lq
    else
        if l_check~=l
            error('The number of outputs must be equal for all batches.')
        elseif r_check~=r
            error('The number of inputs must be equal for all batches.')
        end
    end
    
    if l == 0
        error('DVAR4VARX requires an output vector y.')
    end
    
    % store the past and future vectors
    m = r+l;
    z = [u; y];
    Z = zeros(p*m,N-p);
    for i = 1:p
        Z((i-1)*m+1:i*m,:) = z(:,i:N+i-p-1);
    end

    % solve VARX problem
    y = y(:,p+1:N);
    u = u(:,p+1:N);
    if size(VARX,2)/m > p
        Z = [Z; u];
    end

    % calculate the innovation sequence
    E = y - VARX*Z;
    sigma = (E*E')/(N-p);

    % (sum of) asymptotic variance matrices (for each batch), compute
    % sample-by-sample for memory efficiency
    for j = 1:size(Zps,1)
        P = P + kron(Zps(j,:)',eye(l))*sigma*kron(Zps(j,:),eye(l));        
    end
end

% Calculate average of asymptotic variance matrices
P = P/batch;

