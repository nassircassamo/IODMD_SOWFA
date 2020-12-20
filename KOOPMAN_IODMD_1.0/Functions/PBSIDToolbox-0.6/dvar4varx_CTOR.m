function [P,sigma] = dvar4varx_CTOR(u,y,p,VARX,Zps)
%DVAR4VARX Asymptotic variance of the VARX estimation
%  P=dvar4varx(u,y,p,VARX,Zps) returns the covariance of the VARX
%  estimation and acts as a pre-processor for dvar2frd. The latter is used
%  to calculate the probalistic error bounds around the identified bode
%  diagrams. The input matrix u, output matrix y, past window p, regression
%  model VARX, and data matrix Zps (obtained from dordvarx) must have the
%  same number of observations.
%
%  [P,sigma]=dvar4varxc(u,y,p,VARX,Zps) also returns the covariance matrix
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
if iscell(y)
    if ~iscell(u)
         error('Output y is batch data. Input u must also be batch data.')
    end
    if length(y) ~= length(u)
        error('The amount of batches for u and y must be the same.')
    end
    batch = length(y);
    yb = y;
    ub = u;
else
    if iscell(u)
        error('Input u is batch data. Output y must also be batch data.')
    else
        batch = 1;
    end
end
P = cell(1,batch);
for i = 1:batch
    clear z Z Y U E N sigma
    if batch > 1
        y = yb{i};
        u = ub{i};
    end
    
    % check dimensions of inputs
    if size(y,2) < size(y,1)
        y = y';
    end
    N = size(y,2);
    
    % check consistency of outputs for consecutive batches
    if i > 1
        l_check = size(y,1);
        if l_check ~= l
            error('The number of outputs must be equal for all batches.')
        else
            clear l_check
        end
    end
    l = size(y,1);
    
    if isempty(u);
        r = 0;
        u = zeros(0,N);
    else
        if size(u,2) < size(u,1)
            u = u';
        end
        
        % check consistency of inputs for consecutive batches
        if i > 1
            r_check = size(u,1);
            if r_check ~= r
                error('The number of inputs must be equal for all batches.')
            else
                clear r_check
            end
        end
        r = size(u,1);
        if ~isequal(N,length(u))
            error('The number of rows of vectors/matrices u and y must be the same.')
        end
    end
    if l == 0
        error('DVAR4VARX requires an output vector y.')
    end
    
    % store the past and future vectors
    m = r+l;
    z = [u; y];
    Z = zeros(p*m,N-p);
    for j = 1:p
        Z((j-1)*m+1:j*m,:) = z(:,j:N+j-p-1);
    end
    
    % solve VARX problem
    Y = y(:,p+1:N);
    U = u(:,p+1:N);
    if size(VARX,2)/m > p
        Z = [Z; U];
    end
    
    % calculate the innovation sequence
    E = Y - VARX*Z;
    sigma = (E*E')/length(E);
    
    % asymptotic variance
    P{i} = zeros(l*size(Zps,2),l*size(Zps,2));
    for j = 1:size(Zps,1)
        P{i} = P{i} + kron(Zps(j,:)',eye(l))*sigma*kron(Zps(j,:)',eye(l))';
    end
end

% average covariance matrices obtained for all data batches
P_new = zeros(size(P{1},1),size(P{1},2));
for i = 1:batch
    P_new = P_new + P{i};
end
P_new = P_new./batch;
clear P
P = P_new;
clear P_new
 