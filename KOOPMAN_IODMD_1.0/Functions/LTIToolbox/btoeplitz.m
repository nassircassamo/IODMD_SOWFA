function T = btoeplitz(C,R)
% btoeplitz(C,R) creates block-toeplitz matrix with 
% C first block-column, and R first block-row
%
% btoeplitz(C) is a symmetric (or Hermitian) block-Toeplitz matrix
% with C its first block-column

% Rufus Fraanje, April 2005, visit ISVR

if nargin<2,
    [nc,mc]=size(C);
    nbc = nc/mc; % number of block-colums, block-rows,
    if floor(nbc)~=nbc,
        error('Number of columns of C should be multiple of number of rows of C.');
    end;
    T = zeros(nc);
    T(:,1:mc) = C;
    for i=2:nbc,
        T(:,(i-1)*mc+1:i*mc) = [ C((i-1)*mc+1:i*mc,:)'; T(1:(nbc-1)*mc,(i-2)*mc+1:(i-1)*mc)];
    end;
    % Ready here
else
    [nc,mc]=size(C);
    [nr,mr]=size(R);

    % block-size is mr x mc
    % hence, nr should be multiple of nc
    %    and mc should be multiple of mr

    nbc = mr/mc; % number of block-columns,
    if floor(nbc)~=nbc
        error('Number of colums of C should be multiple of number of columns of R.');
    end;
    nbr = nc/nr; % number of block-rows,
    if floor(nbr)~=nbr,
        error('Number of rows of R should be multiple of number of rows of C.');
    end;
    % first block element in C and R should be equal
    if norm(C(1:nr,:)-R(:,1:mc),'fro')>=max(size(C,1),size(R,2))*eps,
        error('First block element in C and R should be equal.');
    end;

    T = zeros(nbr*nr,nbc*mc);
    T(:,1:mc) = C;
    T(1:nr,:) = R;
    if nbr>nbc, % try to get smallest for-loop
        for i=2:nbc,
            T(:,(i-1)*mc+1:i*mc) = [ R(:,(i-1)*mc+1:i*mc); T(1:(nbr-1)*nr,(i-2)*mc+1:(i-1)*mc)];
        end;
    else
        for i=2:nbr,
            T((i-1)*nr+1:i*nr,:) = [C((i-1)*nr+1:i*nr,:) T((i-2)*nr+1:(i-1)*nr,1:(nbc-1)*mc)];    
        end
    end;
    
    % ready here
end;
