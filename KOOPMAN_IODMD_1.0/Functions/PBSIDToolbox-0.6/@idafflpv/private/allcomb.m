function B = allcomb(A)
    %% all combinations of columns in A

    [ni,nm] = size(A);
    ii = ni:-1:1;
    A = mat2cell(A, ones(1,ni), nm);
    % flip using ii if last column is changing fastest
    [B{ii}] = ndgrid(A{ii}) ;
    % concatenate
    B = reshape(cat(ni+1,B{:}),[],ni)';
end