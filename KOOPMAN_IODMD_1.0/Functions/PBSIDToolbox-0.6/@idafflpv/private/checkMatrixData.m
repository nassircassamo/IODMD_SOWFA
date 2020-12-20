function M = checkMatrixData(M,MatrixName)
% Checks A,B,C,D,K,E data is of proper type

if ~isnumeric(M)
    error('Invalid model: %s matrix must be a numeric array.',MatrixName)
elseif MatrixName~='D' && ~all(isfinite(M(:)))
    error('Invalid model: %s matrix contains Inf or NaN.',MatrixName)
else
    % Convert to full double
    if issparse(M)
        warning('Sparse matrix %s converted to full.',MatrixName)
        M = full(M);
    end
    M = double(M);
end
