function [xi, xi2, O] = computeCompression(N, K, filelabel)

if(strcmp(filelabel, 'scalar'))

    % Retrieves scalar field projection matrix
    O = getScalarProjectionMatrix(N)';
 
else

    % Diagonalizes coupling matrix
    [V, D] = eig(K);

    if(strcmp(filelabel, 'kernel'))

        O = V(:, abs(diag(D)) < 1e-4);

    elseif strcmp(filelabel, 'kernelComplement')

        O = V(:, abs(diag(D)) > 1e-4);

    end

end


% Retrieves height operator
Z = computeHeightOperator(N);

% Computes compression
xi = O' * Z * O;
xi2 = O' * Z^2 * O;

% Remove the negligible imaginary numerical components
xi = real(xi);
xi2 = real(xi2);

% Rounds the xi and xi2 matrices until they are hermitian
roundToDigit = 18;
tf = false;
while (~tf)
    xi = round(xi, roundToDigit);
    xi2 = round(xi2, roundToDigit);
    tf = ishermitian(xi) & ishermitian(xi2);
    roundToDigit = roundToDigit - 1;
end   