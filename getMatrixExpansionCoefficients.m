function c = getMatrixExpansionCoefficients(A, sphHarmonicMatrices)

N = length(A);

% Itertatively take the inner product of A with each of these matrices, so
% as to compute the coefficients of A's expansion in this basis.
for j = 0:(N-1)    
    
    mValues = -j:1:j;
    
    for m = 1:length(mValues)
                              
        % computing coefficients and constructing values of function space
        % correspondence
        c(j+1,m) = (1/N) * trace(sphHarmonicMatrices(:,:,j+1,m)' * A);
        
    end
    
end