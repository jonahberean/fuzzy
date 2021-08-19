function [xiReduced, xi2Reduced] = reducedmatrices(N, xi, xi2, nullEigvecs)
% reducedmatrices - generates symmetry reduced matrices \xi, and \xi_2

% Syntax:  [xiReduced, xi2Reduced] = reducedmatrices(N)
%
% Inputs:
%    N   - coordinate matrix dimension
%
% Outputs:
%    xiReduced  - projected Z operator, further reduced by symmetry
%    xi2Reduced - projected Z^2 operator, further reduced by symmetry
%
% Other m-files required: blkByBlkDiag.m
% MAT-files required: /projections/xi_N.mat, 
%                     !!! clearly these namings can be improved
%                     /couplingMatrixEigenvectors/couplingMatrixEigenvectors_N.mat
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Jonah Berean-Dutcher
% email: jbd@phas.ubc.ca
% March 2021; Last revision: 16-March-2021
%------------- BEGIN CODE --------------

    % set tolerance for null
    tolerance = 1e-4;

    % generate the matrix representation of so(3) symmetry
    so3Mat = symmetrymatrix(N);
    
    % project this matrix into Ker(K)
    so3Mat = nullEigvecs' * so3Mat * nullEigvecs;
    
    % so3Operator is originally antisymmetric, and we want the Hermitian 
    % version
    so3Mat  = 1i*so3Mat;
    roundToDigit = 18;
    tf = false;
    while (~tf)
        so3Mat = round(so3Mat, roundToDigit);
        tf = ishermitian(so3Mat);
        roundToDigit = roundToDigit - 1;
    end

    % diagonalize the projected so3Operator
    [so3MatEigvecs, so3MatDiag] = eig(so3Mat);
    
    % block diagonalize xi and xi2 using the eigenvectors of the projected
    % so3 operator
    xiBlkDiag  = so3MatEigvecs' * xi  * so3MatEigvecs;
    xi2BlkDiag = so3MatEigvecs' * xi2 * so3MatEigvecs;

    % block by block diagonalize xi and xi2 using the diagonal form of the
    % so3 operator
    [xiEigvecs , xiDiag]   = blkbyblkdiag(xiBlkDiag , so3MatDiag);
    [xi2Eigvecs, xi2Diag]  = blkbyblkdiag(xi2BlkDiag, so3MatDiag);

    % retrieve the subspaces of xi and xi2 for which so3Diag has 
    % non-negative eigenvalues
    constraint        = diag(so3MatDiag) > -tolerance;
    xiEigvecsReduced  = xiEigvecs  (constraint, constraint);
    xiDiagReduced     = xiDiag     (constraint, constraint);
    xi2EigvecsReduced = xi2Eigvecs (constraint, constraint);
    xi2DiagReduced    = xi2Diag    (constraint, constraint);

    % construct the reduced matrices
    xiReduced  = xiEigvecsReduced  * xiDiagReduced  * xiEigvecsReduced';
    xi2Reduced = xi2EigvecsReduced * xi2DiagReduced * xi2EigvecsReduced';
    
    % round the reduced xi and xi2 matrices until they are symmetric
    roundToDigit = 18;
    tf = false;
    while (~tf)
        xiReduced = round(xiReduced, roundToDigit);
        xi2Reduced = round(xi2Reduced, roundToDigit);
        tf = ishermitian(xiReduced) & ishermitian(xi2Reduced);
        roundToDigit = roundToDigit - 1;
    end

end