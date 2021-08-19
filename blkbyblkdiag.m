function [eigvecs, diagMat] = blkbyblkdiag(blkDiagMat, diagCommutingMat)
% blkbyblkdiag.m - Diagonalizes a block diagonal matrix by way of a 
%                  different, commuting matrix,  

% Syntax:  [eigvecs, diagMat] = blkbyblkdiag(blkDiagMat, diagCommutingMat)
%
% Inputs:
%    blkDiagMat       - Block diagonal matrix
%    diagCommutingMat - Matrix that commutes with the first, given in
%                       diagonal form
%
% Outputs:
%    eigvecs          - resulting eigenvectors
%    diagMat          - resulting diagonal form of the matrix
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: Jonah Berean-Dutcher
% email: jbd@phas.ubc.ca
% March 2021; Last revision: 16-March-2021
%------------- BEGIN CODE --------------

% create matrices to hold the eigenvectors and eigenvalues of the fully
% diagonalized coupling matrix
eigvecs = zeros(size(blkDiagMat));
diagMat = zeros(size(blkDiagMat));

% getting the eigenvalues of the diagonal commuting matrix
diagCommutingEigenvalues = diag(diagCommutingMat);

% get the unique eigenvalues, and their indices
% note: uniquetol takes only real valued input.
[~, ia] = uniquetol(diagCommutingEigenvalues,1e-4);
ia = sort(ia, 'ascend');

% Iterate over each degenerate block, diagonalizing block by block 
for i=1:(length(ia))
    
    % get the index of the block start, and block end
    start = ia(i);
    if(i < length(ia))
        stop  = ia(i+1)-1;
    else
        stop = (length(blkDiagMat));
    end
    
    % diagonalize the block and store the results
    [blockEigvecs, blockDiag] = eig(blkDiagMat(start:stop, start:stop));
    eigvecs(start:stop, start:stop) = blockEigvecs;
    diagMat(start:stop, start:stop) = blockDiag;
    
end
