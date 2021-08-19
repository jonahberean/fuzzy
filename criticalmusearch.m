function [delta, ZBar] = criticalmusearch(N, nPoints, criticalMu)
% criticalmusearch.m - Uses critical mu values to search for min(delta) k's
%
% Syntax:  [delta,ZBar] = criticalmusearch(N)
%
% Inputs:
%    N          - coordinate matrix dimension
%    nPoints    - number of random points ~ number of subspace vectors
%                 computed
%    criticalMu - mu values at which to perform the generation of subspace
%                 vectors
%
% Outputs:
%    delta   - array of spread of the optimized vectors
%    Zbar    - array of heights of the optimized vectors
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Jonah Berean-Dutcher
% email: jbd@phas.ubc.ca
% March 2021; Last revision: 16-March-2021
%------------- BEGIN CODE --------------

    % Producing projection of Z into Ker(K)
    M = load([fileparts(pwd), '/savedcouplingmatrices/',num2str(N), '.mat']);
    [V, D] = eig(real(M.K));
    V = V(:, abs(diag(D)) < 1e-4);
    xi  = V' * heightoperator(N)     * V;
    xi2 = V' * (heightoperator(N))^2 * V;

    % round the reduced xi and xi2 matrices until they are symmetric
    roundToDigit = 18;
    tf = false;
    while (~tf)
        xi = round(xi, roundToDigit);
        xi2 = round(xi2, roundToDigit);
        tf = ishermitian(xi) & ishermitian(xi2);
        roundToDigit = roundToDigit - 1;
    end
    
    % generate array of points on the quarter circle, in the positive
    % quadrant, for generating vectors in the 2D degenerate subspace
    theta  = linspace(0,pi/2, nPoints)';
    points = [cos(theta), sin(theta)];

    
    % initialize arrays to hold delta, Zbar, but separated by critical mu
    delta     = zeros(nPoints, length(criticalMu));
    ZBar      = zeros(nPoints, length(criticalMu));
   
    % iterate over the critical mu
    for i = 1:length(criticalMu)

        % diagonalize the search matrix, S, at this mu
        [SEigvecs, SDiag] = eig(xi2 - criticalMu(i)*xi);
        SEigvals = diag(SDiag);
        
        % retrieve the lowest three eigenvalues
        SEigvalsDegen = SEigvals(1:3);
        
        % Instead of generating random points, we just generate equally
        % spaced points in the angular coordinate of the quarter circle in
        % the positive quadrant
        % retrieve the two eigenvalue indices we need.
        [~, minI] = min(SEigvalsDegen);
        [~, maxI] = max(SEigvalsDegen);
        
        % !!! It should be possible to remove this loop by vectorizing
        for j = 1:length(points)
            
            % generate subspace vector
            subspaceVec = ...
                [SEigvecs(:, minI), SEigvecs(:, maxI)] * points(j,:)';
            
            % compute height and delta
            ZBar(j,i) = subspaceVec' * xi * subspaceVec;
            delta(j,i) = subspaceVec' * xi2 * subspaceVec - ...
                (subspaceVec' * xi * subspaceVec)^2;
            
        end
        
    end
    
end
