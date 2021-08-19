function [Z, delta, subspaceVectors, criticalMuValues] = musearch(N, xi, xi2, numMuValues, numSubspaceVectors)

%region - doc
%{
musearch.m - Returns Z, delta, and the vectors themselves from a search routine upon the compressed heigh operators.

Inputs:
    N   - coordinate matrix dimension
    xi  - compressed height operator
    xi2 - compression of the square of the height operator
    numMuValues - number of values of mu to search over when looking for
    critical points
    numSubspaceVectors - number of vectors to generate in each subspace

Outputs:
    Z               - height of vectors
    delta           - spread of vectors
    subspaceVectors - array containing all vectors computed

Other m-files required: none
Subfunctions: none
MAT-files required: none 

Author: Jonah Berean-Dutcher
email: jbd@phas.ubc.ca
July 2021; Last revision: 18-Aug-2021 
%}
%endregion - doc

fprintf('Generating minimal eigenvectors.\n')

% Initializes the range of mu values to be searched over
muValues            = linspace(-1, 1, numMuValues);

% Initializes an array to hold variances of successive spectra during the search
lowestEigvaluesVar  = zeros(1, length(muValues));

% Iterates over muValues and computes the variance of the bottom three elements of the spectrum of the search matrix
for i = 1:length(muValues)

    eigenvalues = eig(xi2 - muValues(i)*xi);
    lowestEigvaluesVar(i) = var(eigenvalues(1:3));

end

% By finding the local minimums in variance, we identify the critical values of mu where the spectrum has intersections.
criticalMuValues = muValues(islocalmin(lowestEigvaluesVar));

% initializes an array to hold all of the subspace vectors (minimal eigenvectors) generated.
subspaceVectors = zeros(length(xi), length(criticalMuValues), numSubspaceVectors);


% Our search for vectors within each subspace can be limited to a search over the quarter circle, and still produce the full range of results.
theta  = linspace(0.01,pi/2-0.01, numSubspaceVectors)';
points = [cos(theta), sin(theta)];


% Initializes an array to hold the integrated projection operator.
P = zeros(length(xi));

% Initializes arrays to hold delta, Zbar, but separated by critical mu
delta = zeros(length(criticalMuValues), size(points,1));
Z     = zeros(length(criticalMuValues), size(points,1));

% Iterates over the critical values of mu, and generates vectors in the relevant subspace.
for i = 1:length(criticalMuValues)

    [V,D] = eig(xi2 - criticalMuValues(i)*xi);
    eigenvalues = diag(D);

    % Retrieves the two eigenvalue indices we need.
    [~, minI] = min(eigenvalues(1:3));
    [~, maxI] = max(eigenvalues(1:3));
        
        for j = 1:size(points,1)
            
            % Generates a subspace vector.
            subspaceVectors(:, i, j) = [V(:, minI), V(:, maxI)] * points(j,:)';
            
            % Computes height and delta.
            Z(i,j) = subspaceVectors(:, i, j)' * xi * subspaceVectors(:, i, j);
            delta(i,j) = ...
                subspaceVectors(:, i, j)' * xi2 * subspaceVectors(:, i, j) ...
                - (Z(i,j))^2;

        end

end

%     % save the results
%     if(saveFlag)        
%         save(saveFilename,'Z', 'delta', 'kArr');
%     end

% else

%     % load existing results
%     M = load(saveFilename);
%     Z = M.Z;
%     delta = M.delta;
%     kArr = M.kArr;

% end

% % % plot the results
% if(plotFlag)
%     % Plots a figure of delta v. Z
%     f = figure('visible','off');

%     sz = 25;
%     scatter(Z, delta, sz, 'filled');

%     xlabel('$\bar{Z}$','interpreter','latex')
%     ylabel('$\Delta$','interpreter','latex')

%     if (strcmp(projFlag, 'scalar'))

%         projText = 'scalar';

%     elseif(strcmp(projFlag, 'nullK'))

%         projText = 'Null($K$)';

%     elseif(strcmp(projFlag, 'none'))

%         projText = 'full space';

%     end

%     legend(['minimization in ', projText],'interpreter','latex', 'Location', 'east')

%     ax = gca;
%     set(gca,'FontSize',22)

%     saveas(f, strcat('savedconstrainedminimization/figures/minEigvecs_', num2str(N),'_', projFlag, '_numHeights_', num2str(numHeights), '.png'))
% end