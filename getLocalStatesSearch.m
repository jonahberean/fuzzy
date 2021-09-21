function [varZ, Z, localStateVectors, filename] = getLocalStatesSearch(N, xi, xi2, heightGridSize, overwriteFlag, filelabel)

% Checks to see if the file exists, will make if not
filename = ['savedLocalStatesSearch/',num2str(N), '_heightGridSize_', num2str(heightGridSize), '_', filelabel];
if (isfile([filename,'.mat']) & ~overwriteFlag)
    fprintf(['Loading ', filename, '\n'])
    M                 = load(filename);
    varZ              = M.varZ;
    Z                 = M.Z;
    localStateVectors = M.localStateVectors;

else

    % Initializes the range of mu values to be searched over
    muValues = linspace(0, 1.5, heightGridSize);

    % Initializes an array to hold variances of successive spectra during the search
    lowestEigvaluesVar = zeros(1, length(muValues));

    % Iterates over muValues and computes the variance of the bottom three elements of the spectrum of the search matrix
    for i = 1:length(muValues)

        eigenvalues = eig(xi2 - muValues(i)*xi);
        lowestEigvaluesVar(i) = var(eigenvalues(1:3));

    end

    numLocalStateVectors = 50;

    % By finding the local minimums in variance, we identify the critical values of mu where the spectrum has intersections.
    criticalMuValues = muValues(islocalmin(lowestEigvaluesVar));

    % initializes an array to hold all of the subspace vectors (minimal eigenvectors) generated.
    localStateVectors = zeros(length(xi), length(criticalMuValues), numLocalStateVectors);

    % Our search for vectors within each subspace can be limited to a search over the quarter circle, and still produce the full range of results.
    theta  = linspace(0.01,pi/2-0.01, numLocalStateVectors)';
    points = [cos(theta), sin(theta)];

    % Initializes arrays to hold varZ, Zbar, but separated by critical mu
    varZ = zeros(length(criticalMuValues), size(points,1));
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
                localStateVectors(:, i, j) = [V(:, minI), V(:, maxI)] * points(j,:)';
                
                % Computes height and varZ.
                Z(i,j) = localStateVectors(:, i, j)' * xi * localStateVectors(:, i, j);
                varZ(i,j) = ...
                    localStateVectors(:, i, j)' * xi2 * localStateVectors(:, i, j) ...
                    - (Z(i,j))^2;

            end

    end

    Z = Z(:);
    varZ = varZ(:);

    % save the results    
    save([filename,'.mat'],'varZ', 'Z', 'localStateVectors');
end