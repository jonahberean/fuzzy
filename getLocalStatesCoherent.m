function [varZ, Z, localStateVectors] = getLocalStatesCoherent(N, projectionMatrix, xi, xi2, heightGridSize, numVectorsPerHeight, overwriteFlag, filelabel)
    


% Checks to see if the file exists, will make if not
filename = ['savedLocalStatesCoherent/',num2str(N), '_heightGrigSize_', num2str(heightGridSize), '_numVectorsPerHeight', num2str(numVectorsPerHeight), '_', filelabel,'.mat'];
if (isfile(filename) & ~overwriteFlag)
    M                 = load(filename);
    varZ              = M.varZ;
    Z                 = M.Z;
    localStateVectors = M.localStateVectors;
else
            
    % Returns the L3 generator of su(2)
    [~,~,L3,~,~] = computeSu2(N);

    % Generates a list of angles from 0 to pi/2
    listOfAngles = acos(linspace(1,0,heightGridSize));

    % Generates a table of random coefficients
    y = -1 + 2*randn(numVectorsPerHeight,3);
    listOfCoefficients = bsxfun(@rdivide,y,sqrt(sum(y.^2,2)));

    % Initializes arrays to store the Z and delta values
    varZ = zeros(heightGridSize, numVectorsPerHeight);
    Z    = zeros(heightGridSize, numVectorsPerHeight);
    localStateVectors = zeros(length(xi),heightGridSize, numVectorsPerHeight);

    % In the (2j + 1) = N irrep. of su(2), using as a basis the eigenvectors of 
    % L3 angular momentum |m>. The coherent state at the north pole is the 
    % state with the largest angular momentum in the 3-direction, which is 
    % explicitly the N-dimensional vector: (1, 0, 0,..., 0) 
    zKet    = zeros(N,1);
    zKet(1) = 1;

    % Iterates over the number of angles at which we'd like to perform the computation
    for i = 1:heightGridSize

        fprintf('Working with height %d of %d\n', i, heightGridSize)

        % Sets angle value
        beta = listOfAngles(i);

        % Performs a rotation on the state: |n> = D(R)|z>, such that the state 
        % is rotated by an angle \beta about the y-axis. 
        % D implements this rotation - from Wigner formula
        D    = computeWignerMatrix(N, beta);
        nKet = D * zKet;

        % Generates the projection operator |n><n|
        nProjectionOperator = nKet * nKet';

        % Iterates over the coefficients we generated randomly earlier
        for ii = 1:numVectorsPerHeight       

            c = listOfCoefficients(ii,:);

            % Reformulates the projection operator into a column vector (Y) in 
            % the 3N^2 dimensional space, using the c_i coefficients as weights 
            % for the three matrices
            Y = computeVectorFromMatrices(...
                c(1) * nProjectionOperator, ...
                c(2) * nProjectionOperator, ...
                c(3) * nProjectionOperator);
            
            % Project Y 
            Y = projectionMatrix'*Y';
            
            % normalize Y
            Y = Y/norm(Y);

            % compute the delta and the Z
            Z(i, ii)     = Y' * xi  * Y;
            varZ(i, ii) = Y' * xi2 * Y - (Z(i,ii))^2;
            localStateVectors(:,i,ii) = Y;

        end
        
    end

    % save the results    
    save(filename,'varZ', 'Z', 'localStateVectors');
end