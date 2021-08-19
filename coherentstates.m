function [delta, Z] = coherentstates(N, numAngles, numCoeffs, plotFlag, saveFlag, runFlag, projFlag)
% coherentstates.m -
%
% Syntax: coherentstates(N)
%
% Inputs:
%     N          - coordinate matrix dimension
%     numAngles  -
%     numCoeffs - 
%     plotFlag   - boolean to plot or not
%     saveFlag   - boolean to save results or not
%     runFlag    - boolean to run or if false, just load and plot existing results
%     projFlag   - string flag that determines what subspace the height operator will be compressed into
%                      - 'nullKT' - the [null(k)]^T
%                      - 'nullK'  - the null space of the K matrix
%                      - 'none'   - no projection at all 
%       
% Outputs:
%
% Other m-files required: wignermatrix.m, matricestovector.m
% Subfunctions:
% MAT-files required: none
%
% Author: Jonah Berean-Dutcher
% email: jbd@phas.ubc.ca
% March 2021; Last revision: 10-Aug-2021

saveFilename = ['savedcoherentstates/', num2str(N), '_numAngles_', ...
num2str(numAngles), '_numCoeffs_',num2str(numCoeffs), '_', projFlag, '.mat'];

if (runFlag)

    % Loads and diagonalizes the K matrix
    M = load(['savedcouplingmatrices/',num2str(N), '.mat']);
    [V, D] = eig(real(M.K));

    if (strcmp(projFlag, 'nullK'))
        
        % Projects into Null(K)
        U = V(:, abs(diag(D)) < 1e-4);
        
    elseif (strcmp(projFlag, 'nullKT'))
        
        % Projects into [Null(K)]^T 
        U = V(:, abs(diag(D)) > 1e-4);
            
    else
        
        % No projection will be carried out
        U = eye(3*N^2);

    end

    % Carries out the projection.
    xi   = U' * heightoperator(N)     * U;
    xi2  = U' * (heightoperator(N))^2 * U;

    % Returns the L3 generator of su(2)
    [~,~,L3,~,~] = su2generators(N);

    % Generates a list of angles from 0 to pi/2
    listOfAngles = linspace(0.01, pi/2, numAngles);

    % Generates a table of random coefficients
    y = -1 + 2*randn(numCoeffs,3);
    listOfCoefficients = bsxfun(@rdivide,y,sqrt(sum(y.^2,2)));

    % Initializes arrays to store the Z and delta values
    delta = zeros(numAngles, numCoeffs);
    Z     = zeros(numAngles, numCoeffs);

    % In the (2j + 1) = N irrep. of su(2), using as a basis the eigenvectors of 
    % L3 angular momentum |m>. The coherent state at the north pole is the 
    % state with the largest angular momentum in the 3-direction, which is 
    % explicitly the N-dimensional vector: (1, 0, 0,..., 0) 
    zKet = zeros(N,1);
    zKet(1) = 1;

    % Iterates over the number of angles at which we'd like to perform the computation
    for i = 1:numAngles

        % Sets angle value
        beta = listOfAngles(i);

        % Prints status
        fprintf('\nbeta = %.2f (%d out of %d)\n', beta, i, numAngles)

        % Performs a rotation on the state: |n> = D(R)|z>, such that the state 
        % is rotated by an angle \beta about the y-axis. 
        % D implements this rotation - from Wigner formula
        D    = wignermatrix(N, beta);
        nKet = D * zKet;

        % Generates the projection operator |n><n|
        nProjectionOperator = nKet * nKet';

        % Iterates over the coefficients we generated randomly earlier
        for ii = 1:numCoeffs       

            c = listOfCoefficients(ii,:);

            % Reformulates the projection operator into a column vector (Y) in 
            % the 3N^2 dimensional space, using the c_i coefficients as weights 
            % for the three matrices
            Y = matricestovector(...
                c(1) * nProjectionOperator, ...
                c(2) * nProjectionOperator, ...
                c(3) * nProjectionOperator, N);
            
            % Project Y 
            Y = U'*Y;
            
            % normalize Y
            Y = Y/norm(Y);

            disp(size(xi))
            disp(size(Y))

            % compute the delta and the Z
            Z(i, ii)     = Y' * xi  * Y;
            delta(i, ii) = Y' * xi2 * Y - (Z(i,ii))^2;

        end
        
    end

    % save the results
    if(saveFlag)        
        save(saveFilename,'Z', 'delta');
    end

else
    % load existing results
    M = load(saveFilename);
    Z = M.Z;
    delta = M.delta;

end

% Generating plots.
if(plotFlag)
    % Plots a figure of delta v. Z
    f = figure('visible','off');
    set(gca,'FontSize',22) 
    hold on
    for i = 1:numAngles

        sz = 25;
        scatter(Z(i,:), delta(i,:), sz, 'filled');

    end
    hold off

    xlabel('$\bar{Z}$','interpreter','latex')
    ylabel('$\Delta$','interpreter','latex')

    saveas(f, strcat('savedcoherentstates/figures/minEigvecs_', num2str(N), '_numAngles_', ...
    num2str(numAngles), '_numCoeffs_',num2str(numCoeffs), '_', projFlag, '.png'))
end



% Testing.
if (runFlag)
    % Computes j, the highest weight of the rep.
    j = (1/2)*(N-1);

    % Generates a vector pointing in the z-direction
    z = [0; 0; 1];

    % Retrieves the su2 generator matrices for this N
    [L1, L2, L3, ~, ~] = su2generators(N);

    % 3D rotation matrix about the y-axis
    R = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];

    % Rotates the vector z
    n = R*z; 

    % Generates the L_n operator
    LnOperator = n(1)*L1 + n(2)*L2 + n(3)*L3;

    % Evaluating whether or not the correct eigenvalue property is observed
    coherentStateTest = norm(abs((LnOperator*nKet) - j*nKet));
    if (coherentStateTest > 0.001)

        disp('Warning, L_n|n> ~= j|n>')
        disp(coherentStateTest)

    end

    % Evaluating whether or not the coherent state is properly normalized
    unity = nKet'*nKet;
    if(abs(unity - 1) > 0.001)

        disp('Warning, <n|n> ~= 1')
        disp(unity)

    end
end

