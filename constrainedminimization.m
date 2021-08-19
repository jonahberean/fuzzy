function [delta, Z, kArr] = constrainedminimization(N, numHeights, plotFlag, saveFlag, runFlag, projFlag)

%{
constrainedMin.m - Constrained minimization on the spread problem

Inputs:
   N          - coordinate matrix dimension
   numHeights - number of height values to optimize at
   plotFlag   - boolean to plot or not
   saveFlag   - boolean to save results or not
   runFlag    - boolean to run or if false, just load and plot existing results
   projFlag   - string flag that determines what subspace the height operator will be compressed into
                    - 'scalar' - the subspace associated with the scalar field
                    - 'nullK'  - the null space of the K matrix
                    - 'none'   - no projection at all 

Outputs:
   delta     - array of spread of the optimized vectors
   Z         - array of heights of the optimized vectors
   kArr - array of the optimized vectors themselves

Other m-files required: none
Subfunctions: none
MAT-files required: none

Author: Jonah Berean-Dutcher
email: jbd@phas.ubc.ca
March 2021; Last revision: 29-July-2021 
%}

%------------- BEGIN CODE --------------
% Initializes a filename for saving results
saveFilename = ['savedconstrainedminimization/', num2str(N), '_', projFlag, '_numHeights_', num2str(numHeights), '.mat'];

if(runFlag)

    % Generates the height operators in the full 3N^2 dimensional space
    xi = heightoperator(N);
    xi2 = xi^2;

    if (strcmp(projFlag, 'scalar'))

        % Returns the projection matrix to be used in compressing the height operators
        filename = ['savedscalarfieldprojection/',num2str(N), '.mat'];
        if ~isfile(filename)
            A = scalarfieldprojection(N);
        else
            M = load(filename);
            A = M.A;
        end

        % Compress the height operators into the scalar field subspace
        xi  = A * xi    * A';
        xi2 = A * xi2    * A';

    elseif(strcmp(projFlag, 'nullK'))

        M = load(['savedcouplingmatrices/',num2str(N), '.mat']);
        [V, D] = eig(real(M.K));
        V = V(:, abs(diag(D)) < 1e-4);
        xi  = V' *  xi    * V;
        xi2 = V' * xi2 * V;

    end

    % Rounds the xi and xi2 matrices until they are hermitian
    fprintf('Rounding the xi and xi2 matrices.\n');
    roundToDigit = 18;
    tf = false;
    while (~tf)
        xi = round(xi, roundToDigit);
        xi2 = round(xi2, roundToDigit);
        tf = ishermitian(xi) & ishermitian(xi2);
        roundToDigit = roundToDigit - 1;
    end   
    
    % define array of height bins to optimize in, the defined heightMax is given to be the upper limit of the height operator's spectral range
    heightMin = 0.01;
    heightMax = max(eig(xi));
    height = linspace(heightMin,heightMax,numHeights);
    heightBinWidth = height(2) - height(1);

    % initialize arrays to hold spread, Z, and the vectors, k 
    delta     = zeros(1,         length(height));
    Z         = zeros(1,         length(height));
    kArr = zeros(length(xi),length(height));

    % define an intial guess as a normalized vector of ones
    guess = ones(1,length(xi));
    x0.k = guess / norm(guess);

    % iterate over each height bin
    for i = 1:length(height)

        % print status to console
        statement = '\nMinimizing in bin = %d out of %d\n';
        fprintf(statement, i, length(height))    

        % define the optimization problem
        prob = optimproblem("ObjectiveSense","min");
        k = optimvar('k',length(xi));
        prob.Objective = k'*xi2*k;
        options = optimoptions(@fmincon, 'MaxFunctionEvaluations',1e10);

        % define the normalization constraint
        prob.Constraints.cons1 = k'*k == 1;

        % define the height bin constraint
        prob.Constraints.cons2 = k'*xi*k >= height(i);
        prob.Constraints.cons3 = ...
            k'*xi*k <= height(i) + heightBinWidth;

        % solve the problem
        [sol, ~, exitflag, ~] = solve(prob, x0, 'Options', options);

        % retrieve the optimised k vector
        kOpt = vertcat(sol.k);

        % if optimization was successful, record the result
        if (exitflag == 1)
            
            delta(i) = kOpt' * xi2 * kOpt - (kOpt' * xi * kOpt)^2;
            Z(i)  = kOpt' * xi  * kOpt;
            kArr(:,i)   = kOpt;

        % if optimisation was not successful, record nan (and just height value for Z)    
        else

            delta(i) = nan;
            Z(i) = height(i);
            kArr(:,i) = nan(size(kOpt));
            
        end

    end

    % save the results
    if(saveFlag)        
        save(saveFilename,'Z', 'delta', 'kArr');
    end

else

    % load existing results
    M = load(saveFilename);
    Z = M.Z;
    delta = M.delta;
    kArr = M.kArr;

end

% % plot the results
if(plotFlag)
    % Plots a figure of delta v. Z
    f = figure('visible','off');

    sz = 25;
    scatter(Z, delta, sz, 'filled');

    xlabel('$\bar{Z}$','interpreter','latex')
    ylabel('$\Delta$','interpreter','latex')

    if (strcmp(projFlag, 'scalar'))

        projText = 'scalar';

    elseif(strcmp(projFlag, 'nullK'))

        projText = 'Null($K$)';

    elseif(strcmp(projFlag, 'none'))

        projText = 'full space';

    end

    legend(['minimization in ', projText],'interpreter','latex', 'Location', 'east')

    ax = gca;
    set(gca,'FontSize',22)

    saveas(f, strcat('savedconstrainedminimization/figures/minEigvecs_', num2str(N),'_', projFlag, '_numHeights_', num2str(numHeights), '.png'))
end

