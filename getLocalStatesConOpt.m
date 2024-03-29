function [varZ, Z, localStateVectors, filename] = getLocalStatesConOpt(N, xi, xi2, heightGridSize, overwriteFlag, filelabel)

%region
%{
makeLocalStatesConOpt.m - Computes the localized states on the sphere,
throughout Z.
%}
%endregion - doc

% Checks to see if the file exists, will make if not
filename = ['savedLocalStatesConOpt/',num2str(N), '_heightGridSize_', num2str(heightGridSize), '_', filelabel];
if (isfile([filename,'.mat']) & ~overwriteFlag)
    fprintf(['Loading ', filename, '\n'])
    M                 = load(filename);
    varZ              = M.varZ;
    Z                 = M.Z;
    localStateVectors = M.localStateVectors;
else
    fprintf(['Generating ', filename, '\n'])
    heightGrid = linspace(0.01, max(eig(xi)), heightGridSize);
    % heightGrid = linspace(0.18, 0.22, heightGridSize);
    heightGridSpacing = heightGrid(2) - heightGrid(1);

    % initialize arrays to hold spread, Z, and the vectors
    varZ      = zeros(1,         length(heightGrid));
    Z         = zeros(1,         length(heightGrid));
    localStateVectors = zeros(length(xi),length(heightGrid));

    % define an intial guess as a normalized vector of ones
    guess = ones(1,length(xi));
    x0.k = guess / norm(guess);

    % iterate over each height bin
    for i = 1:length(heightGrid)

        % print status to console
        statement = '\nMinimizing in bin = %d out of %d\n';
        fprintf(statement, i, length(heightGrid))    

        % define the optimization problem
        prob = optimproblem("ObjectiveSense","min");
        k = optimvar('k',length(xi));
        prob.Objective = k'*xi2*k;
        options = optimoptions(@fmincon, 'MaxFunctionEvaluations',1e5);

        % define the normalization constraint
        prob.Constraints.cons1 = k'*k == 1;

        % define the height bin constraint
        prob.Constraints.cons2 = k'*xi*k >= heightGrid(i);
        prob.Constraints.cons3 = ...
            k'*xi*k <= heightGrid(i) + heightGridSpacing;
        % prob.Constraints.cons3 = ...
        %     k'*xi*k <= heightGrid(i) + 0.1;

        % solve the problem
        [sol, ~, exitflag, ~] = solve(prob, x0, 'Options', options);

        % retrieve the optimised k vector
        kOpt = vertcat(sol.k);

        % if optimization was successful, record the result
        if (exitflag == 1)
            
            varZ(i) = kOpt' * xi2 * kOpt - (kOpt' * xi * kOpt)^2;
            Z(i)  = kOpt' * xi  * kOpt;
            localStateVectors(:,i)   = kOpt;

        % if optimisation was not successful, record nan (and just height value for Z)    
        else

            varZ(i) = nan;
            Z(i) = heightGrid(i);
            localStateVectors(:,i) = nan(size(kOpt));
            
        end

    end

    % save the results    
    save([filename,'.mat'],'varZ', 'Z', 'localStateVectors');

end

