% function investigateazimuth(startN, stopN)

%region - doc
%{
investigateazimuth.m - For a single N value, we'll generate a high density of
minimal eigenvectors, and do some tests to try and get at what's happening with
the azimuthal dependency.

Inputs:
   N         - size of the coordinate matrices

Outputs:

Other m-files required: compressedheightoperators.m, musearch.m, 
vectortomatrices.m

Subfunctions: none
MAT-files required: none

Author: Jonah Berean-Dutcher
email: jbd@phas.ubc.ca
Aug 2021; Last revision: 18-Aug-2021 
%}
%endregion - doc

N = 5;

f = figure('visible','off');
% hold on

% [Zbar, delta, subspaceVectors, compressingMatrix, L1, L2, L3, criticalMuValues, numFieldGrid] = scalarfieldprocedure(N);
Npoints = 20;
[delta, Z, kArr] = constrainedminimization(N, Npoints, false, false, false, 'scalar');
% [delta, Z, kArr] = constrainedminimization(N, Npoints, true, true, true, 'scalar');

for i = 1:size(kArr,2)

    % Lifts the minimal eigenvector to the 3N^2 dimensional space
    Y = kArr(:,i)' * compressingMatrix;

    % Changes representation of the column vector into three matrices. Note
    % that these matrices are still a representation within the real vector space.
    [Rx1, Rx2, Rx3] = vectortomatrices(Y, N);

    % Maps the three real matrices to complex hermitian matrices,
    x1 = 1/2*(Rx1 + Rx1.') - 1/(2*sqrt(-1))*(Rx1 - Rx1.');
    x2 = 1/2*(Rx2 + Rx2.') - 1/(2*sqrt(-1))*(Rx2 - Rx2.');
    x3 = 1/2*(Rx3 + Rx3.') - 1/(2*sqrt(-1))*(Rx3 - Rx3.');

    % Computes the scalar field matrix, or radial projection, of these matrices.
    phi = (1/2) * (L1 * x1 + x1 * L1 + L2 * x2 + x2 * L2 + L3 * x3 + x3 * L3);

    % Returns the function that the phi matrix is mapped to.
    [field, theta, phi, c] = matrixtofield(phi, numFieldGrid);

    % Iterating over theta and taking the variance over phi of each
    % lattitudinal slice.
    aziVar = zeros(1,size(field,1));
    for t = 1:size(field,1)

        aziVar(t) = var(field(t,:));

    end

    % Plotting
    if (~isnan(kArr(1,i)))
        aziVarPeak(i) = max(aziVar);
        
    end

   
end
% hold off

scatter(Z(1:19), aziVarPeak, 'DisplayName', num2str(Z(i)))

saveas(f, strcat('figures/aziTestNew2.png'));

   
% aziVarMean = zeros(size(subspaceVectors,2), size(subspaceVectors,3) );

% % Iterate over the subspace vectors of one critical mu
% % fprintf('\nIterating over minimal eigenvectors.\n')
% % for i = 1:size(subspaceVectors,2)     
   
% %    fprintf('criticalMuValue %d of %d\n', i, length(criticalMuValues))

% i = 3;

% for j = 1:size(subspaceVectors,3)

%     fprintf('subSpaceVector %d of %d\n', j, size(subspaceVectors,3))

%     % Lifts the minimal eigenvector to the 3N^2 dimensional space
%     Y = subspaceVectors(:,i,j)' * compressingMatrix;

%     % Changes representation of the column vector into three matrices. Note
%     % that these matrices are still a representation within the real vector space.
%     [Rx1, Rx2, Rx3] = vectortomatrices(Y, N);

%     % Maps the three real matrices to complex hermitian matrices,
%     x1 = 1/2*(Rx1 + Rx1.') - 1/(2*sqrt(-1))*(Rx1 - Rx1.');
%     x2 = 1/2*(Rx2 + Rx2.') - 1/(2*sqrt(-1))*(Rx2 - Rx2.');
%     x3 = 1/2*(Rx3 + Rx3.') - 1/(2*sqrt(-1))*(Rx3 - Rx3.');

%     % Computes the scalar field matrix, or radial projection, of these matrices.
%     phi = (1/2) * (L1 * x1 + x1 * L1 + L2 * x2 + x2 * L2 + L3 * x3 + x3 * L3);

%     % Returns the function that the phi matrix is mapped to.
%     [field, theta, phi, c] = matrixtofield(phi, numFieldGrid);

%     % Removes negligible imaginary component from field
%     if(max(imag(field), [], 'all') > 1e-10)
%         fprintf('Warning, removing non-negligible imaginary component! Magnitude: %.10f', max(imag(field), [], 'all'))
%     end
%     field = real(field);

%     % Iterating over theta and taking the variance over phi of each
%     % lattitudinal slice.
%     aziVar = zeros(1,size(field,1));
%     for t = 1:size(field,1)

%         aziVar(t) = var(field(:,t));

%     end

%     % Plotting

%     scatter(theta(:,3), sum(abs(field),2), 'DisplayName', num2str(Zbar(i,j)))
%     legend;

%     saveas(f, strcat('figures/aziTest12.png'));

% end

% % end

% hold off
% % Plotting
% f = figure('visible','off');
% hold on
% scatter(Zbar(:), delta(:), 'DisplayName', num2str(N))
% legend;
% hold off
% saveas(f, strcat('figures/aziTest1.png'));