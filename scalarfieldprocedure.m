function [Zbar, delta, subspaceVectors, compressingMatrix, L1, L2, L3, criticalMuValues, numFieldGrid] = scalarfieldprocedure(N)

%region - doc
%{
scalarfieldprocedure.m - Runs through various steps to associate a function to
the degrees of freedom of a scalar field

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

% Initializing parameters for testing
% N = 6;
numMuValues = 10000;
numSubspaceVectors = 6;
numFieldGrid = 50;

% Returns generators of the su(2) Lie algebra
[L1, L2, L3, ~, ~] = su2generators(N);

% rescaling the su(2) generators so as to normalize the fuzzy sphere radius to 1.
nu     = 2 / sqrt(N^2 - 1);
L1 = nu*L1;
L2 = nu*L2;
L3 = nu*L3;

% Returns compressed height operators for this value of N
[compressingMatrix, xi, xi2] = compressedheightoperators(N, 'scalar');

% Returns minimal eigenvectors, with their heights and variances for a given
% compressed matrix.
[Zbar, delta, subspaceVectors, criticalMuValues] = musearch(...
                        N, xi, xi2, numMuValues, numSubspaceVectors);

% % Iterating over the minimal eigenvectors, each is lifted into it's
% % representation in the full 3N^2 dimensional theory space.
% fprintf('\nIterating over minimal eigenvectors.\n')
% for i = 1:length(criticalMuValues)     
   
%    fprintf('criticalMuValue %d of %d\n', i, length(criticalMuValues))

%    for j = 1:numSubspaceVectors

%       fprintf('subSpaceVector %d of %d\n', j, numSubspaceVectors)

%       % Lifts the minimal eigenvector to the 3N^2 dimensional space
%       Y = subspaceVectors(:,i,j)' * compressingMatrix;

%       % Changes representation of the column vector into three matrices. Note
%       % that these matrices are still a representation within the real vector space.
%       [Rx1, Rx2, Rx3] = vectortomatrices(Y, N);

%       % Maps the three real matrices to complex hermitian matrices,
%       x1 = 1/2*(Rx1 + Rx1.') - 1/(2*sqrt(-1))*(Rx1 - Rx1.');
%       x2 = 1/2*(Rx2 + Rx2.') - 1/(2*sqrt(-1))*(Rx2 - Rx2.');
%       x3 = 1/2*(Rx3 + Rx3.') - 1/(2*sqrt(-1))*(Rx3 - Rx3.');

%       % Computes the scalar field matrix, or radial projection, of these matrices.
%       phi = (1/2) * (L1 * x1 + x1 * L1 + L2 * x2 + x2 * L2 + L3 * x3 + x3 * L3);

%       % Returns the function that the phi matrix is mapped to.
%       [field, theta, phi, c] = matrixtofield(phi, numFieldGrid);

%       % Removes negligible imaginary component from field, and takes absolute magnitude.
%       if(max(imag(field), [], 'all') > 1e-10)
%          fprintf('Warning, removing non-negligible imaginary component! Magnitude: %.10f', max(imag(field), [], 'all'))
%       end
%       field = abs(real(field));

%       % % Redefinition of theta, to conform with MATLAB's conventions for spherical coordinates.
%       % theta = linspace(pi/2, -pi/2,   numFieldGrid);

%       % % Returns the coordinates of a sphere in cartesian
%       % [X,Y,Z]= sphere(length(field)-1);  

%       % % Plotting a 3D figure that shows explicitly where the field is supported on the sphere, and the positioning of the spherical cap.
%       % f = figure('visible','off');

%       % % Plots the field values upon the sphere. 
%       % s = surf(X,Y,Z, field, 'FaceAlpha','0.5');
%       % s.EdgeColor = 'none';
%       % colorbar

%       % % Saves the figure.
%       % saveas(f, strcat('figures/fieldSupport3D_', num2str(N), '_Z_',
%       % num2str(Zbar(i,j)), '.png'));
      
%       % Plotting figures that present the results along either theta or phi
%       f = figure('visible','off');

%       % Computing the vector angle from the height
%       vecAngle = acos(Zbar(i,j));

%       % For each lattitudinal slice, we sum over absolute field values at all azimuthal values. 
%       v = sum(abs(field), 2);
%       v = v/norm(v);
%       hold on
%       scatter(theta(:,1), v);
%       % Plotting a line indicating the edge of the cap.
%       xline(vecAngle,'-',{'\theta^*'});
%       hold off
%       saveas(f, strcat('figures/fieldSupport_Theta_', ...
%       num2str(N), '_', num2str(vecAngle), '.png'))

%       % Plotting figures that present the results along either theta or phi
%       f = figure('visible','off');

%       % For each lattitudinal slice, we sum over absolute field values at all azimuthal values. 
%       v = sum(abs(field), 1);
%       v = v/norm(v);
%       scatter(phi(2,:), v);
%       saveas(f, strcat('figures/fieldSupport_Phi_', ...
%       num2str(N), '_', num2str(vecAngle), '.png'))
      

%    end

% end


