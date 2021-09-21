
%
% Plot a figure of an integral of the function along a lattitudinal contour upon
% the sphere, as the polar angle defining the contour is changed.
% f = figure('visible','off');
% sz = 25;
% ax = gca;
% set(gca,'FontSize',22)

% Colors = ax.ColorOrder;
% hold on

filelabel = 'scalar';
[xi, xi2, O] = computeCompression(N, 1, filelabel);

% Test to make sure I'm not rounding more than I think
if (imag(xi2) > 1e-10) | (imag(xi) > 1e-10)

    fprintf('Warning, rounding off non-negligible imaginary components!')

end

overwriteFlag = false;
heightGridSize = 5;
[varZ, Z, localStateVectors] = getLocalStatesConOpt(N, xi, xi2, heightGridSize, overwriteFlag, filelabel);

for i = 1:length(Z)-1

    % Lifts the minimal eigenvector to the 3N^2 dimensional space
    indexChosen = i;
    Y = localStateVectors(:,indexChosen)' * A;

    % Computes the polar angle of the chosen height
    thetaStar = acos(Z(i))

    % Changes representation of the column vector into three matrices. Note
    % that these matrices are still a representation within the real vector space.
    [Rx1, Rx2, Rx3] = computeMatricesFromVector(Y);

    % Maps the three real matrices to complex hermitian matrices,
    x1 = computeComplexMatrixFromReal(Rx1);
    x2 = computeComplexMatrixFromReal(Rx2);
    x3 = computeComplexMatrixFromReal(Rx3);

    [J1, J2, J3, Jminus, Jplus] = computeSu2(N);

    % Computes the scalar field matrix, or radial projection, of these 
    phi = (1/2) * (J1 * x1 + x1 * J1 + J2 * x2 + x2 * J2 + J3 * x3 + x3 * J3);

    sphHarmonicMatrices = getSphericalHarmonicMatrices(N, 'overwrite');

    c = getMatrixExpansionCoefficients(phi, sphHarmonicMatrices);

    polarGrid = linspace(0,pi,25);
    azimuthGrid = linspace(0,2*pi,25);

    sphHarmonics = getSphericalHarmonics(N, polarGrid, azimuthGrid, 'overwrite');
    
    field = getFieldFromCoefficients(c, sphHarmonics);

    % Plotting a 3D figure that shows explicitly where the field is supported on the sphere, and the positioning of the spherical cap.
    f = figure('visible','off');
    hold on


    % thetaStar = acos(Z(indexChosen));
    % plot(polarGrid, sum(abs(field),2), 'color', Colors(i,:), 'LineWidth', 0.8);
    % xline(thetaStar,'-',{'\theta^*'}, 'color', Colors(i,:),'LineWidth', 1);
    
    % Redefinition of theta, to conform with MATLAB's conventions for spherical coordinates.
    % theta = linspace(pi/2, -pi/2,   length(polarGrid));

    % % Returns the coordinates of a sphere in cartesian
    [XSphere,YSphere,ZSphere]= sphere(length(field)-1);  

    % % Plots the field values upon the sphere. 
    % s = surf(XSphere,YSphere,ZSphere, abs(flip(field,1)), 'FaceAlpha','0.5');
    % s.EdgeColor = 'none';
    % colorbar

    [X,Y] = meshgrid(1:0.5:10,1:20);
    Z = sin(X) + cos(Y);
    surf(X,Y,Z)

    % Plot an ring indicating the position of theta_*
    % azimuthal_line = 0:0.1:2*pi;
    % X_line = sin(azimuthal_line).*cos(thetaStar);
    % Y_line = sin(azimuthal_line).*sin(thetaStar);
    % Z_line = ones(size(X_line))*Z(indexChosen);
    % plot3(X_line, Y_line, Z_line,'.')

    hold off


    % Saves the figure.
    print(f, strcat('savedScalarPolarDependence/fieldSupport3D_', num2str(N),'_Z_', num2str(thetaStar), '.eps'),'-depsc');
    
    print(f, strcat('savedScalarPolarDependence/fieldSupport3D_', num2str(N),'_Z_', num2str(thetaStar), '.png'),'-dpng');

end
% hold off
% xlabel('$\theta$','interpreter','latex')
% ylabel('$f_{\theta^*}$','interpreter','latex', 'Rotation', 0)

% % Creates a filename for the figure to save
% figureFilename = strcat('figs/savedScalarPolarDependence_', filelabel, '.eps');

% % Prints a vectorized eps image of the figure
% print(f,figureFilename,'-depsc');