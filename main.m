%{
main.m - scripting space for development and running large processes

Other m-files required: varies
Subfunctions: varies
MAT-files required: varies

Author: Jonah Berean-Dutcher
email: jbd@phas.ubc.ca
July 2021; Last revision: 04-Aug-2021 
%}

% AUG. 17 TODO

% Our entire setup appears to be working now, and generating fully real functions. I think first we should work on generating a few first attempts at plotting these real functions, and looking at the result. I want to be able to ask Joanna:

% What are we meant to do with this? We think we're in a sector? 

% shit, we still haven't demonstrated/resolved the whole not-actually-full-rank minimal eigenvector problem. 

% 1) Run this proc. and plot fields.

% 2) Revisit the minimal eigvec problem

% 3) Draft the question about sectors, and what the way forward is for Joanna.

% 4) Meet with Joanna

% 5) Just talk prep from there on.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% N = 4;
% A = complex(rand(N, N), rand(N,N));
% A = (1/2)*(A + A');
% numFieldGrid = 25;
% [f, theta, phi, c, expansionTest, Y, Jplus, Jminus, Ji] = matrixtofield(A, numFieldGrid);

% j = 1;
% p = 2;
% sqrt(p*(2*j+1-p))

% n = 3;
% m = 1;
% th = pi/2+0.01;
% phi = pi/2+0.01;
% Y = harmonicY(n,m,th,phi);
% Y1 = harmonicY(n,-m,th,phi);

% Now we have correctly normalized our su(2) generators, and subsequently our coupling matrices (as the potential is a function of those generators) as well as our height operator of course. 
% [X,Y,Z]= sphere(5)
% [azimuth,elevation,r] = cart2sph(X,Y,Z)

% N = 5;
% numMuValues        = 1000;
% numSubspaceVectors = 3;
% numFieldGrid       = 100;

% [L1, L2, L3, ~, ~] = su2generators(N);

% % Returns compressions of the Z and Z^2 operators into the subspace associated with the scalar field, within the non-null sector.
% fprintf('Generating compressed height operators.\n')
% [A, xi, xi2] = compressedheightoperators(N);
% % xi = heightoperator(N);
% % xi2 = xi^2;
% % Returns Z, delta, and the vectors themselves from a search routine upon the compressed heigh operators.
% fprintf('Generating minimal eigenvectors.\n')
% [Zbar, delta, subspaceVectors, criticalMuValues] = musearch(...
%                         N, xi, xi2, numMuValues, numSubspaceVectors);

% for i = 1:length(criticalMuValues)

%     for j = 1:numSubspaceVectors

%         height = Zbar(i,j);

%         % Produce a 3N^2 dimensional vector that corresponds to the minimal eigenvector from the phi subspace. 

%         Yvec = subspaceVectors(:,i,j)'*A;

%         % Produce the representation of this vector in terms of the real-valued matrices.
%         [Rx1, Rx2, Rx3] = vectortomatrices(Yvec, N);

%         % And then get the representation in terms of complex-valued matrices
%         % Substituting in the real valued matrix elements
%         x1 = 1/2*(Rx1 + Rx1.') + 1/(2*sqrt(-1))*(Rx1 - Rx1.');
%         x2 = 1/2*(Rx2 + Rx2.') + 1/(2*sqrt(-1))*(Rx2 - Rx2.');
%         x3 = 1/2*(Rx3 + Rx3.') + 1/(2*sqrt(-1))*(Rx3 - Rx3.');

%         % Act with the scalar field map that takes these three matrices into a single matrix
%         P = (1/2) * (L1 * x1 + x1 * L1 + L2 * x2 + x2 * L2 + L3 * x3 + x3 * L3);

%         % Now P is in the space of complex, hermitian matrices.

%         [field, theta, phi] = matrixtofield(P, numFieldGrid);

%         % Redefinition of theta, to conform with MATLAB's conventions for spherical coordinates.
%         % theta = linspace(pi/2, -pi/2,   numFieldGrid);

%         % Plotting a 3D figure that shows explicitly where the field is supported on the sphere, and the positioning of the spherical cap.
%         f = figure('visible','off');

%         % Returns the coordinates of a sphere in cartesian
%         [X,Y,Z]= sphere(length(field)-1);       
        
%         if(max(imag(field), [], 'all') > 1e-10)
%             fprintf('Warning, removing non-negligible imaginary component! Magnitude: %.10f', max(imag(field), [], 'all'))
%         end

%         % Plots the field values upon the sphere. 
%         s = surf(X,Y,Z, flip(real(field),1), 'FaceAlpha','0.5');
%         s.EdgeColor = 'none';
%         colorbar

%         % hold on;
%         % % Plotting a ring corresponding to the cap edge.
%         % capAngle = acos(height) - pi/2;
%         % [X,Y,Z] = sph2cart(phi(2,:), capAngle*ones(1,numFieldGrid), ones(1, numFieldGrid));
%         % plot3(X,Y,Z, '-k', 'LineWidth', 2.5)
%         % hold off;

%         % Saves the figure.
%         saveas(f, strcat('figures/fieldSupport3D_', ...
%         num2str(N), '_', num2str(capAngle), '.png'))
%     end
% end


% % [delta, Z, kArr] = constrainedminimization(N, 100, true, false, true, 'none')
% % Plots a figure of delta v. Z
% f = figure('visible','off');
% set(gca,'FontSize',22) 
% hold on

% for i = 1:length(criticalMuValues)

%     sz = 25;
%     scatter(Zbar(i,:), delta(i,:), sz, 'filled');

% end


% hold off

% xlabel('$\bar{Z}$','interpreter','latex')
% ylabel('$\Delta$','interpreter','latex')

% % % legend('minimization', 'coherent states')

% saveas(f, strcat('figures/minEigvec_', num2str(N), '.png'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Generate A
% 
% % Iterate through Z

% for i = 1:length(criticalMuValues)

%     for j = 1:numSubspaceVectors

%         height = Zbar(i,j);

%         Y = subspaceVectors(:,i,j)'*A;

%         [x1, x2, x3] = vectortomatrices(Y, N);

%         P = (1/2) * (L1 * x1 + x1 * L1 + L2 * x2 + x2 * L2 + L3 * x3 + x3 * L3);

%         [field, theta, phi] = matrixtofield(P, numFieldGrid);

%         % Redefinition of theta, to conform with MATLAB's conventions for spherical coordinates.
%         theta = linspace(pi/2, -pi/2,   numFieldGrid);

%         % Plotting a 3D figure that shows explicitly where the field is supported on the sphere, and the positioning of the spherical cap.
%         f = figure('visible','off');

%         % Returns the coordinates of a sphere in cartesian
%         [X,Y,Z]= sphere(length(field)-1);        

%         % Plots the field values upon the sphere. We need to flip the field matrix along the 1st axis so that the values correspond correctly to MATLAB's spherical coordinates conventions.
%         s = surf(X,Y,Z, abs(flip(field,1)), 'FaceAlpha','0.5');
%         s.EdgeColor = 'none';
%         colorbar

%         hold on;
%         % Plotting a ring corresponding to the cap edge.
%         capAngle = acos(height/(N-1)) - pi/2;
%         [X,Y,Z] = sph2cart(phi, capAngle*ones(1,numFieldGrid), ones(1,length(phi)));
%         plot3(X,Y,Z, '-k', 'LineWidth', 2.5)
%         hold off;
        
%         % Saves the figure.
%         saveas(f, strcat('figures/fieldSupport3D_', ...
%         num2str(N), '_', num2str(capAngle), '.png'))


%     end

% end

% For each subspace vector, represent it in the 3N^2 dimensional space

% Represent the 3N^2 dimensional vectors as 3 N x N matrices

% Compute \phi from these three matrices

% Map the N x N matrix \phi to function space

% Plot it!



% % Iterates over cap defining inclination angle
% for t = 1:length(angleArr)

%         capAngle = angleArr(t);
%         capHeight = (N-1) * cos(-capAngle + pi/2);
        
%         fprintf('\nangle = %.2f (%d of %d)\n', capAngle, t, length(angleArr));

%         % Returns a field function upon the non-commutative sphere that is in correspondence with the given projection operator.
%         fprintf('Mapping projection operator to a field.\n')
%         [field, theta, phi] = matrixtofield(P, numFieldGrid);

%         fprintf('Producing figures.\n')

        % % Redefinition of theta, to conform with MATLAB's conventions for spherical coordinates.
        % theta = linspace(pi/2, -pi/2,   numFieldGrid);

        % % Plotting a 3D figure that shows explicitly where the field is supported on the sphere, and the positioning of the spherical cap.
        % f = figure('visible','off');

        % % Returns the coordinates of a sphere in cartesian
        % [X,Y,Z]= sphere(length(field)-1);        

        % % Plots the field values upon the sphere. We need to flip the field matrix along the 1st axis so that the values correspond correctly to MATLAB's spherical coordinates conventions.
        % s = surf(X,Y,Z, abs(flip(field,1)), 'FaceAlpha','0.5');
        % s.EdgeColor = 'none';
        % colorbar

        % hold on;
        % % Plotting a ring corresponding to the cap edge.
        % [X,Y,Z] = sph2cart(phi, capAngle*ones(1,numFieldGrid), ones(1,length(phi)));
        % plot3(X,Y,Z, '-k', 'LineWidth', 2.5)
        % hold off;
        
        % % Saves the figure.
        % saveas(f, strcat('figures/fieldSupport3D_', ...
        % num2str(N), '_', num2str(capAngle), '.png'))
        
        % % Plots a figure that shows the same data but along only the theta coordinate. 
        % f = figure('visible','off');
        % % For each lattitudinal slice, we sum over absolute field values at all azimuthal values. 
        % v = sum(abs(field), 2);
        % v = v/norm(v);
        % plot(theta, v);
        % % Plotting a line indicating the edge of the cap.
        % xline(capAngle,'-',{'\theta^*'});
        % saveas(f, strcat('figures/fieldSupport_Theta_', ...
        % num2str(N), '_', num2str(capAngle), '.png'))


%     end

%     % Plots a figure of delta v. Z for this N
%     f = figure('visible','off');
%     % For each lattitudinal slice, we sum over absolute field values at all azimuthal values. 
%     hold on
%     for k = 1:size(Zbar,1)

%         scatter(Zbar(k,:), delta(k,:));

%     end
%     hold off
%     % Plotting a line indicating the edge of the cap.
%     saveas(f, strcat('figures/minEigvecs_', num2str(N),'.png'))

% end


% Project the height operator appropriately

% Find {|k>}

% iterating over height:
    % take <k|A to get a state in the full 3N^2 space
    % unpack it into matrices
    % take the scalar field projection
    % map this N x N matrix to function space
    % Plot it!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Debugging sphharmat.m
% N = 3;
% [A, xi, xi2] = compressedheightoperators(N)

% [V,D] = eig(xi);
% vec = V(:,1);


% numFieldGrid = 25;
% [field, theta, phi, c, expansionTest, Y] = matrixtofield(X, numFieldGrid);
% % M = max(f,[],'all');
% % disp(M)

% % Redefinition of theta, to conform with MATLAB's conventions for spherical coordinates.
% theta = linspace(pi/2, -pi/2,   numFieldGrid);

% % Plotting a 3D figure that shows explicitly where the field is supported on the sphere, and the positioning of the spherical cap.
% f = figure('visible','off');

% % Returns the coordinates of a sphere in cartesian
% [X,Y,Z]= sphere(length(field)-1);        

% % Plots the field values upon the sphere. We need to flip the field matrix along the 1st axis so that the values correspond correctly to MATLAB's spherical coordinates conventions.
% s = surf(X,Y,Z, abs(flip(field,1)), 'FaceAlpha','0.5');
% s.EdgeColor = 'none';
% colorbar

% % hold on;
% % % Plotting a ring corresponding to the cap edge.
% % [X,Y,Z] = sph2cart(phi, capAngle*ones(1,numFieldGrid), ones(1,length(phi)));
% % plot3(X,Y,Z, '-k', 'LineWidth', 2.5)
% % hold off;

% % Saves the figure.
% saveas(f, strcat('figures/fieldSupport3D_', ...
% num2str(N), '.png'))

% % Plots a figure that shows the same data but along only the theta coordinate. 
% f = figure('visible','off');
% % For each lattitudinal slice, we sum over absolute field values at all azimuthal values. 
% v = sum(abs(field), 1);
% v = v/norm(v);
% plot(phi, v);
% % Plotting a line indicating the edge of the cap.
% % xline(capAngle,'-',{'\theta^*'});
% saveas(f, strcat('figures/fieldSupport_Theta_', ...
% num2str(N), '_', '.png'))

% Run constrainedminimization
% N = 5;
% numAngles = 25;
% [delta_cm, Z_cm, kArr] = constrainedminimization(N, 500, false, false, false, 'nullK');
% [delta, Z] = coherentstates(N, numAngles, 100, false, false, true, 'nullk');

% % Plots a figure of delta v. Z
% f = figure('visible','off');
% set(gca,'FontSize',22) 
% hold on
% scatter(Z_cm, delta_cm, sz, 'filled');
% for i = 1:numAngles

%     sz = 25;
%     scatter(Z(i,:), delta(i,:), sz, 'black', 'filled');

% end


% hold off

% xlabel('$\bar{Z}$','interpreter','latex')
% ylabel('$\Delta$','interpreter','latex')

% legend('minimization', 'coherent states')

% saveas(f, strcat('figures/coherentstates_', num2str(N), '.png'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Debugging matrixtofield. There are serious problems with the spherical harmonic matrices being constructed. 

% N = 7;
% capAngle = 0;


% numMuValues        = 1000;
% numSubspaceVectors = 10;
% numFieldGrid       = 20;

% % Returns compressions of the Z and Z^2 operators into the subspace associated with the scalar field, within the non-null sector.
% fprintf('Generating compressed height operators.\n')
% [xi, xi2] = compressedheightoperators(N);

% % Returns Z, delta, and the vectors themselves from a search routine upon the compressed heigh operators.
% fprintf('Generating minimal eigenvectors.\n')
% [Zbar, delta, subspaceVectors, criticalMuValues] = musearch(...
%                         N, xi, xi2, numMuValues, numSubspaceVectors);

% % Returns a projection operator constructed by integrating several projection operators, each associated with a minimal eigenvector, over a range in height corresponding to a spherical cap
% fprintf('Generating integrated projection operator.\n')
% P = integratedprojection(N, capAngle, Zbar, subspaceVectors);

% % Returns a field function upon the non-commutative sphere that is in correspondence with the given projection operator.
% fprintf('Mapping projection operator to a field.\n')
% [field, theta, phi, c, expansionTest, Y] = matrixtofield(P, numFieldGrid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Here I want to see how the rank of the integrated projection operator is changing, as a function of the cutoff angle.

% N = 8;

% numMuValues        = 1000;
% numSubspaceVectors = 10;
% numFieldGrid       = 20;

% % Returns compressions of the Z and Z^2 operators into the subspace associated with the scalar field, within the non-null sector.
% fprintf('Generating compressed height operators.\n')
% [xi, xi2] = compressedheightoperators(N);

% % Returns Z, delta, and the vectors themselves from a search routine upon the compressed heigh operators.
% fprintf('Generating minimal eigenvectors.\n')
% [Zbar, delta, subspaceVectors, criticalMuValues] = musearch(...
%                         N, xi, xi2, numMuValues, numSubspaceVectors);

% angleArr = linspace(-pi/2, pi/2, 50);
% angleArr = linspace(0, 0, 1);

% % Iterates over cap defining inclination angle
% for t = 1:length(angleArr)

%     capAngle = angleArr(t); 

%     % Returns a projection operator constructed by integrating several projection operators, each associated with a minimal eigenvector, over a range in height corresponding to a spherical cap
%     fprintf('Generating integrated projection operator.\n')
%     P = integratedprojection(N, capAngle, Zbar, subspaceVectors);

%     fprintf('\nFor N = %d, at capAngle = %.2f, dimP = %d, rank(P) = %d\n', N, capAngle, length(P), rank(P))

% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Here we'll use the faster approach to generating minimal eigenvectors (critial mu) values, and we'll produce 3D plots that show how the mapped functions are supported on the sphere for various cap angles.
% NArr = 5:1:8;

% % angleArr = linspace(0, 0, 1);
% angleArr = linspace(3*pi/8, -3*pi/8, 5);

% numMuValues        = 1000;
% numSubspaceVectors = 10;
% numFieldGrid       = 20;

% % Iterates over N
% for n = 1:length(NArr)

%     N = NArr(n)
%     fprintf('\nN = %d (%d of %d)\n', N, n, length(NArr));


%     % Returns compressions of the Z and Z^2 operators into the subspace associated with the scalar field, within the non-null sector.
%     fprintf('Generating compressed height operators.\n')
%     [xi, xi2] = compressedheightoperators(N);

%     % Returns Z, delta, and the vectors themselves from a search routine upon the compressed heigh operators.
%     fprintf('Generating minimal eigenvectors.\n')
%     [Zbar, delta, subspaceVectors, criticalMuValues] = musearch(...
%                             N, xi, xi2, numMuValues, numSubspaceVectors)

%     % Iterates over cap defining inclination angle
%     for t = 1:length(angleArr)

%         capAngle = angleArr(t);

%         fprintf('\nangle = %.2f (%d of %d)\n', capAngle, t, length(angleArr));
        
%         % Returns a projection operator constructed by integrating several projection operators, each associated with a minimal eigenvector, over a range in height corresponding to a spherical cap
%         fprintf('Generating integrated projection operator.\n')
%         P = integratedprojection(N, capAngle, Zbar, subspaceVectors);

%         % Returns a field function upon the non-commutative sphere that is in correspondence with the given projection operator.
%         fprintf('Mapping projection operator to a field.\n')
%         [field, theta, phi] = matrixtofield(P, numFieldGrid);

%         fprintf('Producing figures.\n')

        % % Redefinition of theta, to conform with MATLAB's conventions for spherical coordinates.
        % theta = linspace(pi/2, -pi/2,   numFieldGrid);

        % % Plotting a 3D figure that shows explicitly where the field is supported on the sphere, and the positioning of the spherical cap.
        % f = figure('visible','off');

        % % Returns the coordinates of a sphere in cartesian
        % [X,Y,Z]= sphere(length(field)-1);        

        % % Plots the field values upon the sphere. We need to flip the field matrix along the 1st axis so that the values correspond correctly to MATLAB's spherical coordinates conventions.
        % s = surf(X,Y,Z, abs(flip(field,1)), 'FaceAlpha','0.5');
        % s.EdgeColor = 'none';
        % colorbar

        % hold on;
        % % Plotting a ring corresponding to the cap edge.
        % [X,Y,Z] = sph2cart(phi, capAngle*ones(1,numFieldGrid), ones(1,length(phi)));
        % plot3(X,Y,Z, '-k', 'LineWidth', 2.5)
        % hold off;
        
        % % Saves the figure.
        % saveas(f, strcat('figures/fieldSupport3D_', ...
        % num2str(N), '_', num2str(capAngle), '.png'))
        
        % % Plots a figure that shows the same data but along only the theta coordinate. 
        % f = figure('visible','off');
        % % For each lattitudinal slice, we sum over absolute field values at all azimuthal values. 
        % v = sum(abs(field), 2);
        % v = v/norm(v);
        % plot(theta, v);
        % % Plotting a line indicating the edge of the cap.
        % xline(capAngle,'-',{'\theta^*'});
        % saveas(f, strcat('figures/fieldSupport_Theta_', ...
        % num2str(N), '_', num2str(capAngle), '.png'))


%     end

%     % Plots a figure of delta v. Z for this N
%     f = figure('visible','off');
%     % For each lattitudinal slice, we sum over absolute field values at all azimuthal values. 
%     hold on
%     for k = 1:size(Zbar,1)

%         scatter(Zbar(k,:), delta(k,:));

%     end
%     hold off
%     % Plotting a line indicating the edge of the cap.
%     saveas(f, strcat('figures/minEigvecs_', num2str(N),'.png'))

% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inspecting the projected height operator in this new space more closely. In this block I have worked with the new N^2 dimensional space. Many of the features of the NullK space are replicated here, and I've developed an analogous method for quickly finding null eigenvectors.
%{
 % Setting N
N = 5;
angle = pi/2;

% Load, diagonalize K. Retrieve vectors spanning [Null(K)]^T
M = load(['savedcouplingmatrices/',num2str(N), '.mat']);
[V, D] = eig(real(M.K));
V = V(:, abs(diag(D)) > 1e-4);

% Constructing a projection operator onto the [Null(K]^T)] space.
T = zeros(size(V,1));
for i = 1:size(V,2)
    
    T = T + V(:,i)*V(:,i)';

end

% Constructs a matrix with orthonormal columns spanning the range of T.
T = orth(T);

% Returns the operator that projects into the space associated to the scalar field. Only runs if this data has not been computed (and saved) previously.
filename = ['savedscalarfieldprojection/',num2str(N), '.mat'];
if ~isfile(filename)
    A = real(scalarfieldprojection(N));
else
    M = load(filename);
    A = real(M.A);
end

% Compresses the operator A (scalar field), and the height operator (Z), into the subspace defined by T (non-gauge d.o.f.)
AinT = orth(T'*A*T);
xi  = T' *  heightoperator(N)    * T;
xi2 = T' * (heightoperator(N))^2 * T;

% Compresses the height operator once more, using AinT as just defined.
xi  = AinT' * xi     * AinT;
xi2 = AinT' * xi2    * AinT;

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

% Iterating over mu, we compute the spectrum of xi2 - mu*xi
mu                  = linspace(-0.87*2*(N-1), +0.87*2*(N-1), 100);
eigenvalues         = zeros(length(xi), length(mu));
lowestEigvalue      = zeros(1, length(mu));
lowestEigvaluesVar  = zeros(1, length(mu));
tolerance           = 0.01;

for i = 1:length(mu)
    eigenvalues(:,i) = eig(xi2 - mu(i)*xi);
    % retrieve the lowest eigenvalue
    lowestEigvalue(i) = eigenvalues(1);
    lowestEigvaluesVar(i) = var(eigenvalues(1:3,i));
end

TF = islocalmin(lowestEigvaluesVar);
criticalPoints = mu(TF);

% f = figure('visible','off');
% hold on
% for k = 1:length(xi)
%     scatter(mu, eigenvalues(k,:))
% end
% xline(criticalPoints(1))
% xline(criticalPoints(2))
% saveas(f,strcat('figures/inspect_', num2str(N)),'png')   

% Now we search for random vectors that lie within the degenerate subspace identified.
degenerateEigenvectors = zeros(length(xi), 2*length(criticalPoints));
subspaceVectors = zeros(numPoints, length(criticalPoints), length(xi));
numPoints = 5000;
theta  = linspace(0,pi/2, numPoints)';
points = [cos(theta), sin(theta)];

intProjOp = zeros(length(xi));
% initialize arrays to hold delta, Zbar, but separated by critical mu
delta     = zeros(numPoints, length(criticalPoints));
ZBar      = zeros(numPoints, length(criticalPoints));
for i = 1:length(criticalPoints)
    [V,D] = eig(xi2 - criticalPoints(i)*xi);
    eigenvalues = diag(D);

    % Instead of generating random points, we just generate equally
    % spaced points in the angular coordinate of the quarter circle in
    % the positive quadrant
    % retrieve the two eigenvalue indices we need.
    [~, minI] = min(eigenvalues(1:3));
    [~, maxI] = max(eigenvalues(1:3));

    degenerateEigenvectors(:, (2*i)-1:(2*i)) = [V(:, minI), V(:, maxI)];
        
        % !!! It should be possible to remove this loop by vectorizing
        for j = 1:length(points)
            
            % generate subspace vector
            subspaceVec = ...
                [V(:, minI), V(:, maxI)] * points(j,:)';
            
            % compute height and delta
            ZBar(j,i) = subspaceVec' * xi * subspaceVec;
            delta(j,i) = subspaceVec' * xi2 * subspaceVec - ...
                (subspaceVec' * xi * subspaceVec)^2;

            subspaceVectors(j, i, :) = subspaceVec;
            
            intProjOp = intProjOp + subspaceVec * subspaceVec';

        end

end
intProjOp = intProjOp / (mean(mean(intProjOp)));

% % % Runs constrained optimization and plots the result.
% [deltaMin, Z, kArr] = constrainedminimization(N, xi, xi2, -0.87*(N-1), +0.87*(N-1), 100);


f = figure('visible','off');
hold on
for k = 1:length(criticalPoints)
    scatter(ZBar(:,k), delta(:,k))
end
% scatter(Z(:), deltaMin(:), 50, 'black')
saveas(f,strcat('figures/inspect_', num2str(N)),'fig')      

% for k = 1:size(degenerateEigenvectors,2)
%     for kk = 1:size(degenerateEigenvectors,2)
%         fprintf('overlap of %d, and %d: %.10f\n', k, kk, dot(degenerateEigenvectors(:, k), degenerateEigenvectors(:,kk)))
%     end
% end

% for j = 1:length(xi)
%     plot(mu, eigenvalues(j,:))
% end
% saveas(f,strcat('figures/eigenvalues_', num2str(N)),'png')      

% % Runs constrained optimization and plots the result.
% [delta, Z, kArr] = constrainedminimization(N, xi, xi2, -0.87*(N-1), +0.87*(N-1), 600);

% % Generates a figure plotting height vs. spread of the eigenvectors found.
% f = figure('visible','off');
% scatter(Z(:), delta(:))
% saveas(f,strcat('figures/minEigvecs_AinT_', num2str(N)),'png')    

% We want to know if these are generated from the same degeneracy points as we saw previously. We can quickly get at that question by just trying to compute a delta v. Z plot with that technique. 

% Next we want to inspect the orthogonality, of the degenerate eigenvectors used to construct these minimal eigenvectors. 
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Iterating over N and cap angle with the integrated projection procedure.
%{
% Defines arrays to hold ranges in N and in cap angle
NArr = 5:1:5;
angleArr = linspace(0.01, pi-0.01, 10);

% Iterates over N
for n = 1:length(NArr)

    % Iterates over cap defining inclination angle
    % for t = 1:length(angleArr)

        % Calls the integratedprojection.m function for this value of N, and this cap angle
        fprintf('N = %d (%d of %d), angle = %.2f (%d of %d)\n', NArr(n), n, length(NArr), angleArr(t), t, length(angleArr));
        integratedprojection(NArr(n), angleArr(t));
        N = NArr(n)

        M = load(['savedcouplingmatrices/',num2str(N), '.mat']);
        [V, D] = eig(real(M.K));
        V = V(:, abs(diag(D)) > 1e-4);

        % Constructing a projection operator onto the [Null(K]^T)] space.
        T = zeros(size(V,1));
        for i = 1:size(V,2)
            
            T = T + V(:,i)*V(:,i)';

        end
        
        % Constructs a matrix with orthonormal columns spanning the range of T.
        T = orth(T);

        % Returns the operator that projects into the space associated to the scalar field. Only runs if this data has not been computed (and saved) previously.
        filename = ['savedscalarfieldprojection/',num2str(N), '.mat'];
        if ~isfile(filename)
            A = scalarfieldprojection(N);
        else
            M = load(filename);
            A = M.A;
        end

        % Compresses the operator A (scalar field), and the height operator (Z), into the subspace defined by T (non-gauge d.o.f.)
        AinT = orth(real(T'*A*T));
        xi  = T' *  heightoperator(N)    * T;
        xi2 = T' * (heightoperator(N))^2 * T;
        
        % Compresses the height operator once more, using AinT as just defined.
        xi  = AinT' *  xi    * AinT;
        xi2 = AinT' * xi2 * AinT;

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

        % Running constrained minimization to find minimal eigenvectors
        [delta, Z, kArr] = constrainedminimization(N, xi, xi2, -0.87*(N-1), +0.87*(N-1), 600);
        fprintf('N = %d\n', N);

        % Generates a figure plotting height vs. spread of the eigenvectors found.
        % f = figure('visible','off');
        % scatter(Z(:), delta(:))
        % saveas(f,strcat('figures/minEigvecs_AinT_', num2str(N)),'png')      

    % end

end 
%}