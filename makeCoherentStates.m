% Defines a figure window for plotting, sets font and marker size
f = figure('visible','off');
ax = gca;
set(gca,'FontSize',22)
sz = 200;
hold on

% Defines N
N = 5;

% Retrieves coupling matrix
K = getCouplingMatrix(N);

% Settings for getLocalStatesCoherent()
% overwriteFlag = true;
overwriteFlag = false;
numVectorsPerHeight = 20;
heightGridSize = 20;
filelabelFull = 'full';
filelabelKernelComplement = 'kernelComplement';

% Retrieves the compressions, and coherent state results for the full space
Z = computeHeightOperator(N);
projectionMatrix = eye(size(K));
[varZ, Z, localStateVectors] = ...
    getLocalStatesCoherent(N, projectionMatrix, Z, Z^2, ...
    heightGridSize, 1, overwriteFlag, filelabelFull);

% Plots Z vs. varZ
legendLabel = ['$Y_{|\hat{n}\rangle \langle \hat{n}|}$'];
scatter(Z(:), varZ(:), sz, '.', 'DisplayName', legendLabel);

% Retrieves the coherent state results for the kernelComplement space
[xi, xi2, O] = computeCompression(N, K, 'kernelComplement');
projectionMatrix = O;
[varZ, Z, localStateVectors] = getLocalStatesCoherent(N, projectionMatrix, ...
    xi, xi2, heightGridSize, ...
    numVectorsPerHeight, overwriteFlag, filelabelKernelComplement);

% Plots Z vs. varZ
legendLabel = ['$Y_{|\hat{n}\rangle \langle \hat{n}|} \in W$'];
scatter(Z(:), varZ(:), sz, '.', 'DisplayName', legendLabel);


hold off
xlabel('$\bar{Z}$','interpreter','latex')
ylabel('$\langle (\Delta Z)^2 \rangle$','interpreter','latex')
legend('interpreter','latex')

% Prints a vectorized eps image of the figure
print(f,'LocalCoherent.eps','-depsc');