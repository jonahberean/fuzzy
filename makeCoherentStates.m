%
% generate varZ v. Z for coherent states in the full theory
N = 10;
overwriteFlag = false;
K = getCouplingMatrix(N);
[V, D] = eig(K);
kernel = V(:, abs(diag(D)) > 1e-4);
Z = getHeightOperator(N);
xi = kernel' * Z * kernel;
xi2 = kernel' * Z^2 * kernel;
heightGridSize = 5;

numVectorsPerHeight = 10;

[varZ, Z, localStateVectors] = getLocalStatesCoherent(N, kernel, xi, xi2, heightGridSize, numVectorsPerHeight, false, 'null');

% Plots a figure of delta v. Z
f = figure('visible','off');

sz = 25;
hold on
scatter(Z(:), varZ(:), sz, 'filled');

[varZ, Z, localStateVectors] = getLocalStatesCoherent(N, eye(3*N^2), getHeightOperator(N), getHeightOperator(N)^2, 20, 2, false, 'full');

scatter(Z(:), varZ(:), sz, 'filled');

[varZ, Z, localStateVectors] = getLocalStatesConOpt(N, xi, xi2, 50, false, 'null')

scatter(Z(:), varZ(:), sz, 'filled');

[varZ, Z, localStateVectors] = getLocalStatesConOpt(N, getHeightOperator(N), getHeightOperator(N)^2, 50, false, 'full')

scatter(Z(:), varZ(:), sz, 'filled');

hold off

xlabel('$\langle Z \rangle$','interpreter','latex')
ylabel('$\Delta_Z$','interpreter','latex')
legend('$|\Psi_{\hat{n}}\rangle \in \mathrm{Null}(K)$', '$|\Psi_{\hat{n}}\rangle \in {R}^{3N^2}$', '$|\Psi_{\mathrm{min}}\rangle\in \mathrm{Null}(K)$','$|\Psi_{\mathrm{min}}\rangle \in R^{3N^2}$', 'interpreter','latex', 'location', 'east')
ax = gca;
set(gca,'FontSize',22)

saveas(f, strcat('savedLocalStatesCoherent/', num2str(N), '_heightGridSize_', num2str(heightGridSize), '_heightGridSize_', num2str(numVectorsPerHeight),'.png'))
saveas(f, strcat('savedLocalStatesCoherent/', num2str(N), '_heightGridSize_', num2str(heightGridSize), '_heightGridSize_', num2str(numVectorsPerHeight),'.fig'))