
% Defines array of matrix sizes for iteration
% NArr = [5, 7, 9];
NArr = 7;

% Defines a figure window for plotting, sets font size
f = figure('visible','off');
ax = gca;
set(gca,'FontSize',22)
hold on

% Iterates over array of sizes
for i = 1:length(NArr)

    N = NArr(i);

    % Retrieves coupling matrix
    K = getCouplingMatrix(N);

    % Retrieves compressions
    filelabel = 'kernelComplement';
    % filelabel = 'full';
    % filelabel = 'scalar';
    [xi, xi2, O] = computeCompression(N, K, filelabel);

    % Retrieves the constraint minimization results 
    % overwriteFlag = true;
    overwriteFlag = true;
    heightGridSize = 5000;
    [varZ, Z, localStateVectors, filename] = ...
        getLocalStatesSearch(N, xi, xi2, heightGridSize, overwriteFlag, filelabel);

    % Plots Z vs. varZ
    legendLabel = ['$ N =', num2str(N), '$'];
    scatter(Z, varZ, 'DisplayName', legendLabel, 'LineWidth',2.5);

end

hold off
grid on
xlabel('$\bar{Z}_W$','interpreter','latex')
ylabel('$\langle (\Delta Z)^2 \rangle_W$','interpreter','latex')
% legend('interpreter','latex')

% Creates a filename for the figure to save
figureFilename = strcat('figs/LocalSearch_', filelabel, '.eps');

% Prints a vectorized eps image of the figure
print(f,figureFilename,'-depsc');


% Creates a filename for the figure to save
figureFilename = strcat('figs/LocalSearch_', filelabel, '.png');

% Prints a png image of the figure
print(f,figureFilename,'-dpng');

% % Iterating over \bar{Z}_W, computing successive dot products of adjacent vectors
% dotProduct = zeros(1,length(localStateVectors)-1);
% for i = 1:length(Z)-1

%     dotProduct(i) = dot(localStateVectors(:,i), localStateVectors(:,i+1));

% end

