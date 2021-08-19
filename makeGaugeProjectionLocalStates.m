
N = 5;
overwriteFlag = true;
filelabel = 'null';
K = getCouplingMatrix(N);
[V, D] = eig(K);
kernel = V(:, abs(diag(D)) < 1e-4);
Z = getHeightOperator(N);
xi = kernel' * Z * kernel;
xi2 = xi^2;

% Remove the negligible imaginary numerical components
xi = real(xi);
xi2 = real(xi2);

% Rounds the xi and xi2 matrices until they are hermitian
roundToDigit = 18;
tf = false;
while (~tf)
    xi = round(xi, roundToDigit);
    xi2 = round(xi2, roundToDigit);
    tf = ishermitian(xi) & ishermitian(xi2);
    roundToDigit = roundToDigit - 1;
end   
% % xi = Z;
% % xi2 = Z^2;
% heightGridSize = 50;
[varZ, Z, localStateVectors] = getLocalStatesConOpt(N, xi, xi2, heightGridSize, true, 'test');
% Generates a table of random coefficients


% [V,D] = eig(xi);
% varZ = zeros(1,size(V,2));
% Z = zeros(1,size(V,2));
% y = -1 + 2*randn(100,size(V,2));
% listOfCoefficients = bsxfun(@rdivide,y,sqrt(sum(y.^2,2)));
% for i = 1:100
%     kOpt = V*listOfCoefficients(i,:)';
%     varZ(i) = kOpt' * xi2 * kOpt - (kOpt' * xi * kOpt)^2;
%     Z(i)  = kOpt' * xi  * kOpt;

% end
% Plots a figure of delta v. Z
f = figure('visible','off');

sz = 25;
scatter(Z, varZ, sz, 'filled');

xlabel('$\langle Z \rangle$','interpreter','latex')
ylabel('$\Delta_Z$','interpreter','latex')

ax = gca;
set(gca,'FontSize',22)

saveas(f, strcat('savedLocalStatesConOpt/', num2str(N), '_heightGridSize_', num2str(heightGridSize), '_null','.png'))
saveas(f, strcat('savedLocalStatesConOpt/', num2str(N), '_heightGridSize_', num2str(heightGridSize), '_null','.fig'))
