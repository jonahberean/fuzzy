
%
N = 3;
[J1, J2, J3,Jminus,Jplus] = getSu2(N);
A = getScalarProjectionMatrix(N, J1, J2, J3);

% Compress the height operators into the scalar field subspace
xi  = A * getHeightOperator(N)     * A';
xi2 = A * getHeightOperator(N)^2    * A';

% Test to make sure I'm not rounding more than I think
if (imag(xi2) > 1e-10) | (imag(xi) > 1e-10)

    fprintf('Warning, rounding off non-negligible imaginary components!')

end

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

gridSize = 20;
[varZ, Z, localStateVectors] = getLocalStatesConOpt(N, xi, xi2, gridSize);

% Plots a figure of delta v. Z
f = figure('visible','off');

sz = 25;
scatter(Z, varZ, sz, 'filled');

xlabel('$\langle Z \rangle$','interpreter','latex')
ylabel('$\Delta_Z$','interpreter','latex')

ax = gca;
set(gca,'FontSize',22)

saveas(f, strcat('savedLocalStatesConOpt/', num2str(N), '_gridSize_', num2str(gridSize), '.png'))
saveas(f, strcat('savedLocalStatesConOpt/', num2str(N), '_gridSize_', num2str(gridSize), '.fig'))

