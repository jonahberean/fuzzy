function [A, xi, xi2] = compressedheightoperators(N, subspaceFlag)

%region - doc
%{
integratedprojection.m - Devloping a method for exploring geometry in the matrix mechanics that involves integrating over many projection operators that we believe to be in correspondence with states on the non-commutative sphere sufficiently localized about a point.

Inputs:
    N            - coordinate matrix dimension
    subspaceFlag - string to instruct which subspace to project into

Outputs:
    xi  - compression of the height operator, Z
    xi2 - compression of Z^2

Other m-files required: scalarfieldprojection.m, heightoperator.m
Subfunctions: none
MAT-files required: none

Author: Jonah Berean-Dutcher
email: jbd@phas.ubc.ca
July 2021; Last revision: 05-Aug-2021 
%}
%endregion - doc

fprintf('Generating compressed height operators for N = %d.\n', N)

%region - deprecated
% % Load, diagonalize K. Retrieve vectors spanning [Null(K)]^T
% M = load(['savedcouplingmatrices/',num2str(N), '.mat']);
% [V, D] = eig(real(M.K));
% V = V(:, abs(diag(D)) > 1e-4);

% % Constructing a projection operator onto the [Null(K]^T)] space.
% T = zeros(size(V,1));
% for i = 1:size(V,2)
    
%     T = T + V(:,i)*V(:,i)';

% end

% Constructs a matrix with orthonormal columns spanning the range of T.
% T = orth(T);
%endregion - deprecated

% Returns the projection matrix to be used in compressing the height operators
filename = ['savedscalarfieldprojection/',num2str(N), '.mat'];
if ~isfile(filename)
    A = scalarfieldprojection(N);
else
    M = load(filename);
    A = M.A;
end

% Compress the height operators into the scalar field subspace
xi  = A * heightoperator(N)     * A';
xi2 = A * heightoperator(N)^2    * A';

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
