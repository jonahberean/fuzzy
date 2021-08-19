function A = scalarfieldprojection(N)

%region - doc
%{
% scalarfieldprojection.m - Constructs a projection matrix onto the vector space
of d.o.f. associated with a scalar field

Syntax: A = scalarfieldprojection(N)

Inputs:
   N - coordinate matrix dimension

Outputs:
   A - projection matrix
          
Other m-files required: su2generators.m, matricestovector.m
Subfunctions: none
MAT-files required: none

Author: Jonah Berean-Dutcher
email: jbd@phas.ubc.ca
June 2021; Last revision: 18-Aug-2021
%}
%endregion - doc

fprintf('Generating scalar field projection.\n')

% Returns generators of the su(2) Lie algebra
[L1, L2, L3, ~, ~] = su2generators(N);

% rescaling the su(2) generators so as to normalize the fuzzy sphere radius to 1.
nu     = 2 / sqrt(N^2 - 1);
L1 = nu*L1;
L2 = nu*L2;
L3 = nu*L3;

% Defines three matrices filled with real, symbolic variables
syms Rx1 [N N];
syms Rx2 [N N];
syms Rx3 [N N];

% Constructs a symbolic column vector from the real matrix degrees of freedom.
Y = matricestovector(Rx1, Rx2, Rx3, N);

% Defines complex hermitian matrices in terms of the real matrices
x1 = (1/2*(Rx1+Rx1.')-1/(2*sqrt(-1))*(Rx1-Rx1.'));
x2 = (1/2*(Rx2+Rx2.')-1/(2*sqrt(-1))*(Rx2-Rx2.'));
x3 = (1/2*(Rx3+Rx3.')-1/(2*sqrt(-1))*(Rx3-Rx3.'));

% Computes phi, the scalar field matrix, or radial component of the matrices,
% using the complex, hermitian matrices.
phi = (1/2) * (L1 * x1 + x1 * L1 + L2 * x2 + x2 * L2 + L3 * x3 + x3 * L3);

% We want to define a linear map between real vector spaces so we map the
% complex, hermitian phi to a real matrix
Rphi = 1/2*(phi+phi.')+1/(2*sqrt(-1))*(phi-phi.');

%region - constructing Y_phi
% Constructing Yphi, the column vector of N^2 length, with just the
% entries of phi, as symbolic expressions of the entries of the three
% matrices.
Yphi = sym(zeros(1,N^2));

% Append the main diagonal to Yphi
Yphi(1:N)       = diag(Rphi);

% index of the first element after the main diagonals
start = N+1;

% placeholder for counting indices
step = 0;

% Looping from 1 to N-1, the number of diagonals of a given sign i.e. 
% upper or lower
for i = 1:N
    start = start + step;
    Yphi(start           : start+(N-1-i))         = diag(Rphi, i);
    Yphi(start+1*(N-i)   : start+2*(N-i)-1)       = diag(Rphi, -i);

    step = 2*(N-i);
end
%endregion - constructing Y_phi

% Construct the matrix that carries out the Y -> Yphi transformation by
% taking partial derivatives
A = zeros(N^2, 3*N^2);
for i = 1:3*N^2

    for j = 1:N^2

        A(j,i) = diff(Yphi(j), Y(i));

    end

end


% save the results
saveFilename = ['savedscalarfieldprojection/', num2str(N), '.mat'];
save(saveFilename,'A');
