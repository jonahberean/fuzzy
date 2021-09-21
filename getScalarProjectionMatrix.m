function A = getScalarProjectionMatrix(N)

filename = ['savedScalarProjectionMatrix/', num2str(N), '.mat'];
if isfile(filename)
    M                 = load(filename);
    A              = M.A;
else

    [J1, J2, J3, Jminus, Jplus] = computeSu2(N);

    % rescaling the su(2) generators so as to normalize the fuzzy sphere radius to 1.
    nu     = 2 / sqrt(N^2 - 1);
    J1 = nu*J1;
    J2 = nu*J2;
    J3 = nu*J3;

    % Defines three matrices filled with real, symbolic variables
    syms Rx1 [N N];
    syms Rx2 [N N];
    syms Rx3 [N N];

    % Constructs a symbolic column vector from the real matrix degrees of freedom.
    Y = computeVectorFromMatrices(Rx1, Rx2, Rx3);

    % Defines complex hermitian matrices in terms of the real matrices
    x1 = computeComplexMatrixFromReal(Rx1);
    x2 = computeComplexMatrixFromReal(Rx2);
    x3 = computeComplexMatrixFromReal(Rx3);

    % Computes phi, the scalar field matrix, or radial component of the matrices,
    % using the complex, hermitian matrices.
    phi = (1/2) * (J1 * x1 + x1 * J1 + J2 * x2 + x2 * J2 + J3 * x3 + x3 * J3);

    % We want to define a linear map between real vector spaces so we map the
    % complex, hermitian phi to a real matrix
    Rphi = computeRealMatrixFromComplex(phi);

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
    save(filename,'A');
end