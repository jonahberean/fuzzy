function [f, theta, phi, c] = matrixtofield(A, numFieldGrid)

%{
matrixtofield.m - maps a matrix into a function using a basis of
spherical harmonics

Inputs:
   A         - matrix
   numPoints - dimension of the matrix grid to be returned

Outputs:
   f     - matrix grid of f(theta, phi) values
   theta - vector of theta values
   phi   - vector of phi values
   c     - matrix of expansion coefficients, c_{jm}

Other m-files required: generatesphericalharmonicmatrices.m, harmonicY.m
Subfunctions:
MAT-files required: none

Author: Jonah Berean-Dutcher
email: jbd@phas.ubc.ca
June 2021; Last revision: 29-July-2021 
%}

% 
% A = P;
% The dimension of the given matrix is the dimension we are working in
N = length(A);

% initialize c, vector of the coefficients of
% the matrix (and counterpart spherical harmonic) expansion.
c = zeros(N-1, (2 * (N-1) + 1));

% At this N, construct the basis matroices \hat{Y}. 
filename = ['savedsphericalharmonicmatrices/',num2str(N), '.mat'];
if ~isfile(filename)
    fprintf('Generating spherical harmonic matrices.\n');
    [Y, Jplus, Jminus, Ji] = sphharmat(N);
else
    M = load(filename);
    Y = M.Y;
end

% define vectors of theta and phi for generating values of f(theta, phi)
% Note that the convention of taking theta values from 0 to pi is NOT the same as MATLAB's convention for working with spherical data (pi/2 to -pi/2)
[x,y,z]= sphere(numFieldGrid-1);
[phi,theta,r] = cart2sph(x,y,z);
theta = theta + pi/2;

% initialize f(theta,phi)
f = zeros(length(theta), length(phi));

% initialize a matrix which will hold the expansion of A in terms of the
% Y_{j,m} matrices, for later testing of the comparison
expansionTest = zeros(size(A));

% Itertatively take the inner product of A with each of these matrices, so
% as to compute the coefficients of A's expansion in this basis.
for j = 0:(N-1)    

    % fprintf('Iterating over j. Currently at j = %.2f (out of %.2f).\n', j, jmax);
    
    mValues = -j:1:j;
    
    for m = 1:length(mValues)
                              
        % computing coefficients and constructing values of function space
        % correspondence
        c(j+1,m) = (1/N) * trace(Y(:,:,j+1,m)' * A);             
        
        % constructing the expansion of A for later testing
        expansionTest = expansionTest + c(j+1,m) * Y(:,:,j+1,m);
        
    end
    
end

% Here we normalize the vector of coefficients.
normFactor = norm(c(:));
c = reshape(c(:)/normFactor, [N,2*N-1]);


% Testing that the expansion of A is equal to the original A.
if (max(abs(expansionTest - A), [], 'all') > 0.001)
    
    statement = '\nWarning: The constructed expansion of A in terms the Y_{j,m} deviates from the original A by %.2f \n';
    fprintf(statement, max(abs(expansionTest - A), [], 'all'))
    
end

% Computing the functioin values
for j = 0:(N-1)    

    % fprintf('Iterating over j. Currently at j = %.2f (out of %.2f).\n', j, jmax);
    
    mValues = -j:1:j;
    
    for m = 1:length(mValues)
        
        
        % Note: The function harmonicY does not evaluate a grid all at once, it evaluates f(theta, phi) once for each value of theta and of phi, as ordered pairs running down each column in tandem.

        for p = 1:length(phi)

            % Adds a column of values into f, that's at every theta value, and at one phi value. 
            % f(:,p) = f(:,p) + c(j+1,m) * ...
            %     harmonicY(j, mValues(m), ...
            %     theta, ones(1,length(theta))*phi(p), 'type', 'real')';

            f(:,p) = f(:,p) + c(j+1,m) * harmonicY(j, mValues(m), theta(:,p), phi(:,p));
        end
        
    end
    
end

% Removes negligible imaginary component from field
if(max(imag(f), [], 'all') > 1e-10)
    fprintf('Warning, removing non-negligible imaginary component! Magnitude: %.10f', max(imag(f), [], 'all'))
end
f = real(f);   


